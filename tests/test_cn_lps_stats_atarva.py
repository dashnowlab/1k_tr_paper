"""Scientific and integration tests for the combined ATaRVa CN/LPS table."""

from __future__ import annotations

import csv
import math
import statistics
import sys
import unittest
from collections import Counter
from pathlib import Path

from scipy.stats import binomtest


ROOT = Path(__file__).resolve().parents[1]
DOWNSTREAM = ROOT / "src" / "downstream"
sys.path.insert(0, str(DOWNSTREAM))

import cn_lps_stats_atarva as stats  # noqa: E402
from lpsSTR_atarva import longest_target_lps  # noqa: E402


OUTPUT = ROOT / "data" / "1000g_output" / "1000g_ont_test.cn_lps_extended.tab"
VCF = ROOT / "data" / "1000g_output" / "1000g_ont_test.vcf"
STAT_OUTPUT = ROOT / "data" / "1000g_output" / "1000g_ont_test.statSTR_atarva.tab"
LPS_OUTPUT = ROOT / "data" / "1000g_output" / "1000g_ont_test.lpsSTR_atarva.tab"


def as_mapping(columns: list[str], values: list[str]) -> dict[str, str]:
    return dict(zip(columns, values, strict=True))


def assert_optional_close(
    case: unittest.TestCase, observed: str, expected: str, tolerance: float
) -> None:
    if observed == "." or expected == ".":
        case.assertEqual(observed, expected)
    else:
        case.assertAlmostEqual(float(observed), float(expected), delta=tolerance)


def linear_percentile(sorted_values: list[float], percentile: float) -> float:
    """Independently implement NumPy's default Hyndman-Fan type-7 percentile."""
    if not sorted_values:
        return math.nan
    rank = (len(sorted_values) - 1) * percentile / 100.0
    lower = math.floor(rank)
    upper = math.ceil(rank)
    fraction = rank - lower
    return sorted_values[lower] + fraction * (sorted_values[upper] - sorted_values[lower])


def expand_histogram(text: str) -> list[float]:
    if text == ".":
        return []
    values: list[float] = []
    for item in text.split(","):
        allele, count = item.split(":", 1)
        allele = allele.removesuffix("x")
        values.extend([float(allele)] * int(count))
    return sorted(values)


class TestLongestPureSegment(unittest.TestCase):
    def test_missing_and_no_target_repeat_are_distinct(self) -> None:
        self.assertIsNone(longest_target_lps(".", "CAG"))
        self.assertEqual(longest_target_lps("TTTT", "CAG"), 0)

    def test_decomposed_runs_interruptions_and_literal_sequence(self) -> None:
        self.assertEqual(longest_target_lps("(CAG)4", "CAG"), 4)
        self.assertEqual(longest_target_lps("(CAG)3-T-(CAG)5", "CAG"), 5)
        self.assertEqual(longest_target_lps("CAGCAGTCAG", "CAG"), 2)
        self.assertEqual(longest_target_lps("(CAGCAG)2", "CAG"), 4)

    def test_rotation_and_reverse_complement_are_equivalent(self) -> None:
        self.assertEqual(longest_target_lps("(AGC)4", "CAG"), 4)
        self.assertEqual(longest_target_lps("(CTG)4", "CAG"), 4)


class TestDistributionStatistics(unittest.TestCase):
    def test_known_distribution(self) -> None:
        result = as_mapping(
            stats.SUMMARY_SUFFIXES,
            stats._distribution_summary([1, 1, 2, 4], rare_threshold=0.30),
        )
        expected = {
            "min": 1.0,
            "median": 1.5,
            "mean": 2.0,
            "99th_percentile": 3.94,
            "max": 4.0,
            "range": 3.0,
            "variance": 1.5,
            "stdev": math.sqrt(1.5),
            "coefficient_of_variation": math.sqrt(1.5) / 2.0,
            "iqr": 1.5,
            "mode": 1.0,
            "major_allele_frequency": 0.5,
            "effective_allele_count": 8.0 / 3.0,
            "heterozygosity": 0.625,
            "entropy": 1.5,
            "rare_allele_fraction": 0.5,
        }
        for field, value in expected.items():
            self.assertAlmostEqual(float(result[field]), value, delta=1e-6, msg=field)
        self.assertEqual(result["unique_allele_lengths"], "3")
        self.assertEqual(result["total_genotyped_alleles"], "4")
        self.assertEqual(result["allele_histogram"], "1x:2,2x:1,4x:1")
        self.assertEqual(result["singleton_allele_count"], "2")
        self.assertEqual(result["rare_allele_count"], "2")

    def test_mode_tie_uses_smallest_value(self) -> None:
        result = as_mapping(
            stats.SUMMARY_SUFFIXES,
            stats._distribution_summary([3, 3, 4, 4], rare_threshold=0.01),
        )
        self.assertEqual(result["mode"], "3")

    def test_empty_distribution(self) -> None:
        result = as_mapping(
            stats.SUMMARY_SUFFIXES,
            stats._distribution_summary([], rare_threshold=0.01),
        )
        self.assertEqual(result["total_genotyped_alleles"], "0")
        self.assertEqual(result["unique_allele_lengths"], "0")
        self.assertEqual(result["allele_histogram"], ".")
        self.assertEqual(result["mean"], ".")

    def test_coefficient_of_variation_is_undefined_at_zero_mean(self) -> None:
        result = as_mapping(
            stats.SUMMARY_SUFFIXES,
            stats._distribution_summary([0, 0], rare_threshold=0.01),
        )
        self.assertEqual(result["coefficient_of_variation"], ".")

    def test_reference_relative_statistics(self) -> None:
        result = as_mapping(
            stats.CN_REFERENCE_COLUMNS,
            stats._cn_reference_summary([8, 10, 10, 13], ref_cn=10),
        )
        expected = {
            "cn_reference_frequency": 0.5,
            "cn_non_reference_fraction": 0.5,
            "cn_mean_difference_from_reference": 0.25,
            "cn_median_difference_from_reference": 0.0,
            "cn_max_expansion_from_reference": 3.0,
            "cn_max_contraction_from_reference": 2.0,
            "cn_expanded_allele_fraction": 0.25,
            "cn_contracted_allele_fraction": 0.25,
        }
        for field, value in expected.items():
            self.assertAlmostEqual(float(result[field]), value, delta=1e-6, msg=field)

    def test_cn_lps_pair_statistics_and_partition(self) -> None:
        result = as_mapping(
            stats.COMPARISON_COLUMNS,
            stats._comparison_summary([(5, 5), (6, 4), (4, 5), (0, 0)]),
        )
        self.assertEqual(result["cn_lps_paired_alleles"], "4")
        self.assertEqual(result["fully_pure_allele_count"], "2")
        self.assertEqual(result["interrupted_allele_count"], "1")
        self.assertEqual(result["lps_exceeds_cn_count"], "1")
        self.assertAlmostEqual(float(result["cn_lps_mean_difference"]), 0.25, delta=1e-6)
        self.assertAlmostEqual(
            float(result["cn_lps_mean_ratio"]), (1 + 4 / 6 + 5 / 4) / 3, delta=1e-6
        )


class TestHweCompatibility(unittest.TestCase):
    def test_two_sided_binomial_matches_scipy(self) -> None:
        cases = [(0, 10, 0.2), (3, 10, 0.2), (7, 12, 0.65), (10, 10, 0.8)]
        for successes, trials, probability in cases:
            with self.subTest(successes=successes, trials=trials, probability=probability):
                expected = binomtest(successes, trials, probability).pvalue
                self.assertAlmostEqual(
                    stats._binomial_two_sided(successes, trials, probability),
                    expected,
                    delta=1e-12,
                )

    def test_hwe_uses_expected_and_observed_heterozygosity(self) -> None:
        genotypes = [(1, 1), (1, 2), (1, 2), (2, 2)]
        # Allele frequencies are 0.5/0.5, hence expected HET=0.5;
        # two of four observed genotypes are heterozygous.
        self.assertAlmostEqual(stats._hwe(genotypes), binomtest(2, 4, 0.5).pvalue)


class TestSampleAndSchemaParsing(unittest.TestCase):
    def test_homozygous_compressed_ds_is_duplicated(self) -> None:
        sample = stats._sample_values(
            "1/1:5,5:(CAG)5:4,6:10", ["GT", "CN", "DS", "SD", "DP"], "CAG"
        )
        self.assertEqual(sample["cn_values"], [5.0, 5.0])
        self.assertEqual(sample["lps_values"], [5, 5])
        self.assertEqual(sample["pairs"], [(5.0, 5), (5.0, 5)])
        self.assertAlmostEqual(sample["allele_balance"], 0.4)

    def test_partial_genotype_keeps_correct_allele_alignment(self) -> None:
        sample = stats._sample_values(
            "./1:.,5:.,(CAG)5:.,6:6", ["GT", "CN", "DS", "SD", "DP"], "CAG"
        )
        self.assertEqual(sample["called_alleles"], 1)
        self.assertEqual(sample["cn_values"], [5.0])
        self.assertEqual(sample["lps_values"], [5])
        self.assertTrue(sample["partial"])

    def test_phantom_lpm_schema_is_repaired(self) -> None:
        format_text = "GT:AL:CN:LPM:AR:SD:DP:SN:SQ:MA:MR:DS:MV"
        sample = "0/1:9,12:3,4:A:4,6:10:0:.:.:.:(CAG)3,(CAG)4:."
        fields = [
            "chr1", "101", ".", "CAGCAGCAG", ".", ".", "PASS",
            "START=100;END=109;MOTIF=CAG;REFCN=3;ID=locus1",
            format_text, sample, "",
        ]
        data, repaired = stats._parse_locus(
            fields, line_number=10, samples=["sample1"], sample_indexes=[0],
            declared_formats={"GT", "AL", "CN", "AR", "SD", "DP", "SN", "SQ", "MA", "MR", "DS", "MV"},
        )
        self.assertTrue(repaired)
        self.assertEqual(data.cn_values, [3.0, 4.0])
        self.assertEqual(data.lps_values, [3, 4])
        self.assertEqual(data.depths, [10.0])
        self.assertEqual(data.allele_balances, [0.4])


@unittest.skipUnless(OUTPUT.is_file(), "full generated output is not present")
class TestFullOutputAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with OUTPUT.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            cls.fieldnames = list(reader.fieldnames or [])
            cls.rows = list(reader)

    def test_shape_and_unique_columns(self) -> None:
        self.assertEqual(len(self.rows), 9_778)
        self.assertEqual(len(self.fieldnames), 99)
        self.assertEqual(len(set(self.fieldnames)), 99)
        self.assertEqual(self.fieldnames, stats.OUTPUT_COLUMNS)

    def test_all_row_accounting_and_partitions(self) -> None:
        statuses: Counter[str] = Counter()
        totals = Counter()
        for row in self.rows:
            called = int(row["vcf_called_alleles"])
            cn = int(row["cn_total_genotyped_alleles"])
            lps = int(row["lps_total_genotyped_alleles"])
            paired = int(row["cn_lps_paired_alleles"])
            pure = int(row["fully_pure_allele_count"])
            interrupted = int(row["interrupted_allele_count"])
            exceeds = int(row["lps_exceeds_cn_count"])
            self.assertEqual(cn + int(row["missing_cn_alleles"]), called)
            self.assertEqual(lps + int(row["missing_lps_alleles"]), called)
            self.assertEqual(paired, lps)
            self.assertEqual(pure + interrupted + exceeds, paired)
            self.assertEqual(
                int(row["diploid_called_samples"])
                + int(row["haploid_called_samples"])
                + int(row["partial_genotype_samples"]),
                int(row["vcf_called_samples"]),
            )
            self.assertEqual(
                int(row["homozygous_samples"]) + int(row["heterozygous_samples"]),
                int(row["diploid_called_samples"]),
            )
            self.assertEqual(int(row["hwe_sample_count"]), int(row["diploid_called_samples"]))
            statuses[row["lps_status"]] += 1
            totals.update(called=called, cn=cn, lps=lps)
        self.assertEqual(totals, Counter(called=1_324_020, cn=1_324_020, lps=1_160_773))
        self.assertEqual(statuses, Counter(complete=8_556, partial=161, unavailable=1_061))

    def test_scientific_domains_and_lps_consistency(self) -> None:
        fractions = (
            "vcf_sample_call_rate", "cn_coverage", "lps_coverage",
            "cn_major_allele_frequency", "cn_heterozygosity",
            "cn_rare_allele_fraction", "cn_reference_frequency",
            "cn_non_reference_fraction", "cn_expanded_allele_fraction",
            "cn_contracted_allele_fraction", "fully_pure_allele_fraction",
            "interrupted_allele_fraction", "lps_exceeds_cn_fraction",
            "low_depth_sample_fraction", "allele_balance_outlier_fraction",
        )
        for row_index, row in enumerate(self.rows, start=1):
            for name in fractions:
                if row[name] != ".":
                    self.assertGreaterEqual(float(row[name]), 0.0, (row_index, name))
                    self.assertLessEqual(float(row[name]), 1.0, (row_index, name))
            if row["hwep"] != ".":
                self.assertGreaterEqual(float(row["hwep"]), 0.0, row_index)
                self.assertLessEqual(float(row["hwep"]), 1.0, row_index)
            for prefix in ("cn", "lps"):
                if row[f"{prefix}_min"] != ".":
                    ordered = [
                        float(row[f"{prefix}_min"]),
                        float(row[f"{prefix}_median"]),
                        float(row[f"{prefix}_99th_percentile"]),
                        float(row[f"{prefix}_max"]),
                    ]
                    self.assertEqual(ordered, sorted(ordered), (row_index, prefix))
                    self.assertGreaterEqual(float(row[f"{prefix}_variance"]), 0.0)
                    self.assertGreaterEqual(float(row[f"{prefix}_stdev"]), 0.0)
            # For a target-motif pure segment, LPS cannot exceed total CN.
            self.assertEqual(int(row["lps_exceeds_cn_count"]), 0, row_index)
            if row["cn_lps_mean_difference"] != ".":
                self.assertGreaterEqual(float(row["cn_lps_mean_difference"]), 0.0)
            missing = int(row["missing_lps_alleles"])
            measured = int(row["lps_total_genotyped_alleles"])
            expected_status = "complete" if missing == 0 else "partial" if measured else "unavailable"
            self.assertEqual(row["lps_status"], expected_status, row_index)

    def test_every_distribution_recomputed_from_histogram(self) -> None:
        tolerance = 5.1e-6
        for row_index, row in enumerate(self.rows, start=1):
            for prefix in ("cn", "lps"):
                with self.subTest(row=row_index, prefix=prefix):
                    values = expand_histogram(row[f"{prefix}_allele_histogram"])
                    self.assertEqual(len(values), int(row[f"{prefix}_total_genotyped_alleles"]))
                    self.assertEqual(len(set(values)), int(row[f"{prefix}_unique_allele_lengths"]))
                    if not values:
                        self.assertEqual(row[f"{prefix}_mean"], ".")
                        continue
                    counts = Counter(values)
                    total = len(values)
                    frequencies = {value: count / total for value, count in counts.items()}
                    mean = statistics.fmean(values)
                    expected = {
                        "min": min(values),
                        "median": statistics.median(values),
                        "mean": mean,
                        "99th_percentile": linear_percentile(values, 99),
                        "max": max(values),
                        "range": max(values) - min(values),
                        "variance": statistics.pvariance(values),
                        "stdev": statistics.pstdev(values),
                        "iqr": linear_percentile(values, 75) - linear_percentile(values, 25),
                        "mode": min(value for value, count in counts.items() if count == max(counts.values())),
                        "major_allele_frequency": max(frequencies.values()),
                        "effective_allele_count": 1 / sum(value * value for value in frequencies.values()),
                        "heterozygosity": 1 - sum(value * value for value in frequencies.values()),
                        "entropy": -sum(value * math.log2(value) for value in frequencies.values()),
                    }
                    for suffix, expected_value in expected.items():
                        self.assertAlmostEqual(
                            float(row[f"{prefix}_{suffix}"]), expected_value,
                            delta=tolerance, msg=f"row {row_index} {prefix}_{suffix}",
                        )
                    mean_value = statistics.fmean(values)
                    if mean_value == 0:
                        self.assertEqual(row[f"{prefix}_coefficient_of_variation"], ".")
                    else:
                        self.assertAlmostEqual(
                            float(row[f"{prefix}_coefficient_of_variation"]),
                            statistics.pstdev(values) / mean_value,
                            delta=tolerance,
                            msg=f"row {row_index} {prefix}_coefficient_of_variation",
                        )
                    singleton_types = sum(count == 1 for count in counts.values())
                    rare_observations = sum(
                        count for value, count in counts.items() if frequencies[value] < 0.01
                    )
                    self.assertEqual(int(row[f"{prefix}_singleton_allele_count"]), singleton_types)
                    self.assertEqual(int(row[f"{prefix}_rare_allele_count"]), rare_observations)
                    self.assertAlmostEqual(
                        float(row[f"{prefix}_rare_allele_fraction"]),
                        rare_observations / total, delta=tolerance,
                    )
                    if prefix == "cn":
                        reference = float(row["ref_cn"])
                        differences = [value - reference for value in values]
                        reference_count = sum(value == reference for value in values)
                        reference_expected = {
                            "cn_reference_frequency": reference_count / total,
                            "cn_non_reference_fraction": 1 - reference_count / total,
                            "cn_mean_difference_from_reference": statistics.fmean(differences),
                            "cn_median_difference_from_reference": statistics.median(differences),
                            "cn_max_expansion_from_reference": max(0.0, max(differences)),
                            "cn_max_contraction_from_reference": max(0.0, max(-value for value in differences)),
                            "cn_expanded_allele_fraction": sum(value > 0 for value in differences) / total,
                            "cn_contracted_allele_fraction": sum(value < 0 for value in differences) / total,
                        }
                        for name, expected_value in reference_expected.items():
                            self.assertAlmostEqual(
                                float(row[name]), expected_value, delta=tolerance,
                                msg=f"row {row_index} {name}",
                            )

    @unittest.skipUnless(VCF.is_file(), "source VCF is not present")
    def test_locus_metadata_matches_vcf_record_by_record(self) -> None:
        output_iter = iter(self.rows)
        with VCF.open(encoding="utf-8") as handle:
            records = (line for line in handle if line and not line.startswith("#"))
            for record_index, (line, row) in enumerate(zip(records, output_iter, strict=True), start=1):
                fields = line.rstrip("\r\n").split("\t")
                info = {
                    item.split("=", 1)[0]: item.split("=", 1)[1]
                    for item in fields[7].split(";") if "=" in item
                }
                self.assertEqual(row["chrom"], fields[0], record_index)
                self.assertEqual(row["start"], info.get("START", str(int(fields[1]) - 1)), record_index)
                self.assertEqual(row["end"], info.get("END", str(int(row["start"]) + len(fields[3]))), record_index)
                self.assertEqual(row["motif"], info["MOTIF"].upper(), record_index)
                expected_ref_cn = float(info.get("REFCN", len(fields[3]) / len(info["MOTIF"])))
                self.assertAlmostEqual(float(row["ref_cn"]), expected_ref_cn, delta=5.1e-6)

    @unittest.skipUnless(STAT_OUTPUT.is_file() and LPS_OUTPUT.is_file(), "standalone outputs are absent")
    def test_matches_previous_independent_tables(self) -> None:
        with STAT_OUTPUT.open(newline="", encoding="utf-8") as handle:
            stat_rows = list(csv.DictReader(handle, delimiter="\t"))
        with LPS_OUTPUT.open(newline="", encoding="utf-8") as handle:
            lps_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(stat_rows), len(self.rows))
        self.assertEqual(len(lps_rows), len(self.rows))
        for combined, old_cn, old_lps in zip(self.rows, stat_rows, lps_rows, strict=True):
            self.assertEqual(
                combined["cn_allele_histogram"].replace("x:", ":"), old_cn["acount"]
            )
            self.assertEqual(combined["lps_allele_histogram"], old_lps["lps_allele_histogram"])
            for combined_name, old_name in (
                ("cn_max", "thresh"), ("hwep", "hwep"),
                ("cn_heterozygosity", "het"), ("cn_entropy", "entropy"),
                ("cn_mean", "mean"), ("cn_mode", "mode"), ("cn_variance", "var"),
            ):
                assert_optional_close(self, combined[combined_name], old_cn[old_name], 5.1e-4)
            self.assertEqual(combined["vcf_called_samples"], old_cn["numcalled"])
            for name in ("lps_min", "lps_median", "lps_99th_percentile", "lps_max", "lps_mode", "lps_stdev"):
                assert_optional_close(self, combined[name], old_lps[name], 5.1e-3)
            for name in (
                "lps_unique_allele_lengths", "lps_total_genotyped_alleles",
                "missing_lps_alleles", "lps_status",
            ):
                self.assertEqual(combined[name], old_lps[name])


if __name__ == "__main__":
    unittest.main()
