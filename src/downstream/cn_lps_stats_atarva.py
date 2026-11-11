"""Create an extended per-locus CN and ATaRVa-DS-derived LPS statistics table.

The VCF is read once.  Total repeat copy-number (CN) and longest pure segment
(LPS) statistics are reported in separate, explicitly prefixed column groups.
LPS is calculated from FORMAT/DS for the locus motif, including cyclic rotations
and reverse-complement rotations; it is not a TRGT-LPS HMM score.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
import time
from collections import Counter
from dataclasses import dataclass
from typing import Optional, Sequence, TextIO

import numpy as np

from lpsSTR_atarva import (
    _called_genotype_alleles,
    _effective_format_keys,
    _iter_vcf,
    _parse_float,
    _parse_info,
    _selected_sample_indexes,
    _split_allele_values,
    longest_target_lps,
)


LOCUS_COLUMNS = [
    "chrom", "start", "end", "locus_id", "motif", "motif_length", "ref_cn",
]

SUMMARY_SUFFIXES = [
    "min", "median", "mean", "99th_percentile", "max", "range", "variance",
    "stdev", "coefficient_of_variation", "iqr", "mode", "unique_allele_lengths",
    "total_genotyped_alleles", "allele_histogram", "major_allele_frequency",
    "effective_allele_count", "heterozygosity", "entropy",
    "singleton_allele_count", "rare_allele_count", "rare_allele_fraction",
]

CN_REFERENCE_COLUMNS = [
    "cn_reference_frequency", "cn_non_reference_fraction",
    "cn_mean_difference_from_reference", "cn_median_difference_from_reference",
    "cn_max_expansion_from_reference", "cn_max_contraction_from_reference",
    "cn_expanded_allele_fraction", "cn_contracted_allele_fraction",
]

COMPARISON_COLUMNS = [
    "cn_lps_paired_alleles", "cn_lps_mean_difference",
    "cn_lps_median_difference", "cn_lps_99th_difference",
    "cn_lps_max_difference", "cn_lps_mean_ratio", "cn_lps_median_ratio",
    "fully_pure_allele_count", "fully_pure_allele_fraction",
    "interrupted_allele_count", "interrupted_allele_fraction",
    "lps_exceeds_cn_count", "lps_exceeds_cn_fraction",
]

QC_COLUMNS = [
    "vcf_total_samples", "vcf_called_samples", "vcf_sample_call_rate",
    "vcf_called_alleles", "diploid_called_samples", "haploid_called_samples",
    "homozygous_samples", "heterozygous_samples", "heterozygous_sample_fraction",
    "partial_genotype_samples", "missing_cn_alleles", "missing_lps_alleles",
    "cn_coverage", "lps_coverage", "lps_status", "hwe_sample_count", "hwep",
    "depth_measured_samples", "mean_depth", "median_depth", "min_depth", "max_depth",
    "low_depth_sample_count", "low_depth_sample_fraction",
    "allele_balance_measured_samples", "mean_allele_balance", "median_allele_balance",
    "allele_balance_outlier_count", "allele_balance_outlier_fraction",
]

OUTPUT_COLUMNS = (
    LOCUS_COLUMNS
    + [f"cn_{suffix}" for suffix in SUMMARY_SUFFIXES]
    + CN_REFERENCE_COLUMNS
    + [f"lps_{suffix}" for suffix in SUMMARY_SUFFIXES]
    + COMPARISON_COLUMNS
    + QC_COLUMNS
)


@dataclass
class LocusData:
    chrom: str
    start: int
    end: int
    locus_id: str
    motif: str
    ref_cn: float
    cn_values: list[float]
    lps_values: list[int]
    paired_values: list[tuple[float, int]]
    called_samples: int
    called_alleles: int
    diploid_called_samples: int
    haploid_called_samples: int
    homozygous_samples: int
    heterozygous_samples: int
    partial_genotype_samples: int
    missing_cn_alleles: int
    missing_lps_alleles: int
    diploid_cn_genotypes: list[tuple[float, float]]
    depths: list[float]
    allele_balances: list[float]


@dataclass
class MissingMotifLocus:
    chrom: str
    start: int
    end: int
    locus_id: str
    motif: str


def _warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def _fmt(value: float, precision: int = 6) -> str:
    if not math.isfinite(float(value)):
        return "."
    if float(value).is_integer():
        return str(int(value))
    return f"{float(value):.{precision}f}".rstrip("0").rstrip(".")


def _canonical(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else f"{value:g}"


def _histogram(values: Sequence[float]) -> str:
    counts = Counter(values)
    return ",".join(f"{_canonical(value)}x:{counts[value]}" for value in sorted(counts)) or "."


def _distribution_summary(
    values: Sequence[float], rare_threshold: float
) -> list[str]:
    if not values:
        return [".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".",
                "0", "0", ".", ".", ".", ".", ".", "0", "0", "."]
    array = np.asarray(values, dtype=float)
    counts = Counter(values)
    total = len(values)
    frequencies = {value: count / total for value, count in counts.items()}
    maximum_count = max(counts.values())
    mode = min(value for value, count in counts.items() if count == maximum_count)
    heterozygosity = 1.0 - sum(frequency**2 for frequency in frequencies.values())
    entropy = -sum(frequency * math.log2(frequency) for frequency in frequencies.values())
    effective = 1.0 / sum(frequency**2 for frequency in frequencies.values())
    rare_count = sum(count for value, count in counts.items() if frequencies[value] < rare_threshold)
    q25, q75 = np.percentile(array, [25, 75])
    mean = float(np.mean(array))
    stdev = float(np.std(array, ddof=0))
    return [
        _fmt(float(np.min(array))),
        _fmt(float(np.median(array))),
        _fmt(mean),
        _fmt(float(np.percentile(array, 99))),
        _fmt(float(np.max(array))),
        _fmt(float(np.max(array) - np.min(array))),
        _fmt(float(np.var(array, ddof=0))),
        _fmt(stdev),
        _fmt(stdev / mean) if mean != 0 else ".",
        _fmt(float(q75 - q25)),
        _fmt(float(mode)),
        str(len(counts)),
        str(total),
        _histogram(values),
        _fmt(max(frequencies.values())),
        _fmt(effective),
        _fmt(heterozygosity),
        _fmt(entropy),
        str(sum(count == 1 for count in counts.values())),
        str(rare_count),
        _fmt(rare_count / total),
    ]


def _cn_reference_summary(values: Sequence[float], ref_cn: float) -> list[str]:
    if not values:
        return ["."] * len(CN_REFERENCE_COLUMNS)
    array = np.asarray(values, dtype=float)
    differences = array - ref_cn
    reference_count = int(np.sum(array == ref_cn))
    total = len(values)
    return [
        _fmt(reference_count / total),
        _fmt(1.0 - reference_count / total),
        _fmt(float(np.mean(differences))),
        _fmt(float(np.median(differences))),
        _fmt(max(0.0, float(np.max(differences)))),
        _fmt(max(0.0, float(np.max(-differences)))),
        _fmt(float(np.sum(differences > 0)) / total),
        _fmt(float(np.sum(differences < 0)) / total),
    ]


def _comparison_summary(pairs: Sequence[tuple[float, int]]) -> list[str]:
    if not pairs:
        return [
            "0", ".", ".", ".", ".", ".", ".",
            "0", ".", "0", ".", "0", ".",
        ]
    cn = np.asarray([pair[0] for pair in pairs], dtype=float)
    lps = np.asarray([pair[1] for pair in pairs], dtype=float)
    differences = cn - lps
    nonzero = cn != 0
    ratios = lps[nonzero] / cn[nonzero]
    total = len(pairs)
    pure = int(np.sum(np.isclose(cn, lps)))
    interrupted = int(np.sum(lps < cn))
    exceeds = int(np.sum(lps > cn))
    return [
        str(total),
        _fmt(float(np.mean(differences))),
        _fmt(float(np.median(differences))),
        _fmt(float(np.percentile(differences, 99))),
        _fmt(float(np.max(differences))),
        _fmt(float(np.mean(ratios))) if len(ratios) else ".",
        _fmt(float(np.median(ratios))) if len(ratios) else ".",
        str(pure),
        _fmt(pure / total),
        str(interrupted),
        _fmt(interrupted / total),
        str(exceeds),
        _fmt(exceeds / total),
    ]


def _binomial_two_sided(successes: int, trials: int, probability: float) -> float:
    if trials == 0:
        return math.nan
    if probability <= 0:
        return 1.0 if successes == 0 else 0.0
    if probability >= 1:
        return 1.0 if successes == trials else 0.0

    def pmf(value: int) -> float:
        return (
            math.comb(trials, value)
            * probability**value
            * (1.0 - probability) ** (trials - value)
        )

    observed = pmf(successes)
    tolerance = observed * 1e-12
    return min(1.0, sum(pmf(value) for value in range(trials + 1)
                        if pmf(value) <= observed + tolerance))


def _hwe(genotypes: Sequence[tuple[float, float]]) -> float:
    if not genotypes:
        return math.nan
    counts: Counter[float] = Counter()
    for genotype in genotypes:
        counts.update(genotype)
    total_alleles = sum(counts.values())
    frequencies = [count / total_alleles for count in counts.values()]
    expected_heterozygosity = 1.0 - sum(frequency**2 for frequency in frequencies)
    observed_heterozygotes = sum(first != second for first, second in genotypes)
    return _binomial_two_sided(observed_heterozygotes, len(genotypes), expected_heterozygosity)


def _sample_values(
    sample_text: str, format_keys: Sequence[str], motif: str
) -> dict[str, object]:
    fields = dict(zip(format_keys, sample_text.split(":"))) if sample_text not in {"", "."} else {}
    genotype_text = fields.get("GT", ".")
    raw_gt = genotype_text.replace("|", "/").split("/") if genotype_text not in {"", "."} else []
    called_gt = _called_genotype_alleles(genotype_text)
    if not called_gt:
        return {"called": False}

    def align_to_raw_gt(value: str) -> list[Optional[str]]:
        items = _split_allele_values(value)
        if len(items) == len(raw_gt):
            aligned = list(items)
        elif len(items) == 1 and len(raw_gt) > 1 and "." not in raw_gt and len(set(raw_gt)) == 1:
            aligned = [items[0]] * len(raw_gt)
        else:
            aligned = [items[index] if index < len(items) else None for index in range(len(raw_gt))]
        return [item for item, genotype_allele in zip(aligned, raw_gt) if genotype_allele != "."]

    cn_fields = align_to_raw_gt(fields.get("CN", "."))
    ds_fields = align_to_raw_gt(fields.get("DS", "."))
    cn_values: list[Optional[float]] = [_parse_float(value) for value in cn_fields]
    lps_values: list[Optional[int]] = []
    pairs: list[tuple[float, int]] = []
    for cn_value, ds_value in zip(cn_values, ds_fields):
        lps = longest_target_lps(ds_value, motif) if ds_value is not None else None
        if lps is None and cn_value == 0:
            lps = 0
        lps_values.append(lps)
        if cn_value is not None and lps is not None:
            pairs.append((cn_value, lps))

    depth = _parse_float(fields.get("DP"))
    support_values = [
        value for value in (_parse_float(item) for item in _split_allele_values(fields.get("SD", ".")))
        if value is not None
    ]
    allele_balance = None
    if len(support_values) >= 2 and sum(support_values[:2]) > 0:
        allele_balance = support_values[0] / sum(support_values[:2])

    complete_gt = bool(raw_gt) and "." not in raw_gt
    complete_cn = complete_gt and len(cn_values) == len(raw_gt) and all(value is not None for value in cn_values)
    cn_genotype = tuple(cn_values) if complete_cn else ()
    return {
        "called": True,
        "called_alleles": len(called_gt),
        "cn_values": [value for value in cn_values if value is not None],
        "lps_values": [value for value in lps_values if value is not None],
        "pairs": pairs,
        "missing_cn": sum(value is None for value in cn_values),
        "missing_lps": sum(value is None for value in lps_values),
        "partial": not complete_cn,
        "cn_genotype": cn_genotype,
        "depth": depth,
        "allele_balance": allele_balance,
    }


def _parse_locus(
    raw_fields: Sequence[str],
    line_number: int,
    samples: Sequence[str],
    sample_indexes: Sequence[int],
    declared_formats: set[str],
) -> tuple[LocusData, bool]:
    expected_columns = 9 + len(samples)
    fields = list(raw_fields)
    while len(fields) > expected_columns and fields[-1] == "":
        fields.pop()
    if len(fields) != expected_columns:
        raise ValueError(f"line {line_number}: expected {expected_columns} columns, found {len(fields)}")
    chrom, pos_text, record_id, ref, _alt, _qual, _filter, info_text, format_text = fields[:9]
    try:
        pos = int(pos_text)
    except ValueError as error:
        raise ValueError(f"line {line_number}: invalid POS {pos_text!r}") from error
    info = _parse_info(info_text)
    motif = info.get("MOTIF", "").upper()
    if not motif:
        raise ValueError(f"line {line_number}: missing INFO/MOTIF")
    try:
        start = int(info.get("START", pos - 1))
        end = int(info.get("END", start + len(ref)))
        ref_cn = float(info.get("REFCN", len(ref) / len(motif)))
    except ValueError as error:
        raise ValueError(f"line {line_number}: invalid START, END or REFCN") from error
    locus_id = info.get("ID") or (record_id if record_id != "." else f"{chrom}-{start}-{end}-{motif}")
    format_keys, repaired = _effective_format_keys(format_text, declared_formats)
    for required in ("GT", "CN", "DS", "DP", "SD"):
        if required not in format_keys:
            raise ValueError(f"line {line_number}: effective FORMAT lacks {required}")

    cn_values: list[float] = []
    lps_values: list[int] = []
    pairs: list[tuple[float, int]] = []
    called_samples = called_alleles = diploid = haploid = homozygous = heterozygous = partial = 0
    missing_cn = missing_lps = 0
    diploid_genotypes: list[tuple[float, float]] = []
    depths: list[float] = []
    balances: list[float] = []
    for sample_index in sample_indexes:
        sample = _sample_values(fields[9 + sample_index], format_keys, motif)
        if not sample.get("called"):
            continue
        called_samples += 1
        called_alleles += int(sample["called_alleles"])
        cn_values.extend(sample["cn_values"])
        lps_values.extend(sample["lps_values"])
        pairs.extend(sample["pairs"])
        missing_cn += int(sample["missing_cn"])
        missing_lps += int(sample["missing_lps"])
        partial += int(sample["partial"])
        genotype = sample["cn_genotype"]
        if len(genotype) == 2:
            diploid += 1
            pair = (float(genotype[0]), float(genotype[1]))
            diploid_genotypes.append(pair)
            if pair[0] == pair[1]:
                homozygous += 1
            else:
                heterozygous += 1
        elif len(genotype) == 1:
            haploid += 1
        if sample["depth"] is not None:
            depths.append(float(sample["depth"]))
        if sample["allele_balance"] is not None:
            balances.append(float(sample["allele_balance"]))

    return LocusData(
        chrom, start, end, locus_id, motif, ref_cn, cn_values, lps_values, pairs,
        called_samples, called_alleles, diploid, haploid, homozygous, heterozygous,
        partial, missing_cn, missing_lps, diploid_genotypes, depths, balances,
    ), repaired


def _parse_info_dot_locus(
    raw_fields: Sequence[str], line_number: int, sample_count: int
) -> Optional[MissingMotifLocus]:
    """Parse the identity fields for the special all-missing INFO=. records."""
    expected_columns = 9 + sample_count
    fields = list(raw_fields)
    while len(fields) > expected_columns and fields[-1] == "":
        fields.pop()
    if len(fields) != expected_columns:
        raise ValueError(
            f"line {line_number}: expected {expected_columns} columns, found {len(fields)}"
        )

    chrom, _pos, record_id, _ref, _alt, _qual, _filter, info_text = fields[:8]
    if info_text != ".":
        return None
    if record_id in {"", "."}:
        raise ValueError(
            f"line {line_number}: INFO=. record has no VCF ID from which to extract "
            "start, end and motif"
        )

    match = re.fullmatch(r"(.+)-(\d+)-(\d+)-([A-Za-z]+)", record_id)
    if match is None:
        raise ValueError(
            f"line {line_number}: cannot extract chrom-start-end-motif from "
            f"VCF ID {record_id!r} for INFO=. record"
        )
    id_chrom, start_text, end_text, motif = match.groups()
    if id_chrom != chrom:
        raise ValueError(
            f"line {line_number}: VCF ID chromosome {id_chrom!r} does not match "
            f"CHROM {chrom!r} for INFO=. record"
        )
    start = int(start_text)
    end = int(end_text)
    if end < start:
        raise ValueError(
            f"line {line_number}: VCF ID end {end} precedes start {start} "
            "for INFO=. record"
        )
    return MissingMotifLocus(chrom, start, end, record_id, motif.upper())


def _qc_fields(data: LocusData, args: argparse.Namespace, total_samples: int) -> list[str]:
    genotype_denominator = data.homozygous_samples + data.heterozygous_samples
    low_depth = sum(depth < args.low_depth_threshold for depth in data.depths)
    balance_outliers = sum(
        balance < args.allele_balance_min or balance > args.allele_balance_max
        for balance in data.allele_balances
    )
    lps_status = (
        "complete" if data.missing_lps_alleles == 0
        else "partial" if data.lps_values else "unavailable"
    )
    return [
        str(total_samples),
        str(data.called_samples),
        _fmt(data.called_samples / total_samples if total_samples else math.nan),
        str(data.called_alleles),
        str(data.diploid_called_samples),
        str(data.haploid_called_samples),
        str(data.homozygous_samples),
        str(data.heterozygous_samples),
        _fmt(data.heterozygous_samples / genotype_denominator if genotype_denominator else math.nan),
        str(data.partial_genotype_samples),
        str(data.missing_cn_alleles),
        str(data.missing_lps_alleles),
        _fmt(len(data.cn_values) / data.called_alleles if data.called_alleles else math.nan),
        _fmt(len(data.lps_values) / data.called_alleles if data.called_alleles else math.nan),
        lps_status,
        str(len(data.diploid_cn_genotypes)),
        _fmt(_hwe(data.diploid_cn_genotypes)),
        str(len(data.depths)),
        _fmt(float(np.mean(data.depths))) if data.depths else ".",
        _fmt(float(np.median(data.depths))) if data.depths else ".",
        _fmt(min(data.depths)) if data.depths else ".",
        _fmt(max(data.depths)) if data.depths else ".",
        str(low_depth),
        _fmt(low_depth / len(data.depths) if data.depths else math.nan),
        str(len(data.allele_balances)),
        _fmt(float(np.mean(data.allele_balances))) if data.allele_balances else ".",
        _fmt(float(np.median(data.allele_balances))) if data.allele_balances else ".",
        str(balance_outliers),
        _fmt(balance_outliers / len(data.allele_balances) if data.allele_balances else math.nan),
    ]


def _write_locus(output: TextIO, data: LocusData, args: argparse.Namespace, total_samples: int) -> None:
    row = [
        data.chrom, str(data.start), str(data.end), data.locus_id, data.motif,
        str(len(data.motif)), _fmt(data.ref_cn),
        *_distribution_summary(data.cn_values, args.rare_allele_threshold),
        *_cn_reference_summary(data.cn_values, data.ref_cn),
        *_distribution_summary(data.lps_values, args.rare_allele_threshold),
        *_comparison_summary(data.paired_values),
        *_qc_fields(data, args, total_samples),
    ]
    if len(row) != len(OUTPUT_COLUMNS):
        raise ValueError(f"internal error: produced {len(row)} fields for {len(OUTPUT_COLUMNS)} columns")
    output.write("\t".join(row) + "\n")


def _write_missing_motif_locus(output: TextIO, data: MissingMotifLocus) -> None:
    row = [
        data.chrom,
        str(data.start),
        str(data.end),
        data.locus_id,
        data.motif,
        str(len(data.motif)),
        "NA",
    ]
    row.extend(["NA"] * (len(OUTPUT_COLUMNS) - len(row)))
    row[OUTPUT_COLUMNS.index("lps_status")] = "missing_motif"
    output.write("\t".join(row) + "\n")


def getargs(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input ATaRVa merged VCF")
    parser.add_argument("--out", required=True, help="Output TSV path or prefix")
    parser.add_argument("--samples", help="Optional exact sample-name list")
    parser.add_argument("--numrecords", type=int, help="Only process this many records")
    parser.add_argument("--rare-allele-threshold", type=float, default=0.01)
    parser.add_argument("--low-depth-threshold", type=float, default=10.0)
    parser.add_argument("--allele-balance-min", type=float, default=0.2)
    parser.add_argument("--allele-balance-max", type=float, default=0.8)
    args = parser.parse_args(argv)
    if args.numrecords is not None and args.numrecords < 0:
        parser.error("--numrecords must be nonnegative")
    if not 0 <= args.rare_allele_threshold <= 1:
        parser.error("--rare-allele-threshold must be between 0 and 1")
    if args.low_depth_threshold < 0:
        parser.error("--low-depth-threshold must be nonnegative")
    if not 0 <= args.allele_balance_min < args.allele_balance_max <= 1:
        parser.error("allele-balance thresholds must satisfy 0 <= min < max <= 1")
    return args


def main(args: argparse.Namespace) -> int:
    if not os.path.isfile(args.vcf):
        _warn(f"{args.vcf} does not exist or is not a file")
        return 1
    output_path = args.out if args.out.lower().endswith((".tab", ".tsv")) else args.out + ".tab"
    output_directory = os.path.dirname(os.path.abspath(output_path))
    if not os.path.isdir(output_directory):
        _warn(f"output directory {output_directory} does not exist")
        return 1

    processed = repairs = missing_motif_records = total_called = total_cn = total_lps = 0
    start_time = time.time()
    try:
        with open(args.vcf, encoding="utf-8") as input_handle:
            declared_formats, samples, records = _iter_vcf(input_handle)
            sample_indexes = _selected_sample_indexes(samples, args.samples)
            with open(output_path, "w", encoding="utf-8", newline="") as output_handle:
                output_handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
                for line_number, fields in records:
                    if args.numrecords is not None and processed >= args.numrecords:
                        break
                    missing_motif = _parse_info_dot_locus(fields, line_number, len(samples))
                    if missing_motif is not None:
                        _warn(
                            f"line {line_number}: INFO=. record retained with "
                            f"lps_status=missing_motif; statistics set to NA "
                            f"(locus_id={missing_motif.locus_id})"
                        )
                        _write_missing_motif_locus(output_handle, missing_motif)
                        processed += 1
                        missing_motif_records += 1
                        continue
                    data, repaired = _parse_locus(
                        fields, line_number, samples, sample_indexes, declared_formats
                    )
                    _write_locus(output_handle, data, args, len(sample_indexes))
                    processed += 1
                    repairs += int(repaired)
                    total_called += data.called_alleles
                    total_cn += len(data.cn_values)
                    total_lps += len(data.lps_values)
                    if processed % 1000 == 0:
                        output_handle.flush()
                        elapsed = time.time() - start_time
                        print(
                            f"Processed {processed:,} records ({elapsed / processed:.5f} sec/record)",
                            end="\r", flush=True,
                        )
    except (IndexError, OSError, ValueError) as error:
        _warn(str(error))
        return 1

    elapsed = time.time() - start_time
    print(
        f"Processed {processed:,} loci and {total_called:,} called alleles in {elapsed:.2f}s. "
        f"Measured CN for {total_cn:,} and LPS for {total_lps:,} alleles; "
        f"retained {missing_motif_records:,} INFO=. missing-motif records; "
        f"repaired {repairs:,} phantom-LPM schemas. Output: {output_path}"
    )
    return 0


def run() -> None:
    raise SystemExit(main(getargs()))


if __name__ == "__main__":
    run()
