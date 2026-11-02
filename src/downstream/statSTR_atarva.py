"""Compute per-locus tandem-repeat statistics from an ATARVA merged VCF.

This is an ATARVA-native adaptation of statSTR.  It reads repeat copy numbers
from the FORMAT/CN field and coordinates from INFO/START and INFO/END, avoiding
TRTools' VCF-type inference (ATARVA is not a TRTools harmonizer type).
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, TextIO


MAXPLOTS = 10


@dataclass(frozen=True)
class Locus:
    chrom: str
    start: int
    end: int
    motif: str
    allele_sequences: tuple[str, ...]
    sample_copy_numbers: tuple[tuple[float, ...], ...]


def _warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def _parse_info(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    if text == ".":
        return result
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        elif item:
            result[item] = ""
    return result


def _parse_copy_number(value: str) -> Optional[float]:
    if value in {"", ".", "None", "nan", "NaN"}:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def _extract_sample_copy_numbers(format_keys: Sequence[str], sample: str) -> tuple[float, ...]:
    if sample in {"", "."}:
        return ()
    values = sample.split(":")
    fields = dict(zip(format_keys, values))
    genotype = fields.get("GT", ".")
    if genotype in {"", ".", "./.", ".|."} or "." in genotype.replace("|", "/").split("/"):
        return ()
    cn_text = fields.get("CN", ".")
    copy_numbers = tuple(
        number
        for number in (_parse_copy_number(item) for item in cn_text.split(","))
        if number is not None
    )
    return copy_numbers


def _canonical_number(value: float) -> str:
    return str(int(value)) if value.is_integer() else f"{value:g}"


def _selected_calls(
    locus: Locus, sample_indexes: Optional[Sequence[int]]
) -> list[tuple[float, ...]]:
    if sample_indexes is None:
        return [call for call in locus.sample_copy_numbers if call]
    return [locus.sample_copy_numbers[index] for index in sample_indexes if locus.sample_copy_numbers[index]]


def _allele_counts(locus: Locus, sample_indexes: Optional[Sequence[int]]) -> Counter[float]:
    counts: Counter[float] = Counter()
    for call in _selected_calls(locus, sample_indexes):
        counts.update(call)
    return counts


def _allele_frequencies(counts: Counter[float]) -> dict[float, float]:
    total = sum(counts.values())
    if total == 0:
        return {}
    return {allele: count / total for allele, count in counts.items()}


def _serialize_distribution(counts: Counter[float], as_counts: bool) -> str:
    if not counts:
        return "."
    if as_counts:
        return ",".join(
            f"{_canonical_number(allele)}:{counts[allele]}" for allele in sorted(counts)
        )
    frequencies = _allele_frequencies(counts)
    return ",".join(
        f"{_canonical_number(allele)}:{frequencies[allele]:.3f}" for allele in sorted(counts)
    )


def _heterozygosity(frequencies: dict[float, float]) -> float:
    if not frequencies:
        return math.nan
    return 1.0 - sum(value * value for value in frequencies.values())


def _entropy(frequencies: dict[float, float]) -> float:
    if not frequencies:
        return math.nan
    return -sum(value * math.log2(value) for value in frequencies.values() if value > 0)


def _mean(frequencies: dict[float, float]) -> float:
    if not frequencies:
        return math.nan
    return sum(allele * frequency for allele, frequency in frequencies.items())


def _mode(frequencies: dict[float, float]) -> float:
    if not frequencies:
        return math.nan
    maximum = max(frequencies.values())
    return min(allele for allele, frequency in frequencies.items() if frequency == maximum)


def _variance(frequencies: dict[float, float]) -> float:
    if not frequencies:
        return math.nan
    mean = _mean(frequencies)
    return sum(frequency * (allele - mean) ** 2 for allele, frequency in frequencies.items())


def _binomial_two_sided(k: int, n: int, probability: float) -> float:
    if n == 0:
        return math.nan
    if probability <= 0:
        return 1.0 if k == 0 else 0.0
    if probability >= 1:
        return 1.0 if k == n else 0.0

    def pmf(successes: int) -> float:
        return (
            math.comb(n, successes)
            * probability**successes
            * (1.0 - probability) ** (n - successes)
        )

    observed = pmf(k)
    tolerance = observed * 1e-12
    return min(1.0, sum(pmf(i) for i in range(n + 1) if pmf(i) <= observed + tolerance))


def _hwe_pvalue(locus: Locus, sample_indexes: Optional[Sequence[int]]) -> float:
    diploid_calls = [call for call in _selected_calls(locus, sample_indexes) if len(call) == 2]
    if not diploid_calls:
        return math.nan
    counts: Counter[float] = Counter()
    for call in diploid_calls:
        counts.update(call)
    frequencies = _allele_frequencies(counts)
    expected_heterozygosity = _heterozygosity(frequencies)
    observed_heterozygotes = sum(call[0] != call[1] for call in diploid_calls)
    return _binomial_two_sided(
        observed_heterozygotes, len(diploid_calls), expected_heterozygosity
    )


def _format_float(value: float, precision: int) -> str:
    return "nan" if math.isnan(value) else f"{value:.{precision}f}"


def _group_headers(name: str, prefixes: Sequence[str]) -> list[str]:
    return [name] if not prefixes else [f"{name}-{prefix}" for prefix in prefixes]


def _read_sample_groups(
    sample_argument: Optional[str],
    prefix_argument: Optional[str],
    vcf_samples: Sequence[str],
) -> tuple[list[Optional[list[int]]], list[str]]:
    if not sample_argument:
        return [None], []
    paths = sample_argument.split(",")
    prefixes = prefix_argument.split(",") if prefix_argument else [str(i + 1) for i in range(len(paths))]
    if len(paths) != len(prefixes):
        raise ValueError("--sample-prefixes must have the same number of entries as --samples")
    sample_to_index = {sample: index for index, sample in enumerate(vcf_samples)}
    groups: list[Optional[list[int]]] = []
    for path in paths:
        with open(path, encoding="utf-8") as handle:
            requested = [line.strip() for line in handle if line.strip()]
        indexes = [sample_to_index[sample] for sample in requested if sample in sample_to_index]
        if not indexes:
            raise ValueError(f"No samples from {path} were found in the VCF")
        groups.append(indexes)
    return groups, prefixes


def _parse_locus(fields: Sequence[str], line_number: int, sample_count: int) -> Locus:
    expected_columns = 9 + sample_count
    # ATARVA merge currently leaves a trailing tab on records.  It is not a
    # sample column, so discard only empty surplus fields at the line end.
    normalized_fields = list(fields)
    while len(normalized_fields) > expected_columns and normalized_fields[-1] == "":
        normalized_fields.pop()
    fields = normalized_fields
    if len(fields) != expected_columns:
        raise ValueError(
            f"line {line_number}: expected {expected_columns} columns, found {len(fields)}"
        )
    chrom, pos_text, _record_id, ref, alt, _qual, _filter, info_text, format_text = fields[:9]
    try:
        pos = int(pos_text)
    except ValueError as error:
        raise ValueError(f"line {line_number}: invalid POS value {pos_text!r}") from error
    info = _parse_info(info_text)
    try:
        start = int(info.get("START", pos - 1))
        end = int(info.get("END", start + len(ref)))
    except ValueError as error:
        raise ValueError(f"line {line_number}: invalid START or END") from error
    motif = info.get("MOTIF", "")
    if not motif:
        raise ValueError(f"line {line_number}: missing INFO/MOTIF")
    format_keys = format_text.split(":")
    if "CN" not in format_keys or "GT" not in format_keys:
        raise ValueError(f"line {line_number}: FORMAT must contain GT and CN")
    calls = tuple(_extract_sample_copy_numbers(format_keys, sample) for sample in fields[9:])
    alleles = (ref,) + (() if alt == "." else tuple(alt.split(",")))
    return Locus(chrom, start, end, motif, alleles, calls)


def _write_locus(
    output: TextIO,
    locus: Locus,
    groups: Sequence[Optional[Sequence[int]]],
    args: argparse.Namespace,
) -> None:
    row = [locus.chrom, str(locus.start), str(locus.end)]
    summaries = []
    for group in groups:
        counts = _allele_counts(locus, group)
        frequencies = _allele_frequencies(counts)
        calls = _selected_calls(locus, group)
        summaries.append((group, counts, frequencies, calls))
    if args.thresh:
        row.extend(_format_float(max(counts) if counts else math.nan, args.precision) for _, counts, _, _ in summaries)
    if args.afreq:
        row.extend(_serialize_distribution(counts, as_counts=False) for _, counts, _, _ in summaries)
    if args.acount:
        row.extend(_serialize_distribution(counts, as_counts=True) for _, counts, _, _ in summaries)
    if args.nalleles:
        row.extend(str(sum(frequency >= args.nalleles_thresh for frequency in frequencies.values())) for _, _, frequencies, _ in summaries)
    if args.hwep:
        row.extend(_format_float(_hwe_pvalue(locus, group), args.precision) for group, _, _, _ in summaries)
    if args.het:
        row.extend(_format_float(_heterozygosity(frequencies), args.precision) for _, _, frequencies, _ in summaries)
    if args.entropy:
        row.extend(_format_float(_entropy(frequencies), args.precision) for _, _, frequencies, _ in summaries)
    if args.mean:
        row.extend(_format_float(_mean(frequencies), args.precision) for _, _, frequencies, _ in summaries)
    if args.mode:
        row.extend(_format_float(_mode(frequencies), args.precision) for _, _, frequencies, _ in summaries)
    if args.var:
        row.extend(_format_float(_variance(frequencies), args.precision) for _, _, frequencies, _ in summaries)
    if args.numcalled:
        row.extend(str(len(calls)) for _, _, _, calls in summaries)
    output.write("\t".join(row) + "\n")


def _build_header(args: argparse.Namespace, prefixes: Sequence[str]) -> list[str]:
    header = ["chrom", "start", "end"]
    for enabled, name in (
        (args.thresh, "thresh"),
        (args.afreq, "afreq"),
        (args.acount, "acount"),
        (args.nalleles, "nalleles"),
        (args.hwep, "hwep"),
        (args.het, "het"),
        (args.entropy, "entropy"),
        (args.mean, "mean"),
        (args.mode, "mode"),
        (args.var, "var"),
        (args.numcalled, "numcalled"),
    ):
        if enabled:
            header.extend(_group_headers(name, prefixes))
    return header


def _iter_vcf(handle: Iterable[str]) -> tuple[list[str], Iterable[tuple[int, list[str]]]]:
    iterator = enumerate(handle, 1)
    for line_number, line in iterator:
        if line.startswith("#CHROM\t"):
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) < 10:
                raise ValueError("VCF contains no samples")
            samples = columns[9:]

            def records() -> Iterable[tuple[int, list[str]]]:
                for record_line_number, record_line in iterator:
                    if not record_line.strip() or record_line.startswith("#"):
                        continue
                    yield record_line_number, record_line.rstrip("\r\n").split("\t")

            return samples, records()
    raise ValueError("VCF #CHROM header was not found")


def getargs(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input ATARVA merged VCF")
    parser.add_argument("--out", required=True, help="Output prefix, or stdout")
    parser.add_argument("--precision", type=int, default=3)
    parser.add_argument("--samples", help="Sample-list file, or comma-separated files")
    parser.add_argument("--sample-prefixes", help="Comma-separated names for sample groups")
    parser.add_argument("--only-passing", action="store_true", help="Process only FILTER=PASS records")
    parser.add_argument("--thresh", action="store_true")
    parser.add_argument("--afreq", action="store_true")
    parser.add_argument("--acount", action="store_true")
    parser.add_argument("--nalleles", action="store_true")
    parser.add_argument("--nalleles-thresh", type=float, default=0.01)
    parser.add_argument("--hwep", action="store_true")
    parser.add_argument("--het", action="store_true")
    parser.add_argument("--entropy", action="store_true")
    parser.add_argument("--mean", action="store_true")
    parser.add_argument("--mode", action="store_true")
    parser.add_argument("--var", action="store_true")
    parser.add_argument("--numcalled", action="store_true")
    parser.add_argument(
        "--use-length",
        action="store_true",
        help="Compatibility flag; ATARVA FORMAT/CN is always used",
    )
    parser.add_argument("--plot-afreq", action="store_true", help="Not supported by this adaptation")
    args = parser.parse_args(argv)
    stat_names = ("thresh", "afreq", "acount", "nalleles", "hwep", "het", "entropy", "mean", "mode", "var", "numcalled")
    if not any(getattr(args, name) for name in stat_names):
        parser.error("select at least one statistic")
    if args.plot_afreq:
        parser.error("--plot-afreq is not supported; the local environment lacks matplotlib")
    if not 0 <= args.nalleles_thresh <= 1:
        parser.error("--nalleles-thresh must be between 0 and 1")
    if args.precision < 0:
        parser.error("--precision must be nonnegative")
    return args


def main(args: argparse.Namespace) -> int:
    if not os.path.isfile(args.vcf):
        _warn(f"{args.vcf} does not exist or is not a file")
        return 1
    if args.out != "stdout":
        output_path = args.out + ".tab"
        output_directory = os.path.dirname(os.path.abspath(output_path))
        if not os.path.isdir(output_directory):
            _warn(f"output directory {output_directory} does not exist")
            return 1
    else:
        output_path = None

    start_time = time.time()
    processed = 0
    skipped = 0
    try:
        with open(args.vcf, encoding="utf-8") as input_handle:
            samples, records = _iter_vcf(input_handle)
            groups, prefixes = _read_sample_groups(
                args.samples, args.sample_prefixes, samples
            )
            output_handle: TextIO
            if output_path is None:
                output_handle = sys.stdout
            else:
                output_handle = open(output_path, "w", encoding="utf-8", newline="")
            try:
                output_handle.write("\t".join(_build_header(args, prefixes)) + "\n")
                for line_number, fields in records:
                    if args.only_passing and fields[6] != "PASS":
                        skipped += 1
                        continue
                    locus = _parse_locus(fields, line_number, len(samples))
                    _write_locus(output_handle, locus, groups, args)
                    processed += 1
                    if processed % 1000 == 0:
                        output_handle.flush()
                        if output_path is not None:
                            elapsed = time.time() - start_time
                            print(
                                f"Processed {processed:,} records ({elapsed / processed:.5f} sec/record)",
                                end="\r",
                                flush=True,
                            )
            finally:
                if output_path is not None:
                    output_handle.close()
    except (IndexError, OSError, ValueError) as error:
        _warn(str(error))
        return 1

    if output_path is not None:
        elapsed = time.time() - start_time
        print(
            f"Processed {processed:,} records in {elapsed:.2f}s; skipped {skipped:,}. "
            f"Output: {output_path}"
        )
    return 0


def run() -> None:
    raise SystemExit(main(getargs()))


if __name__ == "__main__":
    run()
