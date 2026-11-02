"""Compute TRExplorer-style LPS summary statistics from an ATaRVa VCF.

LPS is derived from ATaRVa's FORMAT/DS motif decomposition.  For each called
allele, the script finds the longest uninterrupted run of the locus motif,
accepting cyclic rotations and reverse-complement rotations.  It then applies
the current TRExplorer/str-analysis summary definitions: NumPy median and 99th
percentile, population standard deviation (ddof=0), and the smallest value when
multiple modes are tied.

This produces an ATaRVa-DS-derived LPS, not a TRGT-LPS HMM score.
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
from typing import Iterable, Optional, Sequence, TextIO

import numpy as np


DECOMPOSED_RUN = re.compile(r"^\(([A-Za-z]+)\)(\d+)$")
DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class LocusResult:
    chrom: str
    start: int
    end: int
    locus_id: str
    motif: str
    lps_values: tuple[int, ...]
    vcf_called_alleles: int
    missing_lps_alleles: int


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


def _rotations(sequence: str) -> set[str]:
    sequence = sequence.upper()
    return {sequence[index:] + sequence[:index] for index in range(len(sequence))}


def _equivalent_motifs(motif: str) -> set[str]:
    motif = motif.upper()
    reverse_complement = motif.translate(DNA_COMPLEMENT)[::-1]
    return _rotations(motif) | _rotations(reverse_complement)


def _run_units_in_motif(segment_motif: str, target_motif: str, equivalents: set[str]) -> int:
    """Return target-motif units represented by one decomposed motif occurrence."""
    segment_motif = segment_motif.upper()
    target_length = len(target_motif)
    if len(segment_motif) % target_length:
        return 0
    units = len(segment_motif) // target_length
    for offset in range(target_length):
        rotated = segment_motif[offset:] + segment_motif[:offset]
        if all(
            rotated[index:index + target_length] in equivalents
            for index in range(0, len(rotated), target_length)
        ):
            return units
    return 0


def _longest_literal_run(sequence: str, motif: str, equivalents: set[str]) -> int:
    sequence = sequence.upper()
    motif_length = len(motif)
    maximum = 0
    for offset in range(len(sequence)):
        count = 0
        position = offset
        while (
            position + motif_length <= len(sequence)
            and sequence[position:position + motif_length] in equivalents
        ):
            count += 1
            position += motif_length
        maximum = max(maximum, count)
    return maximum


def longest_target_lps(decomposition: str, motif: str) -> Optional[int]:
    """Return the target motif's longest pure run, or None when DS is unavailable.

    A nonempty decomposition with no target-motif run has LPS zero.  This is
    distinct from an unavailable decomposition ("." or empty), which is None.
    """
    if decomposition in {"", ".", "None"}:
        return None
    equivalents = _equivalent_motifs(motif)
    maximum = 0
    for component in decomposition.split("-"):
        match = DECOMPOSED_RUN.fullmatch(component)
        if match:
            segment_motif, copies_text = match.groups()
            represented_units = _run_units_in_motif(segment_motif, motif, equivalents)
            if represented_units:
                maximum = max(maximum, represented_units * int(copies_text))
        elif component:
            maximum = max(maximum, _longest_literal_run(component, motif, equivalents))
    return maximum


def _called_genotype_alleles(genotype: str) -> list[str]:
    if genotype in {"", ".", "./.", ".|."}:
        return []
    return [allele for allele in genotype.replace("|", "/").split("/") if allele != "."]


def _split_allele_values(value: str) -> list[str]:
    if value in {"", "."}:
        return []
    return value.split(",")


def _align_values_to_genotype(
    values: Sequence[str], genotype_alleles: Sequence[str]
) -> list[Optional[str]]:
    """Align a FORMAT allele list to GT ploidy.

    ATaRVa writes one DS value for a diploid homozygous call, while GT, AL and
    CN contain two alleles.  Duplicate that one value only when all called GT
    indexes are identical.  Other length mismatches remain missing rather than
    inventing allele measurements.
    """
    ploidy = len(genotype_alleles)
    if ploidy == 0:
        return []
    if len(values) == ploidy:
        return list(values)
    if len(values) == 1 and ploidy > 1 and len(set(genotype_alleles)) == 1:
        return [values[0]] * ploidy
    return [values[index] if index < len(values) else None for index in range(ploidy)]


def _parse_float(value: Optional[str]) -> Optional[float]:
    if value in {None, "", ".", "None", "nan", "NaN"}:
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def _sample_lps_values(
    sample_text: str, format_keys: Sequence[str], motif: str
) -> tuple[list[int], int, int]:
    if sample_text in {"", "."}:
        return [], 0, 0
    fields = dict(zip(format_keys, sample_text.split(":")))
    genotype_alleles = _called_genotype_alleles(fields.get("GT", "."))
    if not genotype_alleles:
        return [], 0, 0

    ds_values = _align_values_to_genotype(
        _split_allele_values(fields.get("DS", ".")), genotype_alleles
    )
    cn_values = _align_values_to_genotype(
        _split_allele_values(fields.get("CN", ".")), genotype_alleles
    )
    results: list[int] = []
    missing = 0
    for ds_value, cn_value in zip(ds_values, cn_values):
        lps = longest_target_lps(ds_value, motif) if ds_value is not None else None
        if lps is None:
            copy_number = _parse_float(cn_value)
            # ATaRVa represents a called complete deletion with no decomposition.
            if copy_number == 0:
                lps = 0
        if lps is None:
            missing += 1
        else:
            results.append(lps)
    return results, len(genotype_alleles), missing


def _format_decimal(value: float) -> str:
    if not math.isfinite(float(value)):
        return "."
    if float(value).is_integer():
        return str(int(value))
    return f"{float(value):.2f}"


def _allele_histogram(values: Sequence[int]) -> str:
    counts = Counter(values)
    return ",".join(f"{allele}x:{counts[allele]}" for allele in sorted(counts)) or "."


def _summary_fields(values: Sequence[int]) -> list[str]:
    if not values:
        return [".", ".", ".", ".", "0", "0", ".", ".", "."]
    counts = Counter(values)
    maximum_count = max(counts.values())
    mode = min(allele for allele, count in counts.items() if count == maximum_count)
    array = np.asarray(values, dtype=float)
    return [
        str(min(values)),
        _format_decimal(float(np.median(array))),
        _format_decimal(float(np.percentile(array, 99))),
        str(max(values)),
        str(len(counts)),
        str(len(values)),
        str(mode),
        _format_decimal(float(np.std(array, ddof=0))),
        _allele_histogram(values),
    ]


def _iter_vcf(
    handle: Iterable[str],
) -> tuple[set[str], list[str], Iterable[tuple[int, list[str]]]]:
    declared_formats: set[str] = set()
    iterator = enumerate(handle, 1)
    for line_number, line in iterator:
        if line.startswith("##FORMAT=<ID="):
            declared_formats.add(line.split("##FORMAT=<ID=", 1)[1].split(",", 1)[0])
        elif line.startswith("#CHROM\t"):
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) < 10:
                raise ValueError("VCF contains no samples")
            samples = columns[9:]

            def records() -> Iterable[tuple[int, list[str]]]:
                for record_line_number, record_line in iterator:
                    if not record_line.strip() or record_line.startswith("#"):
                        continue
                    yield record_line_number, record_line.rstrip("\r\n").split("\t")

            return declared_formats, samples, records()
    raise ValueError("VCF #CHROM header was not found")


def _selected_sample_indexes(
    samples: Sequence[str], sample_file: Optional[str]
) -> list[int]:
    if sample_file is None:
        return list(range(len(samples)))
    with open(sample_file, encoding="utf-8") as handle:
        requested = {line.strip() for line in handle if line.strip()}
    indexes = [index for index, sample in enumerate(samples) if sample in requested]
    if not indexes:
        raise ValueError(f"No samples from {sample_file} were found in the VCF")
    return indexes


def _effective_format_keys(
    format_text: str, declared_formats: set[str]
) -> tuple[list[str], bool]:
    keys = format_text.split(":")
    repaired = False
    # ATaRVa merge versions around the LPM schema change could emit LPM in
    # FORMAT without adding LPM values or its ##FORMAT declaration.  Removing
    # this phantom key restores the older GT:AL:CN:AR:...:DS:MV alignment.
    if "LPM" in keys and "LPM" not in declared_formats:
        keys.remove("LPM")
        repaired = True
    return keys, repaired


def _parse_locus(
    fields: Sequence[str],
    line_number: int,
    sample_count: int,
    sample_indexes: Sequence[int],
    declared_formats: set[str],
) -> tuple[LocusResult, bool]:
    expected_columns = 9 + sample_count
    normalized = list(fields)
    while len(normalized) > expected_columns and normalized[-1] == "":
        normalized.pop()
    if len(normalized) != expected_columns:
        raise ValueError(
            f"line {line_number}: expected {expected_columns} columns, found {len(normalized)}"
        )
    chrom, pos_text, record_id, ref, _alt, _qual, _filter, info_text, format_text = normalized[:9]
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
    except ValueError as error:
        raise ValueError(f"line {line_number}: invalid INFO/START or INFO/END") from error
    locus_id = info.get("ID") or (record_id if record_id != "." else f"{chrom}-{start}-{end}-{motif}")
    format_keys, repaired = _effective_format_keys(format_text, declared_formats)
    if "GT" not in format_keys or "CN" not in format_keys or "DS" not in format_keys:
        raise ValueError(f"line {line_number}: effective FORMAT must contain GT, CN and DS")

    lps_values: list[int] = []
    called_alleles = 0
    missing_lps = 0
    for sample_index in sample_indexes:
        values, called, missing = _sample_lps_values(
            normalized[9 + sample_index], format_keys, motif
        )
        lps_values.extend(values)
        called_alleles += called
        missing_lps += missing
    return (
        LocusResult(
            chrom,
            start,
            end,
            locus_id,
            motif,
            tuple(lps_values),
            called_alleles,
            missing_lps,
        ),
        repaired,
    )


def _write_result(output: TextIO, result: LocusResult) -> None:
    status = (
        "complete"
        if result.missing_lps_alleles == 0
        else "partial"
        if result.lps_values
        else "unavailable"
    )
    row = [
        result.chrom,
        str(result.start),
        str(result.end),
        result.locus_id,
        result.motif,
        *_summary_fields(result.lps_values),
        str(result.vcf_called_alleles),
        str(result.missing_lps_alleles),
        status,
    ]
    output.write("\t".join(row) + "\n")


def getargs(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input ATaRVa merged VCF")
    parser.add_argument("--out", required=True, help="Output TSV path or prefix")
    parser.add_argument("--samples", help="Optional file containing exact sample names to include")
    parser.add_argument("--numrecords", type=int, help="Only process this many records")
    args = parser.parse_args(argv)
    if args.numrecords is not None and args.numrecords < 0:
        parser.error("--numrecords must be nonnegative")
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

    header = [
        "chrom",
        "start",
        "end",
        "locus_id",
        "motif",
        "lps_min",
        "lps_median",
        "lps_99th_percentile",
        "lps_max",
        "lps_unique_allele_lengths",
        "lps_total_genotyped_alleles",
        "lps_mode",
        "lps_stdev",
        "lps_allele_histogram",
        "vcf_called_alleles",
        "missing_lps_alleles",
        "lps_status",
    ]
    processed = 0
    schema_repairs = 0
    total_called = 0
    total_lps = 0
    total_missing = 0
    start_time = time.time()
    try:
        with open(args.vcf, encoding="utf-8") as input_handle:
            declared_formats, samples, records = _iter_vcf(input_handle)
            sample_indexes = _selected_sample_indexes(samples, args.samples)
            with open(output_path, "w", encoding="utf-8", newline="") as output_handle:
                output_handle.write("\t".join(header) + "\n")
                for line_number, fields in records:
                    if args.numrecords is not None and processed >= args.numrecords:
                        break
                    result, repaired = _parse_locus(
                        fields,
                        line_number,
                        len(samples),
                        sample_indexes,
                        declared_formats,
                    )
                    _write_result(output_handle, result)
                    processed += 1
                    schema_repairs += int(repaired)
                    total_called += result.vcf_called_alleles
                    total_lps += len(result.lps_values)
                    total_missing += result.missing_lps_alleles
                    if processed % 1000 == 0:
                        output_handle.flush()
                        elapsed = time.time() - start_time
                        print(
                            f"Processed {processed:,} records ({elapsed / processed:.5f} sec/record)",
                            end="\r",
                            flush=True,
                        )
    except (IndexError, OSError, ValueError) as error:
        _warn(str(error))
        return 1

    elapsed = time.time() - start_time
    print(
        f"Processed {processed:,} loci and {total_called:,} called VCF alleles in {elapsed:.2f}s. "
        f"Measured {total_lps:,} LPS alleles; {total_missing:,} lacked usable DS. "
        f"Repaired the phantom-LPM FORMAT schema on {schema_repairs:,} records. Output: {output_path}"
    )
    return 0


def run() -> None:
    raise SystemExit(main(getargs()))


if __name__ == "__main__":
    run()
