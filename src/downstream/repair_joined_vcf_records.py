"""Safely repair VCF records joined by missing newline characters.

The utility derives the expected record width from the #CHROM header. A normal
record is written with exactly that many tab-delimited fields (harmless trailing
empty fields are removed). When two or more records were fused, the last sample
field of one record contains the chromosome name of the next record as a
suffix. The utility splits only when the field count, chromosome suffix, POS,
and reconstructed VCF records are all structurally consistent.

Output is written to a temporary file and moved into place only after the whole
input passes validation. Existing output files are never overwritten.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, TextIO


CONTIG_HEADER = re.compile(r"^##contig=<ID=([^,>]+)")
COMMON_CHROMOSOME_SUFFIX = re.compile(r"((?:chr)?(?:[0-9]+|X|Y|M|MT))$")


@dataclass
class RepairSummary:
    header_columns: int = 0
    samples: int = 0
    physical_data_lines: int = 0
    logical_records: int = 0
    joined_lines: int = 0
    recovered_records: int = 0
    trailing_empty_lines: int = 0
    duplicate_chrom_headers: int = 0


def _open_input(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _open_output(path: Path, compressed: bool) -> TextIO:
    if compressed:
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return path.open("w", encoding="utf-8", newline="")


def _trim_trailing_empty_fields(fields: list[str], expected: int) -> int:
    removed = 0
    while len(fields) > expected and fields[-1] == "":
        fields.pop()
        removed += 1
    return removed


def _validate_record(fields: list[str], expected: int, line_number: int) -> None:
    if len(fields) != expected:
        raise ValueError(
            f"physical line {line_number}: reconstructed record has "
            f"{len(fields)} columns; expected {expected}"
        )
    if not fields[0] or fields[0] == ".":
        raise ValueError(f"physical line {line_number}: invalid CHROM {fields[0]!r}")
    try:
        position = int(fields[1])
    except ValueError as error:
        raise ValueError(
            f"physical line {line_number}: invalid POS {fields[1]!r}"
        ) from error
    if position <= 0:
        raise ValueError(f"physical line {line_number}: POS must be positive")
    if not fields[3] or fields[3] == ".":
        raise ValueError(f"physical line {line_number}: invalid REF {fields[3]!r}")
    if expected > 8 and not fields[8]:
        raise ValueError(f"physical line {line_number}: empty FORMAT field")


def _chromosome_suffix(
    junction: str,
    known_contigs: set[str],
    next_position: str,
    line_number: int,
) -> tuple[str, str]:
    try:
        if int(next_position) <= 0:
            raise ValueError
    except ValueError as error:
        raise ValueError(
            f"physical line {line_number}: the field after a possible join is "
            f"not a valid POS ({next_position!r})"
        ) from error

    candidates = [
        contig
        for contig in known_contigs
        if junction.endswith(contig) and len(junction) > len(contig)
    ]
    if candidates:
        chromosome = max(candidates, key=len)
    else:
        match = COMMON_CHROMOSOME_SUFFIX.search(junction)
        if not match or match.start() == 0:
            raise ValueError(
                f"physical line {line_number}: could not unambiguously find the "
                "next chromosome suffix in the joined sample field"
            )
        chromosome = match.group(1)

    sample_value = junction[:-len(chromosome)]
    if not sample_value:
        raise ValueError(
            f"physical line {line_number}: splitting the chromosome suffix "
            "would leave an empty final sample field"
        )
    return sample_value, chromosome


def _split_joined_records(
    fields: list[str],
    expected: int,
    known_contigs: set[str],
    line_number: int,
) -> list[list[str]]:
    if expected <= 1 or (len(fields) - 1) % (expected - 1) != 0:
        raise ValueError(
            f"physical line {line_number}: found {len(fields)} nonempty columns; "
            f"this is not N joined {expected}-column records"
        )

    record_count = (len(fields) - 1) // (expected - 1)
    if record_count < 2:
        raise ValueError(
            f"physical line {line_number}: unexpected field count {len(fields)}"
        )

    records: list[list[str]] = []
    remaining = list(fields)
    while len(remaining) > expected:
        if len(remaining) < 2 * expected - 1:
            raise ValueError(
                f"physical line {line_number}: incomplete joined VCF record"
            )
        junction = remaining[expected - 1]
        final_sample, next_chromosome = _chromosome_suffix(
            junction, known_contigs, remaining[expected], line_number
        )
        record = remaining[:expected]
        record[-1] = final_sample
        _validate_record(record, expected, line_number)
        records.append(record)
        known_contigs.add(record[0])
        known_contigs.add(next_chromosome)
        remaining = [next_chromosome] + remaining[expected:]

    _validate_record(remaining, expected, line_number)
    records.append(remaining)
    if len(records) != record_count:
        raise ValueError(
            f"physical line {line_number}: expected to reconstruct {record_count} "
            f"records but reconstructed {len(records)}"
        )
    return records


def repair_vcf(
    input_path: Path,
    output_path: Optional[Path],
    check_only: bool = False,
    quiet: bool = False,
) -> RepairSummary:
    if not input_path.is_file():
        raise ValueError(f"input does not exist or is not a file: {input_path}")
    if not check_only:
        if output_path is None:
            raise ValueError("--output is required unless --check-only is used")
        if output_path.exists():
            raise ValueError(f"output already exists; refusing to overwrite: {output_path}")
        if not output_path.parent.is_dir():
            raise ValueError(f"output directory does not exist: {output_path.parent}")

    summary = RepairSummary()
    known_contigs: set[str] = set()
    chrom_header: Optional[list[str]] = None
    saw_data = False
    temporary_path: Optional[Path] = None
    output_handle: Optional[TextIO] = None

    try:
        if not check_only and output_path is not None:
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=f".{output_path.name}.", suffix=".tmp", dir=output_path.parent
            )
            os.close(descriptor)
            temporary_path = Path(temporary_name)
            output_handle = _open_output(
                temporary_path, output_path.name.lower().endswith(".gz")
            )

        with _open_input(input_path) as input_handle:
            for line_number, line in enumerate(input_handle, 1):
                if line.startswith("##contig=<ID="):
                    match = CONTIG_HEADER.match(line)
                    if match:
                        known_contigs.add(match.group(1))

                if line.startswith("#CHROM\t"):
                    fields = line.rstrip("\r\n").split("\t")
                    if chrom_header is None:
                        chrom_header = fields
                        summary.header_columns = len(fields)
                        summary.samples = len(fields) - 9
                        if summary.header_columns < 10:
                            raise ValueError("VCF header contains no sample columns")
                        if output_handle is not None:
                            output_handle.write("\t".join(fields) + "\n")
                    elif fields == chrom_header:
                        summary.duplicate_chrom_headers += 1
                        if not quiet:
                            print(
                                f"Skipping duplicate identical #CHROM header at "
                                f"physical line {line_number}",
                                file=sys.stderr,
                            )
                    else:
                        raise ValueError(
                            f"physical line {line_number}: encountered a second "
                            "#CHROM header with a different sample schema"
                        )
                    continue

                if line.startswith("#"):
                    if saw_data:
                        raise ValueError(
                            f"physical line {line_number}: header/comment line occurs "
                            "after data records"
                        )
                    if output_handle is not None:
                        output_handle.write(line)
                    continue

                if not line.strip():
                    continue
                if chrom_header is None:
                    raise ValueError(
                        f"physical line {line_number}: data occurs before #CHROM header"
                    )

                saw_data = True
                summary.physical_data_lines += 1
                fields = line.rstrip("\r\n").split("\t")
                removed = _trim_trailing_empty_fields(fields, summary.header_columns)
                if removed:
                    summary.trailing_empty_lines += 1

                if len(fields) == summary.header_columns:
                    records = [fields]
                    _validate_record(fields, summary.header_columns, line_number)
                elif len(fields) > summary.header_columns:
                    records = _split_joined_records(
                        fields, summary.header_columns, known_contigs, line_number
                    )
                    summary.joined_lines += 1
                    summary.recovered_records += len(records) - 1
                    if not quiet:
                        coordinates = ", ".join(
                            f"{record[0]}:{record[1]}" for record in records
                        )
                        print(
                            f"Repaired physical line {line_number}: split into "
                            f"{len(records)} records ({coordinates})",
                            file=sys.stderr,
                        )
                else:
                    raise ValueError(
                        f"physical line {line_number}: found {len(fields)} columns; "
                        f"expected {summary.header_columns}"
                    )

                for record in records:
                    known_contigs.add(record[0])
                    summary.logical_records += 1
                    if output_handle is not None:
                        output_handle.write("\t".join(record) + "\n")

        if chrom_header is None:
            raise ValueError("VCF #CHROM header was not found")
        if output_handle is not None:
            output_handle.close()
            output_handle = None
        if not check_only and output_path is not None and temporary_path is not None:
            os.replace(temporary_path, output_path)
            temporary_path = None
    except Exception:
        if output_handle is not None:
            output_handle.close()
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()
        raise

    return summary


def getargs(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Input VCF or VCF.GZ")
    parser.add_argument("--output", type=Path, help="New repaired VCF or VCF.GZ")
    parser.add_argument(
        "--check-only", action="store_true",
        help="Scan and report repairs without writing an output file",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress per-line repair messages")
    args = parser.parse_args(argv)
    if args.check_only and args.output is not None:
        parser.error("--output cannot be combined with --check-only")
    if not args.check_only and args.output is None:
        parser.error("--output is required unless --check-only is used")
    return args


def main(args: argparse.Namespace) -> int:
    try:
        summary = repair_vcf(args.input, args.output, args.check_only, args.quiet)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print(f"Header columns: {summary.header_columns}")
    print(f"Samples: {summary.samples}")
    print(f"Physical data lines: {summary.physical_data_lines}")
    print(f"Validated logical records: {summary.logical_records}")
    print(f"Joined physical lines repaired: {summary.joined_lines}")
    print(f"Additional records recovered: {summary.recovered_records}")
    print(f"Lines with trailing empty fields normalized: {summary.trailing_empty_lines}")
    print(f"Duplicate identical #CHROM headers removed: {summary.duplicate_chrom_headers}")
    if args.check_only:
        print("Check-only scan completed; no output was written.")
    else:
        print(f"Repaired VCF: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(getargs()))
