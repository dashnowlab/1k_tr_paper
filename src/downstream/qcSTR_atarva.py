"""Generate core qcSTR plots and data tables from an ATARVA merged VCF.

The ATARVA adaptation reads repeat copy numbers from FORMAT/CN, the reference
copy number from INFO/REFCN, and the motif from INFO/MOTIF.  It intentionally
does not treat depth or SNP metrics as a generic call-quality score.
"""

from __future__ import annotations

import argparse
import math
import os
import statistics
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, legal
from reportlab.pdfgen import canvas


PAGE_WIDTH, PAGE_HEIGHT = landscape(legal)
LEFT = 72
RIGHT = 50
TOP = 64
BOTTOM = 72
INK = colors.HexColor("#263238")
BLUE = colors.HexColor("#1769AA")
RED = colors.HexColor("#C43D3D")
LIGHT_BLUE = colors.HexColor("#DCEAF5")
GRID = colors.HexColor("#D8DEE3")


@dataclass
class QCData:
    samples: list[str]
    sample_calls: list[int]
    chrom_calls: Counter[str]
    diff_units: Counter[float]
    bias_values: dict[int, list[float]]
    records: int
    allele_calls: int


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
    return result


def _parse_cn(format_keys: Sequence[str], sample_text: str) -> tuple[float, ...]:
    if sample_text in {"", "."}:
        return ()
    fields = dict(zip(format_keys, sample_text.split(":")))
    genotype = fields.get("GT", ".")
    genotype_alleles = genotype.replace("|", "/").split("/")
    if not genotype_alleles or "." in genotype_alleles:
        return ()
    values: list[float] = []
    for value in fields.get("CN", ".").split(","):
        try:
            number = float(value)
        except ValueError:
            continue
        if math.isfinite(number):
            values.append(number)
    return tuple(values)


def _iter_records(handle: Iterable[str]) -> tuple[list[str], Iterable[tuple[int, list[str]]]]:
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


def _select_samples(
    all_samples: Sequence[str], sample_file: Optional[str]
) -> tuple[list[int], list[str]]:
    if sample_file is None:
        return list(range(len(all_samples))), list(all_samples)
    with open(sample_file, encoding="utf-8") as handle:
        requested = {line.strip() for line in handle if line.strip()}
    indexes = [index for index, sample in enumerate(all_samples) if sample in requested]
    if not indexes:
        raise ValueError(f"No samples from {sample_file} were found in the VCF")
    return indexes, [all_samples[index] for index in indexes]


def _collect(args: argparse.Namespace) -> QCData:
    start_time = time.time()
    with open(args.vcf, encoding="utf-8") as handle:
        all_samples, records = _iter_records(handle)
        selected_indexes, selected_samples = _select_samples(all_samples, args.samples)
        sample_calls = [0] * len(selected_samples)
        chrom_calls: Counter[str] = Counter()
        diff_units: Counter[float] = Counter()
        bias_values: dict[int, list[float]] = defaultdict(list)
        processed = 0
        allele_calls = 0
        expected_columns = 9 + len(all_samples)

        for line_number, raw_fields in records:
            if args.numrecords is not None and processed >= args.numrecords:
                break
            fields = list(raw_fields)
            while len(fields) > expected_columns and fields[-1] == "":
                fields.pop()
            if len(fields) != expected_columns:
                raise ValueError(
                    f"line {line_number}: expected {expected_columns} columns, found {len(fields)}"
                )
            chrom, pos, _record_id, ref, _alt, _qual, _filter, info_text, format_text = fields[:9]
            info = _parse_info(info_text)
            motif = info.get("MOTIF", "")
            if not motif:
                raise ValueError(f"line {line_number}: missing INFO/MOTIF")
            if args.period is not None and len(motif) != args.period:
                continue
            try:
                ref_cn = float(info["REFCN"])
            except (KeyError, ValueError) as error:
                raise ValueError(f"line {line_number}: missing or invalid INFO/REFCN") from error
            format_keys = format_text.split(":")
            if "GT" not in format_keys or "CN" not in format_keys:
                raise ValueError(f"line {line_number}: FORMAT must contain GT and CN")

            locus_calls = 0
            motif_length = len(motif)
            reference_bp = ref_cn * motif_length
            reference_bin = math.floor(reference_bp / args.refbias_binsize) * args.refbias_binsize
            for output_index, vcf_index in enumerate(selected_indexes):
                copy_numbers = _parse_cn(format_keys, fields[9 + vcf_index])
                if not copy_numbers:
                    continue
                sample_calls[output_index] += 1
                locus_calls += 1
                for copy_number in copy_numbers:
                    difference_units = copy_number - ref_cn
                    diff_units[difference_units] += 1
                    bias_values[reference_bin].append(difference_units * motif_length)
                    allele_calls += 1
            chrom_calls[chrom] += locus_calls
            processed += 1
            if processed % 1000 == 0:
                elapsed = time.time() - start_time
                print(
                    f"Processed {processed:,} records ({elapsed / processed:.5f} sec/record)",
                    end="\r",
                    flush=True,
                )

    return QCData(
        selected_samples,
        sample_calls,
        chrom_calls,
        diff_units,
        dict(bias_values),
        processed,
        allele_calls,
    )


def _setup_page(pdf: canvas.Canvas, title: str, subtitle: str) -> None:
    pdf.setFillColor(INK)
    pdf.setFont("Helvetica-Bold", 18)
    pdf.drawString(LEFT, PAGE_HEIGHT - 42, title)
    pdf.setFont("Helvetica", 9)
    pdf.setFillColor(colors.HexColor("#58656D"))
    pdf.drawRightString(PAGE_WIDTH - RIGHT, PAGE_HEIGHT - 40, subtitle)


def _draw_axes(
    pdf: canvas.Canvas,
    x0: float,
    y0: float,
    width: float,
    height: float,
    x_label: str,
    y_label: str,
) -> None:
    pdf.setStrokeColor(INK)
    pdf.setLineWidth(1)
    pdf.line(x0, y0, x0 + width, y0)
    pdf.line(x0, y0, x0, y0 + height)
    pdf.setFillColor(INK)
    pdf.setFont("Helvetica", 10)
    pdf.drawCentredString(x0 + width / 2, 24, x_label)
    pdf.saveState()
    pdf.translate(20, y0 + height / 2)
    pdf.rotate(90)
    pdf.drawCentredString(0, 0, y_label)
    pdf.restoreState()


def _draw_y_grid(
    pdf: canvas.Canvas,
    x0: float,
    y0: float,
    width: float,
    height: float,
    minimum: float,
    maximum: float,
    formatter,
    ticks: int = 5,
) -> None:
    pdf.setFont("Helvetica", 8)
    for index in range(ticks + 1):
        fraction = index / ticks
        y = y0 + fraction * height
        value = minimum + fraction * (maximum - minimum)
        pdf.setStrokeColor(GRID)
        pdf.setLineWidth(0.5)
        pdf.line(x0, y, x0 + width, y)
        pdf.setFillColor(INK)
        pdf.drawRightString(x0 - 7, y - 3, formatter(value))


def _short_sample(sample: str) -> str:
    return sample.split("-ONT", 1)[0]


def _sample_call_plot(path: str, data: QCData) -> None:
    pdf = canvas.Canvas(path, pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    per_page = 25
    maximum = max(data.sample_calls) if data.sample_calls else 1
    pages = max(1, math.ceil(len(data.samples) / per_page))
    for page_index in range(pages):
        start = page_index * per_page
        stop = min(len(data.samples), start + per_page)
        page_samples = data.samples[start:stop]
        page_counts = data.sample_calls[start:stop]
        _setup_page(
            pdf,
            "Calls per sample",
            f"{data.records:,} loci | page {page_index + 1} of {pages}",
        )
        x0, y0, width, height = LEFT, 150, PAGE_WIDTH - LEFT - RIGHT, 350
        _draw_axes(pdf, x0, y0, width, height, "Sample (short ID)", "Number of called loci")
        _draw_y_grid(pdf, x0, y0, width, height, 0, maximum, lambda value: f"{value:,.0f}")
        slot = width / max(1, len(page_samples))
        bar_width = max(2, slot * 0.68)
        for index, (sample, count) in enumerate(zip(page_samples, page_counts)):
            x = x0 + index * slot + (slot - bar_width) / 2
            bar_height = 0 if maximum == 0 else count / maximum * height
            pdf.setFillColor(BLUE)
            pdf.rect(x, y0, bar_width, bar_height, fill=1, stroke=0)
            pdf.saveState()
            pdf.translate(x + bar_width / 2 + 2, y0 - 8)
            pdf.rotate(-60)
            pdf.setFillColor(INK)
            pdf.setFont("Helvetica", 7)
            pdf.drawRightString(0, 0, _short_sample(sample))
            pdf.restoreState()
        pdf.showPage()
    pdf.save()


def _chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    stripped = chrom.removeprefix("chr")
    if stripped.isdigit():
        return (0, int(stripped))
    if stripped == "X":
        return (1, 23)
    if stripped == "Y":
        return (1, 24)
    if stripped in {"M", "MT"}:
        return (1, 25)
    return (2, stripped)


def _chrom_call_plot(path: str, data: QCData) -> None:
    chroms = sorted((chrom for chrom, count in data.chrom_calls.items() if count), key=_chrom_sort_key)
    counts = [data.chrom_calls[chrom] for chrom in chroms]
    maximum = max(counts) if counts else 1
    pdf = canvas.Canvas(path, pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    _setup_page(pdf, "Calls per chromosome", f"{sum(counts):,} called genotypes")
    x0, y0, width, height = LEFT, 105, PAGE_WIDTH - LEFT - RIGHT, 395
    _draw_axes(pdf, x0, y0, width, height, "Chromosome", "Number of called genotypes")
    _draw_y_grid(pdf, x0, y0, width, height, 0, maximum, lambda value: f"{value:,.0f}")
    slot = width / max(1, len(chroms))
    for index, (chrom, count) in enumerate(zip(chroms, counts)):
        bar_width = slot * 0.65
        x = x0 + index * slot + (slot - bar_width) / 2
        pdf.setFillColor(BLUE)
        pdf.rect(x, y0, bar_width, count / maximum * height, fill=1, stroke=0)
        pdf.setFillColor(INK)
        pdf.setFont("Helvetica", 7)
        pdf.saveState()
        pdf.translate(x + bar_width / 2 + 2, y0 - 8)
        pdf.rotate(-55)
        pdf.drawRightString(0, 0, chrom)
        pdf.restoreState()
    pdf.showPage()
    pdf.save()


def _diff_histogram_plot(path: str, data: QCData) -> None:
    minimum, maximum = -50, 50
    bins = Counter()
    for difference, count in data.diff_units.items():
        integer = int(round(difference))
        bins[max(minimum, min(maximum, integer))] += count
    log_max = max((math.log10(count + 1) for count in bins.values()), default=1)
    pdf = canvas.Canvas(path, pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    _setup_page(pdf, "Allele difference from reference", f"{data.allele_calls:,} called alleles")
    x0, y0, width, height = LEFT, 80, PAGE_WIDTH - LEFT - RIGHT, 420
    _draw_axes(pdf, x0, y0, width, height, "Difference from reference (repeat units)", "Allele count (log10 scale)")
    _draw_y_grid(pdf, x0, y0, width, height, 0, log_max, lambda value: f"10^{value:.0f}")
    values = list(range(minimum, maximum + 1))
    slot = width / len(values)
    for index, value in enumerate(values):
        count = bins[value]
        if count:
            bar_height = math.log10(count + 1) / log_max * height
            pdf.setFillColor(BLUE if value != 0 else RED)
            pdf.rect(x0 + index * slot, y0, max(1, slot - 0.4), bar_height, fill=1, stroke=0)
    pdf.setFillColor(INK)
    pdf.setFont("Helvetica", 8)
    for value in range(minimum, maximum + 1, 10):
        x = x0 + (value - minimum + 0.5) * slot
        pdf.drawCentredString(x, y0 - 14, str(value))
    pdf.setFont("Helvetica", 7)
    pdf.drawString(x0, 44, "Values beyond +/-50 are included in the endpoint bins.")
    pdf.showPage()
    pdf.save()


def _bias_points(args: argparse.Namespace, data: QCData) -> list[tuple[int, float, int]]:
    points = []
    for reference_bin, values in sorted(data.bias_values.items()):
        if not args.refbias_xrange_min <= reference_bin <= args.refbias_xrange_max:
            continue
        if len(values) < max(1, args.refbias_mingts):
            continue
        metric = statistics.mean(values) if args.refbias_metric == "mean" else statistics.median(values)
        points.append((reference_bin, metric, len(values)))
    return points


def _bias_plot(path: str, args: argparse.Namespace, data: QCData) -> list[tuple[int, float, int]]:
    points = _bias_points(args, data)
    pdf = canvas.Canvas(path, pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    _setup_page(
        pdf,
        "Reference-length bias",
        f"{args.refbias_metric.capitalize()} allele difference | bin size {args.refbias_binsize} bp",
    )
    x0, y0, width, height = LEFT, 80, PAGE_WIDTH - LEFT - RIGHT, 420
    _draw_axes(pdf, x0, y0, width, height, "Reference allele length (bp)", f"{args.refbias_metric.capitalize()} difference from reference (bp)")
    if not points:
        pdf.setFillColor(INK)
        pdf.setFont("Helvetica", 12)
        pdf.drawCentredString(x0 + width / 2, y0 + height / 2, "No bins passed the genotype-count and x-range filters")
        pdf.showPage()
        pdf.save()
        return points
    x_min, x_max = args.refbias_xrange_min, args.refbias_xrange_max
    absolute_max = max(abs(metric) for _, metric, _ in points)
    y_limit = max(1.0, absolute_max * 1.1)
    _draw_y_grid(pdf, x0, y0, width, height, -y_limit, y_limit, lambda value: f"{value:.1f}")
    pdf.setFillColor(INK)
    pdf.setFont("Helvetica", 8)
    for index in range(6):
        fraction = index / 5
        x = x0 + fraction * width
        value = x_min + fraction * (x_max - x_min)
        pdf.drawCentredString(x, y0 - 14, f"{value:g}")
        pdf.setFillColor(RED)
        pdf.drawString(x0 + width + 7, y0 + fraction * height - 3, f"{fraction:.1f}")
        pdf.setFillColor(INK)
    pdf.saveState()
    pdf.setFillColor(RED)
    pdf.translate(PAGE_WIDTH - 13, y0 + height / 2)
    pdf.rotate(90)
    pdf.setFont("Helvetica", 8)
    pdf.drawCentredString(0, 0, "Cumulative allele fraction")
    pdf.restoreState()
    zero_y = y0 + height / 2
    pdf.setStrokeColor(colors.HexColor("#7B858C"))
    pdf.setDash(4, 3)
    pdf.line(x0, zero_y, x0 + width, zero_y)
    pdf.setDash()
    total = sum(count for _, _, count in points)
    cumulative = 0
    previous_metric = None
    previous_cumulative = None
    for reference_bin, metric, count in points:
        x = x0 + (reference_bin - x_min) / (x_max - x_min) * width
        y_metric = y0 + (metric + y_limit) / (2 * y_limit) * height
        cumulative += count
        y_cumulative = y0 + cumulative / total * height
        if previous_metric is not None:
            pdf.setStrokeColor(BLUE)
            pdf.setLineWidth(1.8)
            pdf.line(previous_metric[0], previous_metric[1], x, y_metric)
            pdf.setStrokeColor(RED)
            pdf.setLineWidth(1.2)
            pdf.setDash(3, 2)
            pdf.line(previous_cumulative[0], previous_cumulative[1], x, y_cumulative)
            pdf.setDash()
        pdf.setFillColor(BLUE)
        pdf.circle(x, y_metric, 2.7, fill=1, stroke=0)
        previous_metric = (x, y_metric)
        previous_cumulative = (x, y_cumulative)
    pdf.setFillColor(BLUE)
    pdf.rect(PAGE_WIDTH - 270, PAGE_HEIGHT - 66, 12, 3, fill=1, stroke=0)
    pdf.setFillColor(INK)
    pdf.setFont("Helvetica", 8)
    pdf.drawString(PAGE_WIDTH - 253, PAGE_HEIGHT - 68, f"{args.refbias_metric.capitalize()} difference")
    pdf.setStrokeColor(RED)
    pdf.setDash(3, 2)
    pdf.line(PAGE_WIDTH - 145, PAGE_HEIGHT - 64, PAGE_WIDTH - 133, PAGE_HEIGHT - 64)
    pdf.setDash()
    pdf.setFillColor(INK)
    pdf.drawString(PAGE_WIDTH - 128, PAGE_HEIGHT - 68, "Cumulative allele fraction")
    pdf.setFont("Helvetica", 7)
    pdf.drawRightString(PAGE_WIDTH - RIGHT, 44, "Red dashed line uses a 0-1 scale from the bottom to the top of the plot.")
    pdf.showPage()
    pdf.save()
    return points


def _canonical(value: float) -> str:
    return str(int(value)) if value.is_integer() else f"{value:g}"


def _write_tables(prefix: str, args: argparse.Namespace, data: QCData, bias_points: Sequence[tuple[int, float, int]]) -> None:
    with open(prefix + "-sample-callnum.tab", "w", encoding="utf-8", newline="") as handle:
        handle.write("sample\tcalled_loci\tcall_rate\n")
        for sample, count in zip(data.samples, data.sample_calls):
            rate = count / data.records if data.records else math.nan
            handle.write(f"{sample}\t{count}\t{rate:.6f}\n")
    with open(prefix + "-chrom-callnum.tab", "w", encoding="utf-8", newline="") as handle:
        handle.write("chrom\tcalled_genotypes\n")
        for chrom in sorted(data.chrom_calls, key=_chrom_sort_key):
            handle.write(f"{chrom}\t{data.chrom_calls[chrom]}\n")
    with open(prefix + "-diffref-histogram.tab", "w", encoding="utf-8", newline="") as handle:
        handle.write("difference_repeat_units\tallele_count\n")
        for difference in sorted(data.diff_units):
            handle.write(f"{_canonical(difference)}\t{data.diff_units[difference]}\n")
    with open(prefix + "-diffref-bias.tab", "w", encoding="utf-8", newline="") as handle:
        handle.write(f"reference_length_bin_bp\t{args.refbias_metric}_difference_bp\tallele_count\n")
        for reference_bin, metric, count in bias_points:
            handle.write(f"{reference_bin}\t{metric:.6f}\t{count}\n")


def getargs(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="ATARVA merged VCF")
    parser.add_argument("--out", required=True, help="Output prefix")
    parser.add_argument("--samples", help="File containing samples to include")
    parser.add_argument("--period", type=int, help="Only include motifs of this length")
    parser.add_argument("--quality", action="append", default=[], help="Unsupported: ATARVA has no generic scalar call-quality field")
    parser.add_argument("--quality-ignore-no-call", action="store_true")
    parser.add_argument("--refbias-metric", choices=("mean", "median"), default="mean")
    parser.add_argument("--refbias-mingts", type=int, default=100)
    parser.add_argument("--refbias-xrange-min", type=int, default=0)
    parser.add_argument("--refbias-xrange-max", type=int, default=100)
    parser.add_argument("--refbias-binsize", type=int, default=5)
    parser.add_argument("--numrecords", type=int)
    args = parser.parse_args(argv)
    if args.quality:
        parser.error("quality plots are unavailable: ATARVA does not provide a generic scalar call-quality field")
    if args.refbias_binsize < 1:
        parser.error("--refbias-binsize must be at least 1")
    if args.refbias_mingts < 0:
        parser.error("--refbias-mingts must be nonnegative")
    if args.refbias_xrange_min >= args.refbias_xrange_max:
        parser.error("--refbias-xrange-min must be less than --refbias-xrange-max")
    if args.period is not None and args.period < 1:
        parser.error("--period must be positive")
    if args.numrecords is not None and args.numrecords < 0:
        parser.error("--numrecords must be nonnegative")
    return args


def main(args: argparse.Namespace) -> int:
    if not os.path.isfile(args.vcf):
        _warn(f"{args.vcf} does not exist or is not a file")
        return 1
    output_directory = os.path.dirname(os.path.abspath(args.out))
    if not os.path.isdir(output_directory):
        _warn(f"output directory {output_directory} does not exist")
        return 1
    start_time = time.time()
    try:
        data = _collect(args)
        if data.records == 0:
            raise ValueError("No records passed the requested filters")
        if data.allele_calls == 0:
            raise ValueError("No called alleles were found")
        bias_path = args.out + "-diffref-bias.pdf"
        sample_path = args.out + "-sample-callnum.pdf"
        chrom_path = args.out + "-chrom-callnum.pdf"
        histogram_path = args.out + "-diffref-histogram.pdf"
        bias_points = _bias_plot(bias_path, args, data)
        _sample_call_plot(sample_path, data)
        _chrom_call_plot(chrom_path, data)
        _diff_histogram_plot(histogram_path, data)
        _write_tables(args.out, args, data, bias_points)
    except (OSError, ValueError) as error:
        _warn(str(error))
        return 1
    elapsed = time.time() - start_time
    print(
        f"Processed {data.records:,} records and {data.allele_calls:,} called alleles "
        f"in {elapsed:.2f}s. Wrote 4 PDFs and 4 data tables with prefix {args.out}"
    )
    return 0


def run() -> None:
    raise SystemExit(main(getargs()))


if __name__ == "__main__":
    run()
