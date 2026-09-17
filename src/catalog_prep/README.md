# TR Catalog Overlap Merge
 
Merges overlapping TR regions and picks one consensus record per cluster.
 
## What it does
 
1. Cluster overlapping regions with `bioframe.cluster()`.
2. Rank loci in each cluster with `rank_key` (source priority → motif purity →
   interruptions → allele frequency → length → `_orig_id` tiebreak).
3. Keep the top-ranked locus per cluster, using the cluster's **unioned**
   coordinates as `START`/`END`.

## Usage
 
```bash
python filter_duplicates <catalog_path> <output_dir> --format 
```
 
| Flag | Description |
|---|---|
| `cat_path` | Path to the input TR catalog |
| `out_dir` | Directory to write output files to |
| `--format` | Output format (currently just artarva) |
| `--stat_only` | Only compute/report overlap stats, skip writing/performing any merge actions |
| `--rem_dupes` | Remove exact-duplicate regions (`ReferenceRegion` dedup) |
| `--no_merge` | Skip overlap clustering/merging |


## Output files
 
Normal runs produce three files:
 
- **Format-specified file** — the standard output (e.g. `#CHROM`/`START`/`END`/`MOTIF`/`MOTIF_LEN`) for downstream tools.
- **Merged TSV** — the merged catalog with the added `merged_from` / `merged_motifs_other` columns.
- **Annotated input catalog** — a copy of the original input with `_orig_id` added to each row, so any ID in `merged_from` can be traced back to its original record.

## Output columns
 
- `merged_from` — original locus IDs merged into this record
- `merged_motifs_other` — motifs of merged loci not chosen as consensus
## Post-merge QC (bedtools)
 
**No overlaps remain:**
```bash
sort -k1,1 -k2,2n merged.bed > sorted.bed
bedtools merge -i sorted.bed -c 1 -o count > mrgd_check.bed
awk '{sum += $NF > 1} END {print sum}' mrgd_check.bed   # expect 0
```
 
**No coverage lost:**
```bash
cut -f1-3 original.bed > og_coords.bed
cut -f1-3 merged.bed > mrgd_coords.bed
bedtools subtract -a og_coords.bed -b mrgd_coords.bed | wc -l   # expect 0
```

## Stats from TRExplorer 
Original row count       : 5599658
Bedtools overlap OG count: 1355802 (952980 when accounting for 1 kept per cluster)
Post-merge row count     : 4646678
Post-merge overlap count : 0
Bedtools subtract count  : 0

 
## Requirements
 
`bioframe`, `pandas`, `bedtools` (QC only)
