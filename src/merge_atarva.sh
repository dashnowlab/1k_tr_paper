#!/bin/bash -l
#SBATCH --job-name=atarva_merge
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/pl/active/dashnowlab/projects/1k_tr_paper/data/1000g-output/merged//atarva_merge.log
#SBATCH --error=/pl/active/dashnowlab/projects/1k_tr_paper/data/1000g-output/merged//atarva_merge.log

set -euo pipefail

VCF_LIST="${1:?VCF list required}"
OUTPUT_DIR="${2:?Output directory required}"

REFERENCE="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
CATALOG="/pl/active/dashnowlab/projects/1k_tr_paper/data/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.atarva.bed.gz"
ATARVA_ENV="/projects/ealiyev@xsede.org/software/anaconda/envs/atarva-0.9.1"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

module load miniforge
conda activate "$ATARVA_ENV"

mkdir -p "$OUTPUT_DIR"

OUTPUT_VCF="$OUTPUT_DIR/1000g_ont_merged.vcf"

atarva merge \
    -i "$VCF_LIST" \
    -r "$CATALOG" \
    -f "$REFERENCE" \
    -o "$OUTPUT_VCF" \
    -t "$THREADS"

echo "Done: $OUTPUT_VCF"