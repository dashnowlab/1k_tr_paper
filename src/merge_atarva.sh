#!/bin/bash -l
#SBATCH --job-name=atarva_merge
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/pl/active/dashnowlab/work/ealiyev/Sandbox/HPRC_atarva_40k_merged/atarva_merge.log
#SBATCH --error=/pl/active/dashnowlab/work/ealiyev/Sandbox/HPRC_atarva_40k_merged/atarva_merge.log

set -euo pipefail

VCF_LIST="${1:?VCF list required}"
OUTPUT_DIR="${2:?Output directory required}"

REFERENCE="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
CATALOG="/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/benchmark-catalog-v2.atarva.bed.gz"
ATARVA_ENV="/projects/ealiyev@xsede.org/software/anaconda/envs/atarva_0.7.1"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

module load miniforge
module load htslib/1.16

conda activate "$ATARVA_ENV"

mkdir -p "$OUTPUT_DIR/compressed"

COMPRESSED_LIST="$OUTPUT_DIR/atarva_vcfs_compressed.txt"
OUTPUT_VCF="$OUTPUT_DIR/HPRC_ATaRVa_multisample.vcf"

> "$COMPRESSED_LIST"

while read -r VCF
do
    [[ -z "$VCF" ]] && continue

    SAMPLE=$(basename "$VCF" .atarva.vcf)
    OUT="$OUTPUT_DIR/compressed/${SAMPLE}.atarva.vcf.gz"

    if [[ ! -s "$OUT" ]]; then
        bgzip -@ "$THREADS" -c "$VCF" > "$OUT"
    fi

    if [[ ! -s "${OUT}.tbi" ]]; then
        tabix -p vcf "$OUT"
    fi

    echo "$OUT" >> "$COMPRESSED_LIST"

done < "$VCF_LIST"

atarva merge \
    -i "$COMPRESSED_LIST" \
    -r "$CATALOG" \
    -f "$REFERENCE" \
    -o "$OUTPUT_VCF" \
    -t "$THREADS"

echo "Done: $OUTPUT_VCF"