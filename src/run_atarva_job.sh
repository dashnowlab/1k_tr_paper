#!/bin/bash -l
set -eo pipefail

CRAM="${1:?CRAM required}"
OUTPUT_VCF="${2:?Output VCF required}"

REFERENCE="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
CATALOG="/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/benchmark-catalog-v2.atarva.bed.gz"
ATARVA_ENV="/projects/ealiyev@xsede.org/software/anaconda/envs/atarva_0.7.1"

module load miniforge

conda activate "$ATARVA_ENV"

echo "Node: $(hostname)"
echo "Sample: $(basename "$CRAM" .cram)"
echo "ATaRVa: $(which atarva)"
atarva --version

atarva genotype \
    -t "${SLURM_CPUS_PER_TASK:-8}" \
    -f "$REFERENCE" \
    -b "$CRAM" \
    -r "$CATALOG" \
    --format cram \
    -o "$OUTPUT_VCF"