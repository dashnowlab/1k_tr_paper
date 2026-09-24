#!/bin/bash
set -eo pipefail

CRAM_LIST="${1:?Usage: $0 cram.list}"

OUT="/pl/active/dashnowlab/projects/1k_tr_paper/data/1000g-output/"
LOG="${OUT}/logs"
RUN_SCRIPT="$(cd "$(dirname "$0")" && pwd)/run_atarva_job.sh"

echo "CRAM list: $CRAM_LIST"
echo "RUN_SCRIPT: $RUN_SCRIPT"

mkdir -p "$OUT" "$LOG"

while IFS= read -r CRAM; do
    [[ -z "$CRAM" || "$CRAM" == \#* ]] && continue

    SAMPLE=$(basename "$CRAM" .cram)
    VCF="${OUT}/${SAMPLE}.atarva.vcf"

    if [[ -s "$VCF" ]]; then
        echo "$SAMPLE: exists, skipping"
        continue
    fi

    JOBID=$(sbatch --parsable \
        -J "atarva_${SAMPLE}" \
        -p acpu \
        --qos=cpu-normal \
        --cpus-per-task=8 \
        --mem=16G \
        --time=23:59:00 \
        --output="${LOG}/${SAMPLE}_%j.out" \
        --error="${LOG}/${SAMPLE}_%j.err" \
        "$RUN_SCRIPT" "$CRAM" "$VCF")

    echo "$SAMPLE: submitted $JOBID"

done < "$CRAM_LIST"