main_list=/pl/active/dashnowlab/projects/1k_tr_paper/src/1000g-ont-phased-cram.list
completed_list=/pl/active/dashnowlab/projects/1k_tr_paper/data/1000g-output/completed.list

1. Move completed ATaRVa VCFs

cd /pl/active/dashnowlab/projects/1k_tr_paper/data/1000g-output/

for f in logs/*.err; do
    if grep -q "Concatenation complete." "$f"; then
        base=$(basename "$f" | sed -E 's/_[0-9]+\.err$//')
        vcf="${base}.atarva.vcf"

        if [[ -f "$vcf" ]]; then
            mv "$vcf" completed/
        fi
    fi
done

2. Sort, compress, and index completed VCFs

cd /pl/active/dashnowlab/projects/1k_tr_paper/data/1000g-output/completed/

for f in *.vcf; do
    out="${f}.gz"

    if [[ -s "$out" && -s "${out}.tbi" ]]; then
        echo "Skipping completed: $f"
        continue
    fi

    echo "Processing: $f"

    bcftools sort "$f" -Oz -o "$out" &&
        tabix -f -p vcf "$out"
done


3. Generate remaining BAM/CRAM list

awk '
NR == FNR {
    x = $0
    sub(/^.*\//, "", x)
    sub(/\.atarva\.vcf$/, "", x)
    done[x] = 1
    next
}

{
    x = $0
    sub(/^.*\//, "", x)

    if (!(x in done)) {
        print
    }
}
' \
/pl/active/dashnowlab/projects/1k_tr_paper/data/1000g-output/completed.list \
/pl/active/dashnowlab/projects/1k_tr_paper/src/1000g-ont-phased-cram.list \
> remaining_bams.list