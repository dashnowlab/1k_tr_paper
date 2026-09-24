# ATaRVa Workflow

This directory contains scripts and input lists used to genotype tandem repeats with **ATaRVa** and merge per-sample ATaRVa VCFs into a multisample VCF.

## Workflow

The workflow consists of two main steps:

1. Run ATaRVa genotyping independently for each CRAM.
2. Compress/index the resulting VCFs and merge them into a multisample ATaRVa VCF.

```text
CRAM list
   |
   v
submit_atarva.sh
   |
   v
run_atarva_job.sh
   |
   v
Individual *.atarva.vcf files
   |
   v
merge_atarva.sh
   |
   v
HPRC_ATaRVa_multisample.vcf

 ./submit_atarva.sh cram.list

sbatch \
  -J atarva_merge \
  -p acpu \
  --qos=cpu-normal \
  --time=24:00:00 \
  --cpus-per-task=8 \
  --mem=16G \
  --output=/pl/active/dashnowlab/work/ealiyev/Sandbox/HPRC_atarva_40k_merged/atarva_merge_%j.out \
  --error=/pl/active/dashnowlab/work/ealiyev/Sandbox/HPRC_atarva_40k_merged/atarva_merge_%j.err \
  /pl/active/dashnowlab/work/ealiyev/Sandbox/HPRC_atarva_40k_merged/merge.sh
