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

sbatch merge_atarva.sh completed_vcf.list /pl/active/dashnowlab/projects/1k_tr_paper/data/1000g-output/merged/
