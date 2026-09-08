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