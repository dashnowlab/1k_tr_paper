# Data

This directory contains input datasets, metadata, reference catalogues,
and processed datasets used in the analyses.

## Documentation

- [Data Description](data_description.md)
- [Input Datasets](input_datasets.md)
- [Metadata](metadata.md)

## TR explorer catalog

Latest release of [TRexplorer](https://trexplorer.broadinstitute.org/) catalog as of 9/16/2026 is downloaded. Converted to atarva input format using `../src/python/trx2atarva.py`.

```bash
wget -c https://github.com/broadinstitute/trexplorer-catalog/releases/download/v2.0/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz
python trx2atarva.py TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz
bgzip TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.atarva.bed
tabix -p bed TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.atarva.bed.gz
```