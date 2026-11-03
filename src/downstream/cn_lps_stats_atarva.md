# Scientific definitions and validation

`cn_lps_stats_atarva.py` summarizes two different allele measurements:

- **CN** is the total repeat copy number reported in ATaRVa `FORMAT/CN`.
- **LPS** is the longest uninterrupted target-motif run reconstructed from
  ATaRVa `FORMAT/DS`. Cyclic motif rotations and reverse-complement rotations
  are treated as equivalent. This is an ATaRVa-DS-derived measurement, not the
  TRGT-LPS HMM score.

Missing DS is not imputed, except that a called allele with CN zero is assigned
LPS zero. Consequently, CN and LPS denominators can differ. The table reports
those denominators and missing counts explicitly.

## Distribution definitions

All CN and LPS distribution statistics use observed alleles, rather than
samples, as their denominator.

- Percentiles use NumPy's default linear interpolation (Hyndman-Fan type 7).
- Variance and standard deviation describe the observed allele population and
  therefore use `ddof=0` (division by N).
- Coefficient of variation is the unitless population standard deviation
  divided by the arithmetic mean. It is missing (`.`) when the mean is zero.
- A tied mode is resolved to the smallest allele length.
- Heterozygosity is `1 - sum(p_i^2)`.
- Entropy is Shannon entropy in bits: `-sum(p_i * log2(p_i))`.
- Effective allele count is the inverse Simpson concentration:
  `1 / sum(p_i^2)`.
- `singleton_allele_count` counts distinct allele-length categories observed
  exactly once.
- `rare_allele_count` counts allele observations belonging to categories with
  frequency strictly below `--rare-allele-threshold` (default 0.01).
- CN/LPS ratios exclude CN-zero alleles because their ratio is undefined. They
  remain included in the paired-allele and purity counts.

Values are written with at most six digits after the decimal point. Tests use a
tolerance consistent with that output rounding.

## HWE definition and interpretation

`hwep` follows statSTR/TRTools semantics: it is a two-sided binomial test of
the observed number of heterozygous diploid genotypes against the expected
heterozygosity computed from the observed CN allele frequencies. Only complete
diploid CN genotypes contribute. `hwe_sample_count` reports that denominator.

This is not a full multi-allelic exact genotype-table test. HWE is also
sensitive to population mixture, relatedness, selection, technical artifacts,
and sex-chromosome ploidy. For biological interpretation, calculate it within
appropriate ancestry/population strata and treat loci with small sample counts
carefully.

## QC definitions

- Call rate uses selected samples as the denominator.
- CN and LPS coverage use called alleles as the denominator.
- Low depth is `DP < --low-depth-threshold` (default 10).
- Allele balance is `SD[0] / (SD[0] + SD[1])` when at least two support values
  are present and their sum is positive. The default outlier interval is
  `[0.2, 0.8]`; values exactly on either boundary are not outliers.
- `lps_exceeds_cn_count` is retained as a diagnostic inconsistency count. An
  LPS longer than total CN is biologically unexpected and should be inspected.
