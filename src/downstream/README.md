# ATaRVa CN and LPS statistics

`cn_lps_stats_atarva.py` reads an ATaRVa merged VCF once and writes one
extended per-locus table containing CN statistics, DS-derived LPS statistics,
CN/LPS comparisons, coefficient of variation for both measurements, genotype
summaries, HWE, depth, and allele-balance QC.

## Requirements

- Python 3.10 or newer
- NumPy

From PowerShell, install the Python dependencies if they are not already
available:

```powershell
python -m pip install numpy
```

## Run from the repository root

Open PowerShell and change to the project directory:

```powershell
Set-Location "C:\Users\aelbay\OneDrive - The University of Colorado Denver\Desktop\Work\Projects\1k_tr_paper"
```

Run the complete CN/LPS analysis:

```powershell
python .\src\downstream\cn_lps_stats_atarva.py `
  --vcf ".\data\1000g_output\1000g_ont_test.vcf" `
  --out ".\data\1000g_output\1000g_ont_test.cn_lps_extended.tab"
```

The resulting table is:

```text
data\1000g_output\1000g_ont_test.cn_lps_extended.tab
```

The `data/1000g_output/` directory is ignored by Git, so generated tables and
the large input VCF will not be committed accidentally.

## Run from the downstream directory

The script can also be run while the current directory is `src/downstream`:

```powershell
Set-Location .\src\downstream

python .\cn_lps_stats_atarva.py `
  --vcf "..\..\data\1000g_output\1000g_ont_test.vcf" `
  --out "..\..\data\1000g_output\1000g_ont_test.cn_lps_extended.tab"
```

Keep `lpsSTR_atarva.py` in the same directory as
`cn_lps_stats_atarva.py`, because the combined script imports its shared VCF
and LPS parsing functions.

## Command-line options

Show the complete help message:

```powershell
python .\src\downstream\cn_lps_stats_atarva.py --help
```

Important options:

| Option | Default | Description |
|---|---:|---|
| `--vcf PATH` | required | Input ATaRVa merged VCF. |
| `--out PATH` | required | Output `.tab` or `.tsv` file. If neither extension is supplied, `.tab` is added. |
| `--samples PATH` | all samples | Text file containing one exact VCF sample name per line. |
| `--numrecords N` | all loci | Process only the first `N` VCF records. Useful for a smoke test. |
| `--rare-allele-threshold X` | `0.01` | An allele category is rare when its frequency is strictly below `X`. |
| `--low-depth-threshold X` | `10` | A measured sample is low depth when `DP < X`. |
| `--allele-balance-min X` | `0.2` | Lower bound of the accepted allele-balance interval. |
| `--allele-balance-max X` | `0.8` | Upper bound of the accepted allele-balance interval. |

## Quick smoke test

Process only the first ten loci and write a temporary output table:

```powershell
python .\src\downstream\cn_lps_stats_atarva.py `
  --vcf ".\data\1000g_output\1000g_ont_test.vcf" `
  --out ".\data\1000g_output\cn_lps_smoke_test.tab" `
  --numrecords 10
```

The terminal summary reports the number of processed loci, called alleles, CN
measurements, LPS measurements, and repaired phantom-LPM schemas. A successful
run exits with code zero.

## Input assumptions

The input must be an ATaRVa VCF with sample-level `GT`, `CN`, `DS`, `DP`, and
`SD` FORMAT values. The script repairs the known merged-VCF case where `LPM`
appears in the record FORMAT but has no header declaration or sample value.

LPS is calculated from the ATaRVa motif-decomposition field (`DS`). It is not a
TRGT-LPS HMM score. Missing DS is retained as missing rather than imputed,
except that a called complete deletion with CN zero receives LPS zero.

See [cn_lps_stats_atarva.md](cn_lps_stats_atarva.md) for formulas, denominators,
HWE interpretation, and scientific limitations.

## Repair records joined by missing newlines

If a malformed merged VCF contains two or more records joined on one physical
line, scan the entire file without writing output:

```powershell
python .\src\downstream\repair_joined_vcf_records.py `
  --input ".\data\1000g_output\1000g_ont_test.vcf" `
  --check-only
```

Create a validated repaired copy:

```powershell
python .\src\downstream\repair_joined_vcf_records.py `
  --input ".\data\1000g_output\1000g_ont_test.vcf" `
  --output ".\data\1000g_output\1000g_ont_test.fixed.vcf"
```

The repair utility derives the expected width from the VCF header, fixes every
unambiguous missing-newline join, normalizes trailing empty fields, and validates
all records. It writes through a temporary file, never overwrites an existing
output, and removes the temporary output if any record cannot be repaired
safely.
