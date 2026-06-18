# pgscorer

An R package for computing Polygenic Risk Scores (PRS) from [PGS Catalog](https://www.pgscatalog.org/) scoring files against a BGZF-compressed VCF file.

The package unions SNP positions across all models, queries the VCF once per position in batches, and scores every model from the shared dosage table. On **Windows** batches run in-process; on **Linux/macOS** each batch runs in an isolated `Rscript` subprocess to keep peak memory bounded to one batch at a time.

## Installation

```r
remotes::install_github("USCbiostats/pgscorer")
```

### Dependencies

- [tabixr](https://github.com/USCbiostats/tabixr) (Bioconductor) — BGZF/Tabix VCF queries
- [data.table](https://cran.r-project.org/package=data.table) — fast PGS file reading

## Usage

```r
library(pgscorer)

results <- compute_prs(
  vcf_path   = "path/to/genotypes.vcf.gz",
  pgs_files  = c("PGS002164_hmPOS_GRCh37.txt.gz",
                 "PGS002863_hmPOS_GRCh37.txt.gz"),
  batch_size = 10000L,
  output_dir = "."       # writes one .rds per model; NULL to skip
)
```

If `pgs_files` is omitted, all files matching `^PGS.*\.txt\.gz$` in the current directory (or `pgs_dir`) are used automatically.

`compute_prs()` returns a named list of numeric vectors — one element per model, named by sample identifier.

### Arguments

| Argument | Default | Description |
|---|---|---|
| `vcf_path` | — | Path to `.vcf.gz` file (`.tbi` index must exist alongside it) |
| `pgs_files` | `NULL` | Character vector of PGS Catalog scoring file paths; `NULL` = auto-discover in `pgs_dir` |
| `pgs_dir` | `"."` | Directory searched for PGS files when `pgs_files` is `NULL` |
| `chrom` | `NULL` | Chromosome label (e.g. `"21"`); `NULL` = auto-detect from the `.tbi` index |
| `batch_size` | `10000L` | Positions per VCF query batch |
| `output_dir` | `"."` | Directory for output `.rds` files; `NULL` = do not save |
| `verbose` | `TRUE` | Print progress and summary |

### Output files

When `output_dir` is not `NULL`, one RDS file is written per model:

```
pgs_chr<chrom>_<model>_prs_combined.rds
```

Each file contains a named numeric vector of PRS values indexed by sample ID.

## PGS Catalog file format

Scoring files should be downloaded directly from [pgscatalog.org](https://www.pgscatalog.org/) in the standard tab-delimited `.txt.gz` format. Harmonised position columns (`hm_chr`, `hm_pos`) are preferred over raw columns (`chr_name`, `chr_position`) when both are present.

## VCF requirements

- BGZF-compressed (`.vcf.gz`) with a Tabix index (`.vcf.gz.tbi`)
- Must contain a `DS` (dosage) field in the `FORMAT` column

## Allele effect formula

For SNPs where the effect allele matches `REF`:

```
PRS contribution = weight × (2 − dosage)
```

For SNPs where the effect allele matches `ALT`:

```
PRS contribution = weight × dosage
```

Rows where the effect allele matches neither `REF` nor `ALT` are excluded with a warning.

## License

MIT © John Morrison
