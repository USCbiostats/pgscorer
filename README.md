# pgscorer

An R package for computing Polygenic Risk Scores (PRS) from [PGS Catalog](https://www.pgscatalog.org/) scoring files against per-chromosome genotype files, in either BGZF/Tabix VCF or [BinaryDosage](https://cran.r-project.org/package=BinaryDosage) format.

Each PGS model may span any number of chromosomes. `compute_prs()` reads all requested models, determines which chromosomes they need, and — for each chromosome — unions positions across every model that touches it, queries the matching genotype file once, and scores all of those models from the shared dosage table. On **Windows**, VCF batches run in-process; on **Linux/macOS** each VCF batch runs in an isolated `Rscript` subprocess to keep peak memory bounded to one batch at a time. BinaryDosage files are always queried in-process (random-access reads, no large sequential scan to isolate).

## Installation

```r
remotes::install_github("USCbiostats/pgscorer")
```

### Dependencies

- [tabixr](https://github.com/USCbiostats/tabixr) (Bioconductor) — BGZF/Tabix VCF queries
- [data.table](https://cran.r-project.org/package=data.table) — fast PGS file reading
- [BinaryDosage](https://cran.r-project.org/package=BinaryDosage) (CRAN) — only required when scoring against `.bdose` files

## Usage

```r
library(pgscorer)

results <- compute_prs(
  geno_dir   = "path/to/genotypes",   # directory of indexed VCFs (.vcf.gz + .tbi) and/or BinaryDosage files (.bdose + .bdi)
  format     = NULL,                  # NULL = autodetect; "vcf" or "bdose" to force
  pgs_files  = c("PGS002164_hmPOS_GRCh37.txt.gz",
                 "PGS002863_hmPOS_GRCh37.txt.gz"),
  batch_size = 10000L,
  output_dir = "."       # writes one .rds per model; NULL to skip
)
```

If `pgs_files` is omitted, all files matching `^PGS.*\.txt\.gz$` in the current directory (or `pgs_dir`) are used automatically.

`compute_prs()` returns a named list (invisibly), one element per model, each itself a list:

```r
results$PGS002164$prs                  # named numeric vector of PRS values, one per sample
results$PGS002164$unmatched_by_chr     # named integer vector: unmatched SNP count per chromosome
results$PGS002164$unmatched_rsIDs      # named list: unmatched rsIDs per chromosome
results$PGS002164$excluded_chr_counts  # named integer vector: allele-mismatch exclusions per chromosome
```

### Arguments

| Argument | Default | Description |
|---|---|---|
| `geno_dir` | `"."` | Directory containing the genotype files: BGZF VCFs (`*.vcf.gz` with `*.vcf.gz.tbi` alongside) and/or BinaryDosage files (`*.bdose` with `*.bdose.bdi` alongside). File names are arbitrary |
| `format` | `NULL` | `"vcf"` or `"bdose"` to force the genotype file type; `NULL` = autodetect from what's in `geno_dir`. If both types are present, BinaryDosage is used and a message is printed |
| `pgs_files` | `NULL` | Character vector of PGS Catalog scoring file paths; `NULL` = auto-discover in `pgs_dir` |
| `pgs_dir` | `"."` | Directory searched for PGS files when `pgs_files` is `NULL` |
| `batch_size` | `10000L` | Positions per VCF query batch (ignored for BinaryDosage input) |
| `output_dir` | `"."` | Directory for output `.rds` files; `NULL` = do not save |
| `verbose` | `TRUE` | Print progress and summary |

A model may reference a chromosome for which no matching genotype file exists in `geno_dir`. Rather than erroring, every SNP on that chromosome is reported as unmatched (in `unmatched_by_chr` / `unmatched_rsIDs`) and a warning is printed.

### Output files

When `output_dir` is not `NULL`, one RDS file is written per model:

```
pgs_<model>_prs.rds
```

Each file contains that model's full result list (`prs`, `unmatched_by_chr`, `unmatched_rsIDs`, `excluded_chr_counts`), not just the PRS vector.

## PGS Catalog file format

Scoring files should be downloaded directly from [pgscatalog.org](https://www.pgscatalog.org/) in the standard tab-delimited `.txt.gz` format. Harmonised position columns (`hm_chr`, `hm_pos`) are preferred over raw columns (`chr_name`, `chr_position`) when both are present.

## Genotype file requirements

Genotype files are found through their index files, so file names don't matter (`chr1.vcf.gz`, `cohortA_part1.vcf.gz`, ... all work). Each `*.vcf.gz.tbi` identifies a VCF and each `*.bdose.bdi` identifies a BinaryDosage file. The chromosome(s) in each file are read from the index, and a leading `chr` is ignored when matching against the PGS files (`chr1` in a VCF matches `1` in a scoring file). An index with no data file next to it is skipped with a warning.

If two files of the format in use contain the same chromosome, `compute_prs()` stops with an error naming the files.

**VCF:**
- BGZF-compressed with the extension `.vcf.gz` and a Tabix index (`.vcf.gz.tbi`); other extensions such as `.vcf.bgz` are not recognised
- Chromosomes are taken from the contig names in the Tabix index
- Must contain a `DS` (dosage) field in the `FORMAT` column

**BinaryDosage:**
- Format 5 (`.bdose` + companion `.bdose.bdi`), e.g. as produced by `BinaryDosage::vcftobd()`
- Chromosomes are taken from the `.bdi`

## Allele effect formula

For SNPs where the effect allele matches `REF`:

```
PRS contribution = weight × (2 − dosage)
```

For SNPs where the effect allele matches `ALT`:

```
PRS contribution = weight × dosage
```

Rows where the effect allele matches neither `REF` nor `ALT` are excluded from the score and counted per chromosome in `excluded_chr_counts`.

## License

MIT © John Morrison
