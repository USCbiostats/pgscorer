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
  geno_dir   = "path/to/genotypes",   # directory of indexed VCFs (.vcf.gz + .tbi) and/or BinaryDosage .bdose files (Format 4 or 5)
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
| `geno_dir` | `"."` | Directory containing the genotype files: BGZF VCFs (`*.vcf.gz` with `*.vcf.gz.tbi` alongside) and/or BinaryDosage files (`*.bdose`; Format 5 needs `*.bdose.bdi` alongside, Format 4 does not). File names are arbitrary |
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

### Checking a genotype directory first

`find_geno_files()` takes the same `geno_dir` and `format` arguments as `compute_prs()` and reports what would be used, without running any scoring:

```r
find_geno_files("path/to/genotypes")                   # autodetect, prints a report
files <- find_geno_files("path/to/genotypes", format = "vcf", verbose = FALSE)   # just the data.frame
```

```
Genotype files in path/to/genotypes
  Format used : BinaryDosage (autodetected)
  Files       : 4, covering 3 chromosome(s)

 chrom contig             file version   info n_samples duplicated
     1   chr1       chr1.bdose       4 parsed         3       TRUE
     1   chr1 chr1_again.bdose       4 parsed         3       TRUE
     2   chr2       chr2.bdose       4 bdinfo         3      FALSE
    19  chr19      chr19.bdose       4 parsed         3      FALSE

Potential problems:
  - Chromosome 1 is in more than one file (chr1.bdose, chr1_again.bdose); compute_prs() would stop.
  - 3 Format 4 file(s) have no .bdinfo and are parsed in full each time ...
```

It returns (invisibly) a data.frame with one row per file and chromosome: `chrom` (the form matched against PGS files), `contig` (as stored in the file), `file`, `format`, `version` (BinaryDosage format 4 or 5), `info` (`tbi`, `bdi`, `bdinfo`, or `parsed` for a full parse of the data file), `n_samples`, `duplicated`, and `path`. The problems it lists are: a chromosome present in more than one file (`compute_prs()` would stop), Format 4 files without a `.bdinfo` (slow), and files whose sample IDs differ. A duplicated chromosome is reported here instead of raising an error.

## PGS Catalog file format

Scoring files should be downloaded directly from [pgscatalog.org](https://www.pgscatalog.org/) in the standard tab-delimited `.txt.gz` format. Harmonised position columns (`hm_chr`, `hm_pos`) are preferred over raw columns (`chr_name`, `chr_position`) when both are present.

## Genotype file requirements

Genotype files are found by extension, so file names otherwise don't matter (`chr1.vcf.gz`, `cohortA_part1.vcf.gz`, ... all work). Every `*.vcf.gz` that has a `.tbi` is a VCF, and every `*.bdose` is a BinaryDosage file. The chromosome(s) in each file are read from the file's own metadata, and a leading `chr` is ignored when matching against the PGS files (`chr1` in a VCF matches `1` in a scoring file). An index file (`.tbi` / `.bdose.bdi`) with no data file next to it is skipped with a warning.

If two files of the format in use contain the same chromosome, `compute_prs()` stops with an error naming the files.

**VCF:**
- BGZF-compressed with the extension `.vcf.gz` and a Tabix index (`.vcf.gz.tbi`); other extensions such as `.vcf.bgz` are not recognised
- Chromosomes are taken from the contig names in the Tabix index
- Must contain a `DS` (dosage) field in the `FORMAT` column

**BinaryDosage:**
- Extension `.bdose`. Both **Format 5** (`.bdose` + companion `.bdose.bdi`, e.g. as produced by `BinaryDosage::vcftobd()`) and **Format 4** (a single self-contained `.bdose` file, no `.bdi`) are supported, and the two can be mixed in one directory. Formats 1-3 (separate family/map files) are not supported
- Chromosomes are taken from the `.bdi` (Format 5) or the file header (Format 4)

**Speeding up Format 4 files with a `.bdinfo` file.** Format 4 has no information file, so `getbdinfo()` has to parse the entire `.bdose` file to find its chromosomes, samples and SNP offsets, which is slow for large files. `pgscorer` reads that information from `<name>.bdinfo` instead when it exists next to `<name>.bdose`. Create it once per file:

```r
library(BinaryDosage)
saveRDS(getbdinfo("chr1.bdose"), "chr1.bdinfo")   # same name, .bdose replaced by .bdinfo
```

Notes:
- Format 5 files ignore any `.bdinfo` (they use their `.bdi`).
- The file path stored inside a `.bdinfo` doesn't matter; `pgscorer` always substitutes the real location, so files can be moved or copied after the `.bdinfo` is made.
- A `.bdinfo` that can't be read, or whose SNP offsets don't fit the data file, is ignored with a warning and the `.bdose` is parsed as usual.
- Regenerate the `.bdinfo` whenever its `.bdose` is rewritten. A `.bdinfo` from an older version of the data can't be detected beyond those consistency checks.

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
