# batch_worker.R — dispatched by compute_prs() on non-Windows systems.
# Not intended to be run directly.
# Arguments: <positions_rds> <vcf_path> <chrom> <out_rds>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L)
  stop("Usage: Rscript batch_worker.R <positions_rds> <vcf_path> <chrom> <out_rds>")

positions_rds <- args[1L]
VCF_PATH      <- args[2L]
CHROM         <- args[3L]
out_rds       <- args[4L]

library(pgscorer)

batch_pos <- readRDS(positions_rds)
vcf_hits  <- tabixr::query_vcf_positions(VCF_PATH, CHROM, batch_pos)
matched   <- vcf_hits[vcf_hits$POS %in% batch_pos, ]
rm(vcf_hits)

if (nrow(matched) == 0L) {
  saveRDS(NULL, out_rds)
} else {
  # tabixr's data.frame can rename sample columns ("1001" -> "X1001"); the header has the true IDs.
  saveRDS(extract_dosage(matched, tabixr::vcf_samples(VCF_PATH)), out_rds)
}
