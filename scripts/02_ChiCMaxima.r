#!/usr/bin/env Rscript

# Generate CIS IBED from BEDPE for ChiCMaxima using align2ibed.pl
# - Globs *.hicup.bed in the sample dir
# - Calls: perl align2ibed.pl <bed> <baitmap> <rmap> <out.ibed> 1 2 4 5

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--dir",        type="character", help="Sample output dir (e.g., results/SampleA)"),
  make_option("--sample",     type="character", help="Sample name"),
  make_option("--rmap",       type="character", help="Design .rmap file"),
  make_option("--baitmap",    type="character", help="Design .baitmap file"),
  make_option("--perl_bin",   type="character", default="perl", help="Perl binary"),
  make_option("--align2ibed", type="character", help="Path to align2ibed.pl"),
  make_option("--ibed_out",   type="character", help="Output CIS .ibed path (e.g., results/<sample>/<sample>_Cis.ibed)")
)

opt <- parse_args(OptionParser(option_list=option_list))

req <- c("dir","sample","rmap","baitmap","perl_bin","align2ibed","ibed_out")
missing <- req[!nzchar(unlist(opt[req]))]
if (length(missing) > 0) {
  stop("Missing required options: ", paste(missing, collapse=", "))
}

sample_dir <- opt$dir
if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

# Find the BEDPE produced earlier (*.hicup.bed)
bed_candidates <- Sys.glob(file.path(sample_dir, "*.hicup.bed"))
if (length(bed_candidates) == 0) {
  stop("No *.hicup.bed found in: ", sample_dir, ". Ensure upstream step produced a BEDPE.")
}
# Prefer a file containing the sample name
bed_file <- bed_candidates[grep(opt$sample, basename(bed_candidates))]
if (length(bed_file) == 0) bed_file <- bed_candidates[1]
if (length(bed_file) > 1) bed_file <- bed_file[1]

message("=== ChiCMaxima CIS IBED generation ===")
message("  Sample:     ", opt$sample)
message("  BEDPE:      ", bed_file)
message("  baitmap:    ", opt$baitmap)
message("  rmap:       ", opt$rmap)
message("  out IBED:   ", opt$ibed_out)

# Run align2ibed.pl
cmd <- sprintf('"%s" "%s" "%s" "%s" "%s" 1 2 4 5',
               opt$perl_bin, opt$align2ibed, bed_file, opt$baitmap, opt$rmap, opt$ibed_out)
message("[align2ibed] ", cmd)
status <- system(cmd)
if (status != 0) stop("align2ibed.pl failed with status: ", status)

# Sanity check
if (!file.exists(opt$ibed_out) || file.info(opt$ibed_out)$size == 0) {
  stop("IBED not created or empty: ", opt$ibed_out)
}

message("Done. IBED written to: ", opt$ibed_out)
