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
  make_option("--ibed_out",   type="character", help="Output CIS .ibed path (e.g., results/<sample>/<sample>_Cis.ibed)"),
  # ChiCMaxima loop caller (optional). ChiCMaxima is a separate tool
  # (github.com/yousra291987/ChiCMaxima); point --chicmaxima_caller at its
  # ChiCMaxima_Caller.r to call loops here. Empty -> produce the .ibed only.
  make_option("--chicmaxima_caller", type="character", default="",
              help="Path to ChiCMaxima_Caller.r; empty = produce the .ibed input only"),
  make_option("--loops_prefix", type="character", default="",
              help="Output prefix for the loop caller; writes <prefix>_interactions.ibed [default: <dir>/<sample>_ChiCMaxima]"),
  make_option("--window_size", type="integer", default=50L,
              help="ChiCMaxima -w: sliding-window size in restriction fragments [default %default]"),
  make_option("--loess_span", type="double", default=0.05,
              help="ChiCMaxima -s: LOESS span for the virtual-4C smoothing [default %default]"),
  make_option("--cis_window", type="integer", default=1500000L,
              help="ChiCMaxima -c: cis window in bp from the bait start [default %default]")
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

# Run align2ibed.pl:
#   align2ibed.pl <bedpe> <baitmap> <rmap> <output> <r1 chr> <r1 coord> <r2 chr> <r2 coord>
# NOTE: SIX "%s" are required (perl, script, bedpe, baitmap, rmap, output); an
# earlier version had only five, which silently dropped the output path and made
# align2ibed.pl receive "1" as its output argument.
cmd <- sprintf('"%s" "%s" "%s" "%s" "%s" "%s" 1 2 4 5',
               opt$perl_bin, opt$align2ibed, bed_file, opt$baitmap, opt$rmap, opt$ibed_out)
message("[align2ibed] ", cmd)
status <- system(cmd)
if (status != 0) stop("align2ibed.pl failed with status: ", status)

# Sanity check
if (!file.exists(opt$ibed_out) || file.info(opt$ibed_out)$size == 0) {
  stop("IBED not created or empty: ", opt$ibed_out)
}

message("Done. IBED written to: ", opt$ibed_out)

# -----------------------------------------------------------------------------
# Optional: run the ChiCMaxima loop caller on the .ibed.
# ChiCMaxima is a separate tool (github.com/yousra291987/ChiCMaxima). When
# tools.chicmaxima_caller is set we call loops here; when it is empty we stop
# after the .ibed input, with a clear note, so the step still succeeds.
# -----------------------------------------------------------------------------
caller <- opt$chicmaxima_caller
if (!nzchar(caller)) {
  message("[ChiCMaxima] loop caller not configured (tools.chicmaxima_caller is empty); ",
          "produced the .ibed input only.")
  message("  To call loops in the pipeline, set tools.chicmaxima_caller to the ",
          "path of your ChiCMaxima_Caller.r.")
  quit(save = "no", status = 0)
}
if (!file.exists(caller)) {
  stop("tools.chicmaxima_caller is set but the file was not found: ", caller)
}

prefix <- opt$loops_prefix
if (!nzchar(prefix)) prefix <- file.path(sample_dir, paste0(opt$sample, "_ChiCMaxima"))

# Use the same Rscript that is running this wrapper, so the caller runs in the
# same R environment (no reliance on a second Rscript being on PATH).
rscript_bin <- file.path(R.home("bin"), "Rscript")

loop_cmd <- sprintf('"%s" "%s" -i "%s" -o "%s" -w %d -s %s -c %d',
                    rscript_bin, caller, opt$ibed_out, prefix,
                    as.integer(opt$window_size),
                    format(opt$loess_span, scientific = FALSE),
                    as.integer(opt$cis_window))
message("=== ChiCMaxima loop calling ===")
message("[ChiCMaxima_Caller] ", loop_cmd)
status <- system(loop_cmd)
if (status != 0) stop("ChiCMaxima_Caller.r failed with status: ", status)

loops_out <- paste0(prefix, "_interactions.ibed")
if (!file.exists(loops_out) || file.info(loops_out)$size == 0) {
  stop("ChiCMaxima ran but produced no interactions file: ", loops_out)
}
n_lines <- tryCatch(length(readLines(loops_out)), error = function(e) NA_integer_)
message("Done. ChiCMaxima loops: ", loops_out,
        if (!is.na(n_lines)) paste0(" (", n_lines, " lines incl. header)") else "")
