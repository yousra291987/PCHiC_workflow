#!/usr/bin/env Rscript

# Calculate insulation scores with FAN-C:
# fanc insulation <hic>@<resolution>@<balancing> <out> -g -w <windows...>

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--dir",        type="character", help="Sample output dir (e.g., results/SampleA)"),
  make_option("--sample",     type="character", help="Sample name"),
  make_option("--hic",        type="character", help="Path to <sample>.hic"),
  make_option("--out",        type="character", help="Output insulation file path"),
  make_option("--resolution", type="character", default="10kb", help="Resolution (e.g., 10kb)"),
  make_option("--balancing",  type="character", default="KR",   help="Balancing scheme (e.g., KR)"),
  make_option("--windows",    type="character", default="0.5mb 1mb 1.5mb 2mb 2.5mb",
              help="Space- or comma-separated window sizes for FAN-C -w (e.g., '0.5mb 1mb 1.5mb')")
)

opt <- parse_args(OptionParser(option_list=option_list))

req <- c("dir","sample","hic","out")
missing <- req[!nzchar(unlist(opt[req]))]
if (length(missing) > 0) stop("Missing required options: ", paste(missing, collapse=", "))

if (!dir.exists(opt$dir)) dir.create(opt$dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(opt$hic)) stop("HIC file not found: ", opt$hic)

# Normalize windows list (allow space or comma separated)
wins <- trimws(opt$windows)
wins <- gsub(",", " ", wins)
wins <- unlist(strsplit(wins, "\\s+"))
wins <- wins[nzchar(wins)]
windows_str <- paste(wins, collapse = " ")

message("=== FAN-C Insulation ===")
message("  Sample:    ", opt$sample)
message("  HIC:       ", opt$hic)
message("  Out:       ", opt$out)
message("  Res:       ", opt$resolution)
message("  Balance:   ", opt$balancing)
message("  Windows:   ", windows_str)

# Build FAN-C command
cmd <- sprintf('fanc insulation "%s@%s@%s" "%s" -g -w %s',
               opt$hic, opt$resolution, opt$balancing, opt$out, windows_str)
message("[FAN-C] ", cmd)

status <- system(cmd)
if (status != 0) stop("FAN-C returned a non-zero exit status: ", status)

if (!file.exists(opt$out) || file.info(opt$out)$size == 0) {
  stop("Insulation output not created or empty: ", opt$out)
}

message("Done. Insulation written to: ", opt$out)
