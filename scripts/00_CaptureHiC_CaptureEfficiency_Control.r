#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)   # fast & robust
})

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
option_list <- list(
  make_option("--dir",    type="character", help="Sample output dir (e.g. results/SampleA)"),
  make_option("--sample", type="character", help="Sample name"),
  make_option("--mat",    type="character", help="Path to <sample>.mat"),
  make_option("--stats",  type="character", help="Output stats file: <sample>_CaptureEfficiency.stat")
)
opt <- parse_args(OptionParser(option_list=option_list))

req <- c("dir", "sample", "mat", "stats")
missing <- req[!nzchar(unlist(opt[req]))]
if (length(missing) > 0) {
  stop("Missing required options: ", paste(missing, collapse = ", "))
}

dir.create(dirname(opt$stats), showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load MAT
# The original awk used:
#   - column 5 as counts (sum += $5)
#   - conditions on columns $2 and $4 (0/1 flags)
# So we expect tab-delimited with at least 5 columns.
# -----------------------------------------------------------------------------
if (!file.exists(opt$mat)) stop("MAT file not found: ", opt$mat)
DT <- tryCatch(fread(opt$mat, header = FALSE, sep = "\t", showProgress = FALSE),
               error = function(e) stop("Failed to read MAT: ", e$message))

if (ncol(DT) < 5) {
  stop("MAT file has fewer than 5 columns, cannot compute capture efficiency.")
}

# Aliases for readability (keep 1-based R indexing)
# V2 and V4 correspond to $2 and $4 in awk; V5 to $5 (counts)
flag2 <- DT[[2]]
flag4 <- DT[[4]]
counts <- suppressWarnings(as.numeric(DT[[5]]))
counts[is.na(counts)] <- 0

Total <- sum(counts)

# P0: ($2==0 && $4==0)
P0 <- sum(counts[flag2 == 0 & flag4 == 0])

# P1: ($2==0 && $4==1) || ($2==1 && $4==0)
P1 <- sum(counts[(flag2 == 0 & flag4 == 1) | (flag2 == 1 & flag4 == 0)])

# P2: ($2==1 && $4==1)
P2 <- sum(counts[flag2 == 1 & flag4 == 1])

# -----------------------------------------------------------------------------
# Write stats (same lines as your awk output)
# -----------------------------------------------------------------------------
lines <- c(
  paste0("Total: ", Total),
  paste0("P0: ", P0),
  paste0("P1: ", P1),
  paste0("P2: ", P2)
)
writeLines(lines, con = opt$stats)

# Also print a brief summary for the log
cat(sprintf(
  "[%s] Capture efficiency summary\nTotal=%s | P0=%s | P1=%s | P2=%s\nStats: %s\n",
  opt$sample, format(Total, big.mark=","), format(P0, big.mark=","),
  format(P1, big.mark=","), format(P2, big.mark=","), opt$stats
))
