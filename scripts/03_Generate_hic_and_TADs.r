#!/usr/bin/env Rscript

# Generate Juicer .hic from a CHiC MAT using Bam2Juicer.pl + juicer_tools.
# Replaces ad-hoc awk/sort with R logic; no hard-coded paths or modules.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# -------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------
option_list <- list(
  make_option("--dir",        type="character", help="Sample directory (e.g., results/SampleA)"),
  make_option("--sample",     type="character", help="Sample name"),
  make_option("--mat",        type="character", help="Path to <sample>.mat"),
  make_option("--rmap",       type="character", help="Path to <design>.rmap"),
  make_option("--pre",        type="character", help="Path to write <sample>_pre.txt"),
  make_option("--hic",        type="character", help="Path to write <sample>.hic"),
  make_option("--java",       type="character", default="java", help="Java binary"),
  make_option("--juicer_tools", type="character", help="Path to juicer_tools.jar"),
  make_option("--perl",       type="character", default="perl", help="Perl binary"),
  make_option("--bam2juicer", type="character", help="Path to Bam2Juicer.pl"),
  make_option("--genome_id",  type="character", default="mm10", help="Genome ID (mm10, hg19, hg38, ...)"),
  make_option("--mem_gb",     type="integer",   default=10, help="Juicer Java heap in GB"),
  make_option("--q",          type="integer",   default=10, help="Juicer pre -q threshold")
)
opt <- parse_args(OptionParser(option_list=option_list))

req <- c("dir","sample","mat","rmap","pre","hic","juicer_tools","bam2juicer")
missing <- req[!nzchar(unlist(opt[req]))]
if (length(missing) > 0) stop("Missing required options: ", paste(missing, collapse=", "))

if (!dir.exists(opt$dir)) dir.create(opt$dir, recursive=TRUE, showWarnings=FALSE)

message("=== Generate .hic for sample: ", opt$sample, " ===")
message("  MAT:  ", opt$mat)
message("  RMAP: ", opt$rmap)
message("  PRE:  ", opt$pre)
message("  HIC:  ", opt$hic)

# -------------------------------------------------------------------
# 1) Prepare a triplet MAT (col1, col3, col5) without header (R replaces awk/sed)
#    The original pipeline dropped the first row and kept cols 1/3/5.
# -------------------------------------------------------------------
if (!file.exists(opt$mat)) stop("MAT not found: ", opt$mat)
DT <- fread(opt$mat, header = FALSE, sep = "\t", showProgress = FALSE)

if (ncol(DT) < 5) stop("MAT must have at least 5 columns.")

# Drop first row (original sed '1d')
if (nrow(DT) > 0) DT <- DT[-1]

trip <- DT[, .(V1, V3, V5)]
trip_file <- file.path(opt$dir, paste0(opt$sample, "_mat_triplet.tmp"))
fwrite(trip, trip_file, sep = "\t", col.names = FALSE)

# -------------------------------------------------------------------
# 2) Run Bam2Juicer.pl to make PRE
#    Perl Bam2Juicer.pl <mat_triplet> <rmap> <pre_out>
# -------------------------------------------------------------------
cmd_b2j <- sprintf('"%s" "%s" "%s" "%s" "%s"',
                   opt$perl, opt$bam2juicer, trip_file, opt$rmap, opt$pre)
message("[Bam2Juicer] ", cmd_b2j)
status <- system(cmd_b2j)
if (status != 0) stop("Bam2Juicer.pl failed with status: ", status)
if (!file.exists(opt$pre) || file.info(opt$pre)$size == 0) stop("PRE not created or empty: ", opt$pre)

# -------------------------------------------------------------------
# 3) Normalize PRE rows:
#    Original awk ensured for each line that column3 <= column7 by swapping pairs.
#    Implement the same in R, then sort by col3, col7 (numeric, ascending).
# -------------------------------------------------------------------
PRE <- fread(opt$pre, header = FALSE, sep = " ", showProgress = FALSE)
if (ncol(PRE) < 11) {
  # Some Bam2Juicer variants are tab-delimited; try again with tab
  PRE <- fread(opt$pre, header = FALSE, sep = "\t", showProgress = FALSE)
}
if (ncol(PRE) < 11) stop("PRE format unexpected: fewer than 11 columns.")

# Coerce potential numeric cols for comparison/sort
num <- function(x) suppressWarnings(as.numeric(x))
V3 <- num(PRE[[3]]); V7 <- num(PRE[[7]])

swap <- which(V3 > V7)
if (length(swap)) {
  # New order for swapped rows (matches your awk):
  # new = V6 V7 V8 V9 V11  V1 V2 V3 V4 V5  V10
  PRE[swap, `:=`(
    V1 = PRE$V6[swap],  V2 = PRE$V7[swap],  V3 = PRE$V8[swap],  V4 = PRE$V9[swap],  V5 = PRE$V11[swap],
    V6 = PRE$V1[swap],  V7 = PRE$V2[swap],  V8 = PRE$V3[swap],  V9 = PRE$V4[swap],  V10 = PRE$V5[swap]
  )]
  # V11 stays as original V10
  PRE[swap, V11 := PRE$V10[swap]]
}

# Recompute numeric for sorting
V3 <- num(PRE[[3]]); V7 <- num(PRE[[7]])
ord <- order(V3, V7, na.last = TRUE)
PRE <- PRE[ord]

# Write PRE back (space-separated)
fwrite(PRE, opt$pre, sep = " ", col.names = FALSE, quote = FALSE)

# -------------------------------------------------------------------
# 4) Create .hic with juicer_tools (pre)
# -------------------------------------------------------------------
xmx <- paste0("-Xmx", as.integer(opt$mem_gb), "g")
cmd_hic <- sprintf('"%s" %s -jar "%s" pre -q %d "%s" "%s" %s',
                   opt$java, xmx, opt$juicer_tools, as.integer(opt$q), opt$pre, opt$hic, opt$genome_id)
message("[juicer_tools] ", cmd_hic)
status <- system(cmd_hic)
if (status != 0) stop("juicer_tools pre failed with status: ", status)

if (!file.exists(opt$hic) || file.info(opt$hic)$size == 0) stop(".hic not created or empty: ", opt$hic)

# Cleanup temp
unlink(trip_file)

message("Done. HIC written to: ", opt$hic)
