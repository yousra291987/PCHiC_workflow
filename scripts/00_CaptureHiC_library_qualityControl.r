#!/usr/bin/env Rscript

# Capture-HiC QC + mapping wrapper
# - Runs HiCUP with provided config + digest
# - Converts resulting hicup BAM -> BEDPE (bedtools)
# - Generates MAT table via BamToFragments.pl
# - Produces/links a standard mapped.bam for downstream rules

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--output",          type="character", help="Output directory for this sample (e.g., results/SampleA)"),
  make_option("--sample",          type="character", help="Sample name"),
  make_option("--probeinfo",       type="character", help="Probe info file (ProbeInfo.txt) [not used here but kept for future QC]"),
  make_option("--rmap",            type="character", help="Design rmap file"),
  make_option("--baitmap",         type="character", help="Design baitmap file"),
  make_option("--hicup_config",    type="character", help="HiCUP config file"),
  make_option("--digest",          type="character", help="Digest file path for HiCUP"),
  make_option("--hicup_bin",       type="character", default="hicup", help="hicup binary"),
  make_option("--bedtools_bin",    type="character", default="bedtools", help="bedtools binary"),
  make_option("--samtools_bin",    type="character", default="samtools", help="samtools binary"),
  make_option("--perl_bin",        type="character", default="perl", help="perl binary"),
  make_option("--bam_to_fragments",type="character", help="Path to BamToFragments.pl"),
  make_option("--mapped_bam",      type="character", help="Final mapped BAM path to create (results/<sample>/mapped.bam)"),
  make_option("--bed_out",         type="character", help="Output BEDPE path (results/<sample>/<sample>.hicup.bed)"),
  make_option("--mat_out",         type="character", help="Output MAT file path (results/<sample>/<sample>.mat)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Basic checks
req <- c("output","sample","hicup_config","hicup_bin",
         "bedtools_bin","perl_bin","bam_to_fragments","mapped_bam","bed_out","mat_out",
         "rmap","baitmap")
missing <- req[!nzchar(unlist(opt[req]))]
if (length(missing) > 0) {
  stop("Missing required options: ", paste(missing, collapse=", "))
}

outdir  <- opt$output
sample  <- opt$sample

if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

message(">> Running HiCUP QC/mapping for sample: ", sample)

# 1) Run HiCUP
# Note: We pass both --config and --digest explicitly.
# HiCUP usually writes outputs into a subdir per config; we will glob afterwards.
# Build HiCUP cmd WITHOUT forcing --digest (it’s in the config):
hicup_cmd <- sprintf('%s --config "%s"', opt$hicup_bin, opt$hicup_config)
message("   [HiCUP] ", hicup_cmd)
status <- system(hicup_cmd)
if (status != 0) stop("HiCUP failed with status: ", status)

# 2) Collect ALL HiCUP outputs for this sample.
#    HiCUP processes each FASTQ pair independently and writes ONE output per
#    pair (see HiCUP docs: "a file pair generates one output file"). So a sample
#    with several technical replicates yields several *.hicup.bam files, which we
#    merge below. HiCUP writes SAM instead of BAM when samtools is unavailable,
#    so we accept both and convert any SAM to BAM first.
glob_both <- function(root) {
  c(Sys.glob(file.path(root, "*.hicup.bam")),
    Sys.glob(file.path(root, "*.hicup.sam")))
}
hicup_out <- glob_both(outdir)
if (length(hicup_out) == 0) {
  # also search one level down, in case HiCUP used a subdir
  hicup_out <- unlist(lapply(Sys.glob(file.path(outdir, "*")), glob_both))
}
if (length(hicup_out) == 0) {
  stop("No .hicup.bam/.hicup.sam found under: ", outdir,
       ". Ensure HiCUP outputs into the same --output folder or move/symlink results there.")
}
hicup_out <- sort(unique(hicup_out))
message("   Found ", length(hicup_out), " HiCUP output(s): ",
        paste(basename(hicup_out), collapse=", "))

# Convert any SAM to BAM; keep a de-duplicated list of BAMs to merge.
to_bam <- function(path) {
  if (grepl("\\.bam$", path)) return(path)
  bam <- sub("\\.sam$", ".bam", path)
  if (!file.exists(bam)) {
    cmd <- sprintf('%s view -b "%s" -o "%s"', opt$samtools_bin, path, bam)
    message("   [samtools] ", cmd)
    if (system(cmd) != 0) stop("samtools view (SAM->BAM) failed for: ", path)
  }
  bam
}
bam_files <- unique(vapply(hicup_out, to_bam, character(1)))

# 3) Create/overwrite the standard mapped.bam expected downstream.
#    One replicate -> symlink/copy it; several -> merge them with samtools.
#    Note: HiCUP de-duplicates each pair separately, so a merge of technical
#    replicates does not remove duplicates that span pairs (usually negligible).
if (file.exists(opt$mapped_bam)) file.remove(opt$mapped_bam)
if (length(bam_files) == 1) {
  message("   Single HiCUP BAM -> mapped.bam: ", bam_files[1])
  success <- FALSE
  try({
    file.symlink(from = bam_files[1], to = opt$mapped_bam)
    success <- TRUE
  }, silent = TRUE)
  if (!success) file.copy(from = bam_files[1], to = opt$mapped_bam, overwrite = TRUE)
} else {
  message("   Merging ", length(bam_files), " replicate BAMs -> mapped.bam")
  merge_cmd <- sprintf('%s merge -f "%s" %s',
                       opt$samtools_bin, opt$mapped_bam,
                       paste(sprintf('"%s"', bam_files), collapse=" "))
  message("   [samtools] ", merge_cmd)
  if (system(merge_cmd) != 0) stop("samtools merge failed for sample: ", sample)
}
if (!file.exists(opt$mapped_bam)) stop("Could not create mapped.bam at: ", opt$mapped_bam)

# 4) BAM -> BEDPE
bed_cmd <- sprintf('%s bamtobed -bedpe -i "%s" > "%s"',
                   opt$bedtools_bin, opt$mapped_bam, opt$bed_out)
message("   [bedtools] ", bed_cmd)
status <- system(bed_cmd)
if (status != 0) stop("bedtools bamtobed failed with status: ", status)
if (!file.exists(opt$bed_out) || file.info(opt$bed_out)$size == 0) {
  stop("BED output not created or empty: ", opt$bed_out)
}

# 5) Generate MAT table (for capture efficiency)
#    BamToFragments.pl expects: <bed> <baitmap> <rmap> <out> 1 2 4 5
mat_cmd <- sprintf('%s "%s" "%s" "%s" "%s" "%s" 1 2 4 5',
                   opt$perl_bin, opt$bam_to_fragments, opt$bed_out,
                   opt$baitmap, opt$rmap, opt$mat_out)
message("   [perl] ", mat_cmd)
status <- system(mat_cmd)
if (status != 0) stop("BamToFragments.pl failed with status: ", status)
if (!file.exists(opt$mat_out) || file.info(opt$mat_out)$size == 0) {
  stop("MAT output not created or empty: ", opt$mat_out)
}

message(">> Done.")