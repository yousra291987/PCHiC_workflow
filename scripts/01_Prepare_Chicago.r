#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tools)
})

# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
option_list <- list(
  make_option("--dir",        type="character", help="Sample working/output dir (e.g., results/SampleA)"),
  make_option("--sample",     type="character", help="Sample name"),
  make_option("--rmap",       type="character", help="Path to .rmap"),
  make_option("--baitmap",    type="character", help="Path to .baitmap"),
  make_option("--bam",        type="character", help="Input BAM (e.g., results/SampleA/mapped.bam)"),
  make_option("--design_dir", type="character", help="Chicago design directory"),
  make_option("--out",        type="character", help="Output TSV path for interactions"),
  # optional externals
  make_option("--bam2chicago", type="character", default="", help="Path to bam2chicago.sh (optional)"),
  make_option("--make_design", type="character", default="", help="Path to makeDesignFiles_py3.py (optional)"),
  make_option("--python",      type="character", default="python", help="Python binary"),
  make_option("--force_make_design", type="character", default="FALSE",
              help="Whether to force rebuilding the design directory (TRUE/FALSE)")
)
opt <- parse_args(OptionParser(option_list=option_list))

req <- c("dir", "sample", "rmap", "baitmap", "bam", "design_dir", "out")
missing <- req[!nzchar(unlist(opt[req]))]
if (length(missing) > 0) stop("Missing required options: ", paste(missing, collapse=", "))

force_make_design <- toupper(opt$force_make_design) %in% c("TRUE","T","1","YES","Y")

sample_dir <- opt$dir
sample     <- opt$sample
design_dir <- opt$design_dir

if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive=TRUE, showWarnings=FALSE)

message("=== Chicago preparation for sample: ", sample, " ===")
message(" - Work dir: ", sample_dir)
message(" - Design dir: ", design_dir)

# ------------------------------------------------------------------------------
# 1) Optionally (re)build Chicago design directory using chicagoTools
# ------------------------------------------------------------------------------
if (force_make_design) {
  if (!nzchar(opt$make_design)) {
    stop("force_make_design is TRUE but --make_design path was not provided.")
  }
  cmd <- sprintf('"%s" "%s" --designDir="%s" --outfilePrefix="%s" --minFragLen=150 --maxFragLen=40000 --maxLBrownEst=1500000 --binsize=20000 --removeb2b=True --removeAdjacent=True --rmapfile="%s" --baitmapfile="%s"',
                 opt$python, opt$make_design, design_dir, basename(design_dir),
                 opt$rmap, opt$baitmap)
  message("[makeDesignFiles] ", cmd)
  status <- system(cmd)
  if (status != 0) stop("makeDesignFiles_py3.py failed with status: ", status)
}

# ------------------------------------------------------------------------------
# 2) Create .chinput using bam2chicago.sh if provided
#    Otherwise, use any existing *.chinput in sample_dir
# ------------------------------------------------------------------------------
chinput_files <- character(0)

if (nzchar(opt$bam2chicago)) {
  message("[bam2chicago] creating chinput from BAM...")
  # bam2chicago.sh <bam> <baitmap> <rmap> <samplename>
  # It writes into a subdir named <samplename>; move files back to sample_dir.
  cmd <- sprintf('cd "%s" && "%s" "%s" "%s" "%s" "%s"',
                 sample_dir, opt$bam2chicago, opt$bam, opt$baitmap, opt$rmap, sample)
  message(cmd)
  status <- system(cmd)
  if (status != 0) stop("bam2chicago.sh failed with status: ", status)

  # Move generated files up to sample_dir (if they landed in sample subdir)
  generated_dir <- file.path(sample_dir, sample)
  if (dir.exists(generated_dir)) {
    chinput_gen <- Sys.glob(file.path(generated_dir, "*.chinput"))
    bedpe_gen   <- Sys.glob(file.path(generated_dir, "*.bedpe"))
    for (f in c(chinput_gen, bedpe_gen)) {
      file.rename(f, file.path(sample_dir, basename(f)))
    }
    unlink(generated_dir, recursive=TRUE, force=TRUE)
  }
  chinput_files <- Sys.glob(file.path(sample_dir, "*.chinput"))
}

# If not using bam2chicago or if nothing was produced, look for existing .chinput
if (length(chinput_files) == 0) {
  chinput_files <- Sys.glob(file.path(sample_dir, "*.chinput"))
}
if (length(chinput_files) == 0) {
  stop("No .chinput files found. Provide --bam2chicago or ensure .chinput exists in: ", sample_dir)
}
message("Found chinput file(s): ", paste(basename(chinput_files), collapse=", "))

# ------------------------------------------------------------------------------
# 3) Run Chicago (Bioconductor)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Chicago)
})

message("[Chicago] setExperiment with design dir: ", design_dir)
cd <- setExperiment(designDir = design_dir)

message("[Chicago] readAndMerge ...")
cd <- readAndMerge(files = chinput_files, cd = cd)

outprefix <- file.path(sample_dir, paste0(sample, "_Chicago_output"))
message("[Chicago] pipeline ...")
options(scipen=999)
cd <- chicagoPipeline(cd, outprefix = outprefix)

message("[Chicago] exportResults ...")
exportResults(cd, outprefix)

# Try to locate a typical 'significant interactions' export to feed Snakemake's single TSV
# Common Chicago output names include: *_significant_interactions.tsv
sig_candidates <- Sys.glob(paste0(outprefix, "*significant*interaction*.tsv"))
if (length(sig_candidates) == 0) {
  # fallback: any .tsv produced by exportResults
  sig_candidates <- Sys.glob(paste0(outprefix, "*.tsv"))
}
if (length(sig_candidates) == 0) {
  stop("Chicago produced no TSV outputs under: ", outprefix)
}

# Pick the first candidate as our canonical output
sel <- sig_candidates[1]
file.copy(sel, opt$out, overwrite=TRUE)
if (!file.exists(opt$out) || file.info(opt$out)$size == 0) {
  stop("Failed to write final output TSV: ", opt$out)
}

# Save R session for debugging/reproducibility
savefile <- file.path(sample_dir, paste0(sample, "_Chicago_output.RData"))
save(cd, file = savefile)

message("Done. Output TSV: ", opt$out)
