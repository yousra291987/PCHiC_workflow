#!/usr/bin/env Rscript
# Minimal stand-in for ChiCMaxima_Caller.r, used by the downstream smoke test so
# the wrapper's loop-calling path can be exercised without installing the real
# ChiCMaxima tool. It mimics the real caller's interface: reads -i/--input and
# -o/--output, then writes <output>_interactions.ibed (the input columns plus an
# Enrichment column) and <output>_poorlycovered.txt.
args <- commandArgs(trailingOnly = TRUE)
get <- function(flags) {
  for (f in flags) { i <- which(args == f); if (length(i)) return(args[i[1] + 1L]) }
  NA_character_
}
input  <- get(c("-i", "--input"))
output <- get(c("-o", "--output"))
if (is.na(input) || is.na(output)) stop("stub ChiCMaxima_Caller: need -i and -o")

lines <- readLines(input)
if (length(lines) < 1) stop("stub ChiCMaxima_Caller: empty input ", input)
hdr  <- lines[1]
body <- if (length(lines) > 1) lines[-1] else character(0)
keep <- head(body, 3)  # pretend these are the called interactions

writeLines(c(paste0(hdr, "\tEnrichment"),
             if (length(keep)) paste0(keep, "\t1.5") else character(0)),
           paste0(output, "_interactions.ibed"))
writeLines("# baits with too little coverage (stub)",
           paste0(output, "_poorlycovered.txt"))
cat("stub ChiCMaxima_Caller: wrote", length(keep), "interactions\n")
