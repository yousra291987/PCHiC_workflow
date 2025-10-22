#!/bin/env Rscript


#############################
##This script will digest the reference Genome in silico according to the RE used in the HiC/Capture-HiC experiment
######################################



options(warn=-1)


args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <digester.config.file> <Genome> <RE> <output.Dir> \n"))
  q(status=1) 
}

for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
print(digester.config.file)
print(Genome)
print(RE)
print(output.Dir)

cat("Genome digestion..\n")
system(paste("/Path/To/hicup_v0.6.1/hicup_digester --config",digester.config.file,"--genome", Genome,"--outdir",output.Dir))

Digested <- list.files(pattern = "^Digest_*.")
system(paste0("mv ",Digested," Digest_",Genome,"_",RE,".txt"))
date()
sessionInfo()


