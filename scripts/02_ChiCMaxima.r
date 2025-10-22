#!/bin/env Rscript

#############################
## script generates ibed file for ChiCMaxima Cis profile
######################################

options(warn=-1)


args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <Dir>  <Sample.Name> <rmap.file> <baitmap.file> <bam file> \n"))
  q(status=1) 
  
}
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

bed.File=Sys.glob(file.path(Dir, "*.hicup.bed"))
Ibed.Cis.File=paste0(Dir,"/",Sample.Name,"_Cis.ibed")


print(Dir)
print(Sample.Name)
print(rmap.file)
print(baitmap.file)
print(bed.File)
print(Ibed.Cis.File)
#print(bam.file)



cat("Generating ibed file cis interactions ONLY...\n")

system(paste("perl /Path/To/scripts/align2ibed.pl ",bed.File, baitmap.file, rmap.file ,Ibed.Cis.File ," 1 2 4 5"))



