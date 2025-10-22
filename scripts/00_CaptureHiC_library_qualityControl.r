
#!/bin/env Rscript


#############################
##script maps and filter low quality di-tags using HiCup
######################################


options(warn=-1)


args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <hicup config file Path> <output directory Path> <Sample Name> <Digest file> \n"))
  q(status=1) 
}

for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
print(hicup.config.file)
print(output.Dir)
print(sample)
print(digested)

output.Sample=paste(output.Dir,"/",sample)
cat("Quality controle Analysis: running hicup..\n")
system(paste("/Path/To/hicup_v0.6.1/hicup --config", hicup.config.file,"--digest",digested ))

#Converting bam file to bed file
cat("Converting bam to bed ...\n")

bam.File=Sys.glob(file.path(output.Sample, "*.hicup.bam"))
ouptut.bed.File=strsplit(bam.File,"hicup")
system(paste0("bamToBed -bedpe -i ",bam.File,">",ouptut.bed.File[[1]][1],"hicup.bed"))


##Generating mat table to calculate the capture efficiency
cat("Generating mat table ...\n")
bed.File=Sys.glob(file.path(output.Sample, "*.hicup.bed"))
mat.file=paste0(output.Sample,"/",Sample.Name,".mat")
system(paste0("perl /Path/To/scripts/BamToFragments.pl ",bed.File," ",baitmap.file," ",rmap.file," ",mat.file ," 1 2 4 5"))

date()
sessionInfo()


