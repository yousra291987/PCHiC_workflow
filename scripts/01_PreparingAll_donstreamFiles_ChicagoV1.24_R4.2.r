#!/bin/env Rscript

#First Need to generate the Design Dir (one time before running Chicago)
#Using chicagoTools (makeDesignFiles_py3.py) for python >= 3

#system(paste("python /Path/To/chicagoTools/makeDesignFiles_py3.py --designDir=/Path/To/MM10_HindIII_Design  --outfilePrefix=MM10_HindIII --minFragLen=150 --maxFragLen=40000 --maxLBrownEst=1500000 --binsize=20000 --removeb2b=True --removeAdjacent=True --rmapfile=/Path/To/MM10_HindIII.rmap  --baitmapfile=/Path/To/MM10_HindIII.baitmap"))

#############################
##script identify Sig interactions in Chicago
######################################

options(warn=-1)


args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <Dir> <Design.Dir> <Sample.Name> <rmap.file> <baitmap.file> <bam.file>\n"))
  q(status=1) 
  
}
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(Dir)
print(Sample.Name)
print(rmap.file)
print(baitmap.file)
print(Design.Dir)
#print(bam.file)

bed.File=Sys.glob(file.path(Dir, "*.hicup.bed"))
bam.file=Sys.glob(file.path(Dir, "*.hicup.bam"))


#system(paste("cd", Dir))
cat("Generating chinput file ...\n")

system(paste("/Path/To/chicagoTools/bam2chicago.sh",bam.file,baitmap.file,rmap.file,Sample.Name))
system(paste0("mv ",Sample.Name,"/*.chinput"," ."))
system(paste0("mv ",Sample.Name,"/*.bedpe"," ."))
system(paste0("rm -r ",Sample.Name,"/"))

Chinput.File=Sys.glob(file.path(Dir, "*.chinput"))
print(Chinput.File)
cat("Running Chicago...\n")
options(scipen=999)
library(Chicago)
cd <- setExperiment(designDir = Design.Dir)


cd <- readAndMerge(files= Chinput.File, cd=cd)

cd <- chicagoPipeline(cd,outprefix =paste0(Dir,"/",Sample.Name,"_Chicago_output"))

exportResults(cd, file.path(Dir,paste0(Sample.Name,"_Chicago_output")))

save.image(file=paste0(Dir,"/",Sample.Name,"_Chicago_output.RDa"))


