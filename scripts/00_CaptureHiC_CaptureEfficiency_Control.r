
#!/bin/env Rscript


#############################
##script calculates the capture efficiency
######################################

options(warn=-1)


args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <Sample input directory: hicup output directory> <Sample Name> <rmap file> <baitmap file> <bam file>\n"))
  q(status=1) 
  
}

for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
print(Dir)
print(Sample.Name)
print(rmap.file)
print(baitmap.file)
print(bam.file)

mat.file=paste0(Dir,"/",Sample.Name,".mat")
Statics.file=paste0(Dir,"/",Sample.Name,"_CaptureEfficincy.stat")


system(paste('awk \'{ sum += $5} END{ print "Total:",sum}\'', mat.file, '>', Statics.file))

system(paste("awk -F '\t' '($2==0 && $4==0) { sum += $5} END{ print \"P0:\",sum}' ", mat.file, '>>',  Statics.file))

system(paste("awk -F '\t' '($2==0 && $4==1) || ($2==1 && $4==0) { sum += $5} END{ print \"P1:\",sum}' ", mat.file, '>>',  Statics.file))

system(paste("awk -F '\t' '($2==1 && $4==1) { sum += $5} END{ print \"P2:\",sum}' ", mat.file, '>>',  Statics.file))

