
#!/bin/env Rscript


#############################
##This script is for calculating Insulation score
######################################



options(warn=-1)

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <Dir> <Sample.Name> <bam file>\n"))
  q(status=1) 
  
}
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
print(Dir)
print(Sample.Name)
print(bam.file)

#system(paste0("module load fanc/0.9.26"))

cat("Calculating Insulation Score using FAN-C ...\n")
  hic.File=Sys.glob(file.path(Dir, "*.hic"))
  Insul.File=paste0(Dir,"/",Sample.Name,".insulation")
system(paste0("module load fanc/0.9.26; fanc insulation ",hic.File,"@10kb@KR ",Insul.File," -g -w 0.5mb 1mb 1.5mb 2mb 2.5m"))

cat("Done!...\n")



