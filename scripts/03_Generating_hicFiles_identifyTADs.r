
#!/bin/env Rscript


#############################
## script generates hic files for juicebox browser + identifs TADss using Arrowhead
######################################

system("module load Java/1.8.0_212")
system("module load CUDA/9.2.88-GCC-7.3.0-2.30")


options(warn=-1)


args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <Input directory><Sample Name> <rmap file> <baitmap file> <bam file>\n"))
  q(status=1) 
}
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

pre.file=paste0(Dir,"/",Sample.Name,"_pre.txt")
hic.file=paste0(Dir,"/",Sample.Name,".hic")
bed.File=Sys.glob(file.path(Dir, "*.hicup.bed"))
mat.file=paste0(Dir,"/",Sample.Name,".mat")

print(Dir)
print(Sample.Name)
print(bam.file)
print(mat.file)
print(pre.file)
print(hic.file)
print(rmap.file)
print(baitmap.file)


command = sprintf(paste("ls", Dir,"|grep .mat|wc -l"))
done = as.numeric(system(command, intern=T))


mat.file=paste0(Dir,"/",Sample.Name,".mat")

cat("Generating pre hic file ...\n")

system(paste('awk \'{ print $1"\t"$3"\t"$5}\'', mat.file,' |sed \'1d\'' , '>', "x"))
system(paste('mv x', mat.file))
system(paste("perl /path/To/Bam2Juicer.pl ",mat.file," /Path/To/MM10_HindIII.rmap ", pre.file))
system(paste('awk \'{if ($3 > $7){ print $1, $6, $7, $8, $9, $11, $2, $3, $4, $5, $10}else {print}}\' ', pre.file, '>','x'))
system(paste('mv x', pre.file))
cat("Sorting Pre hic file ...\n")
system(paste("sort -k3,3d -k7,7d ", pre.file, ">", "x"))
system(paste('mv x', pre.file))
cat("Generating hic file ...\n")
system(paste("java -Xmx10g -jar /Path/To/juicer_tools.1.6.2_linux_jcuda.0.8.jar pre -q 10 ", pre.file, hic.file," mm10"))
cat("Done!...\n")


