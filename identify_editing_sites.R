#! /usr/bin/Rscript

require(dplyr)
options(stringsAsFactors=FALSE)
options(check.names=FALSE)

cmd_args <- commandArgs(TRUE)
sample = cmd_args[1]
base = cmd_args[2]

infile = paste(sample, ".", base,"_per_pos_sum", sep ="")
infile = read.table(infile, header = TRUE)
infile$pos = as.numeric(infile$pos)
infile$modification_count = as.numeric(infile$modification_count)
infile$total_mir_count = as.numeric(infile$total_mir_count)

infile = infile[(infile$pos >= 0 & infile$pos <= 28),]
infile = arrange(infile, miRNA, pos)

p = c()
for (i in 1:nrow(infile)){
  t = 1 - pbinom(infile[i,6], infile[i,7], p = 0.01)
  p = c(p, t)
}

infile = data.frame(infile, p.value = p, fdr = p.adjust(p, method = "fdr")) #adjust for multiple testing
infile$edited.site = FALSE
infile$edited.site[infile$fdr < 0.05] = TRUE #FDR threshold of 0.05

outfile = paste(sample,".", base, "_tested_sites_per_sample", sep ="")
write.table(infile, file = outfile, sep = "\t", col.names = TRUE, row.names = FALSE)


infile = infile[infile$edited.site,]
outfile = paste(sample,".", base, "_edited_sites_per_sample", sep ="")
write.table(infile, file = outfile, sep = "\t", col.names = TRUE, row.names = FALSE)
