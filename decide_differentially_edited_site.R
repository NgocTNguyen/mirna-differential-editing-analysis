#! /usr/bin/Rscript

cmd_args = commandArgs(TRUE)
options(stringsAsFactors = FALSE)
options(check.names = FALSE)

base = cmd_args[1]
method = cmd_args[2]
group = cmd_args[3]

infile = read.table(paste("result",base,method,group,"txt", sep ="."), header = FALSE)
infile = na.omit(infile)
infile$fdr = p.adjust(infile$V10, method = "fdr")
infile$diff.edited= FALSE
infile$diff.edited[infile$fdr < 0.05] = TRUE
colnames(infile) = c("editing.site", "group1", "group2", "editing.count1", "total.count1", "editing.count2", "total.count2", "mean.editing.level1", "mean.editing.level2", "p.value", "fdr", "is.diff.edited")
write.table(infile, file = paste("decision",base,method,group,"txt", sep ="."), col.names = TRUE, row.names = FALSE, sep = "\t")
