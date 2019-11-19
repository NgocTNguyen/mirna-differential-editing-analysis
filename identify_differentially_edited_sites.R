#! /usr/bin/Rscript

options(stringsAsFactors=FALSE)
options(check.names=FALSE)
cmd_args = commandArgs(TRUE)

site = cmd_args[1]
sgroup = strsplit(site, split = "\\.")[[1]][4]
infile = read.table(site, header = FALSE)
infile$editing_rate = infile$V6/infile$V7
group = strsplit(unique(as.character(infile$V1)), split = "\\.")
group = as.data.frame(group)
group = unique(as.vector(t(group[2,])))

r1 = infile[grepl(infile$V1, pattern = group[1]),"editing_rate"]
r2 = infile[grepl(infile$V1, pattern = group[2]),"editing_rate"]
x1 = infile[grepl(infile$V1, pattern = group[1]),"V6"]
x11 = paste(x1, collapse = ",")
x2 = infile[grepl(infile$V1, pattern = group[2]),"V6"]
x22 = paste(x2, collapse = ",")
n1 = infile[grepl(infile$V1, pattern = group[1]),"V7"]
n11 = paste(n1, collapse = ",")
n2 = infile[grepl(infile$V1, pattern = group[2]),"V7"]
n22 = paste(n2, collapse = ",")

if (length(na.omit(r2)) == 0) { # when second group has no observation
	group[2] = sgroup
	outline = c(site, group[1], group[2], x11, n11, x22, n22, mean(x1/n1), mean(x2/n2), NA)
} else if (length(na.omit(r1)) == 1 | length(na.omit(r2)) == 1) { # when the level of a miRNA is not enough for the test => unreliable to test
	outline = c(site, group[1], group[2], x11, n11, x22, n22, mean(x1/n1), mean(x2/n2), NA)
} else {
	test = t.test(r1, r2, paired = FALSE)
	outline = c(site, group[1], group[2], x11, n11, x22, n22, test$estimate[[1]], test$estimate[[2]], test$p.value)
}

write(outline, file = paste("result",site,sep="."), ncolumns = 10, sep = "\t")

