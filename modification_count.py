#! /usr/bin/python

from __future__ import division
import sys

try:
       indiv = sys.argv[1]
       base = sys.argv[2]
except:
       sys.exit("Usage: please provide a sample id")

infile1 = "mirna_count_per_sample"
infile2 = "%s.modification_count_table" %(base)
path="%s.%s_per_pos_sum" %(indiv, base)

outfile = open(path, "w")

pos_dict = {}
count_dict = {}
group = indiv.strip().split("_")
group = "%s_%s" %(group[0], group[1])
print group

for l in open(infile1):
    if l.find(indiv) != -1:
        l = l.strip().split("\t")
        total_count = int(l[1])
        mirna = l[0]
        count_dict[mirna] = total_count

for l in open(infile2):
    if l.find(indiv) != -1:
        array = l.strip().split('\t')
        mirna = array[0]
        if count_dict[mirna] < 10: continue
        pos = int(array[7])
        if pos == 0: pos = 1
        mod_id = "%s@%d" %(mirna, pos)
        if not pos_dict.has_key(mod_id): pos_dict[mod_id] = int(array[9])
        else:
             count = int(array[9]) + pos_dict[mod_id]
             pos_dict[mod_id] = count

outfile.write("group\tindiv\tmodification_id\tmiRNA\tpos\tmodification_count\ttotal_mir_count\tnucleotide\n")

for item in pos_dict:
    array = item.split("@")
    mirna = array[0]
    pos = int(array[1])
    ol = "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\n" %(group, indiv, item, mirna, pos, pos_dict[item], count_dict[mirna], base)
    outfile.write(ol)

outfile.close()
