#! /bin/sh

###################
# DATA PROCESSING #
###################

awk '$6!~/snp/ && $7!~/-/' subtotal_combined_count_table | sort -u > internal_modification_count_table # select internal modifications only
cut -f1,11,12,13 internal_modification_count_table | sort -u > mirna_count_per_sample

for base in A U C G; do
	awk -v nucleotide="$base" '$6~/-/ && $7~nucleotide' internal_modification_count_table > $mod.modification_count_table # separate internal modifications into each type (A, U, C, G and adar)
done

awk '$6~/adar/' internal_modification_count_table > adar.modification_count_table


###############################################
# IDENTIFICATION OF EDITING SITES PER SAMPLE #
###############################################

for indiv in `cut -f3 data_info.csv -d"," | grep -v Sample`; do
	for base in A U C G adar; do
		echo $indiv $base
		./sumup_modification_counts.py $indiv $base # sum up all counts for each modification position of each miRNA
		./identify_editing_sites.R $indiv $base # identify editing sites for each sample using accumulative binomial distribution

	done
done
rm *per_pos_sum

for base in A U C G adar; do
	cat *.*.${base}_tested_sites_per_sample > $base.tested_sites_all
	cat *.*.${base}_edited_sites_per_sample > $base.edited_sites_all
	cut -f3,8 $base.tested_sites_all | sort -u > $base.list_tested_sites.txt
	grep kit $base.edited_sites_all | cut -f3,8 | sort -u > kit.$base.list_edited_sites.txt
	grep uc $base.edited_sites_all | cut -f3,8 | sort -u > uc.$base.list_edited_sites.txt
	for site in `comm -1 -2 kit.$base.list_edited_sites.txt uc.$base.list_edited_sites.txt | grep mmu | cut -f1`; do
		grep $site $base.tested_sites_all >> $base.comm_edited_sites_all # get list of edited sites reported in both kit and uc data
	done	
	for group in 24h 2w C KA; do
	        grep $group $base.comm_edited_sites_all > $group.$base.comm_edited_sites_all
	       	for site in `cut -f3 $group.$base.comm_edited_sites_all | sort -u | sed s/\"//g`; do # for each miRNA tested in each comparison
                	for method in kit uc; do
                        	grep $site $group.$base.comm_edited_sites_all | grep $method >> $base.$site.$method.$group.txt # write them in a separate file for DE tests
                	done
        	done
	done
done
rm *sites_per_sample

for base in A U C G adar; do
	for method in kit uc; do
		grep -c TRUE $base.mmu*.$method*txt | awk -v FS=":" '{if ($2>1) print $1}' > $base.$method.sites_for_DE_testing # Get the list of sites with at least TWO samples reported with significantly edited sites
	done
done


#################################################
# IDENTIFICATION OF DIFFERENTIALLY EDITED SITES #
#################################################
for base in A U C G adar; do
	for method in kit uc; do
		for site in `cat $base.$method.sites_for_DE_testing`; do
			echo $site
			./identify_differentially_edited_sites.R $site
		done
	done
done

for base in A U C G adar; do
	for method in kit uc; do
		for group in C KA 24h 2w; do
			echo $base $method $group
			cat result.$base.mmu*.$method.$group.txt > result.$base.$method.$group.txt
			./decide_differentially_edited_sites.R $base $method $group
		done
	done
done

rm *mmu*txt


###########
# SUMMARY #
###########

for method in kit uc; do
	for base in A U C G adar; do
		head -n 1 decision.adar.uc.2w.txt > $method.$base.differentially_edited_sites.txt
		cat decision.$base.$method.*.txt | awk '$12~/TRUE/' >> $method.$base.differentially_edited_sites.txt
	done
done


########
# DONE #
########
