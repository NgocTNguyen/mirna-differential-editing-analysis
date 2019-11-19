This folder contains the following files:

1. data_info.csv file is a comma separated file containing information about sequencing files, groups, methods, biological replicates etc

2. subtotal_combined_count_table represents the format required for the input data for the pipeline (i.e. miRNA_editing_pipeline.sh). The file is a miRNA modification count table combined with miRNA expression count table of all samples, providing:
- all the modification types detected in the samples: 3p/5p/internal modifications, single nucleotide modifications (A/U/C/G/adar) and etc
- positions of each modification
- numbers of miRNA copies (counts) detected with each type of modifications
- and the total expression counts of miRNAs carrying the modifications

This file can be generated from the Chimira webserver (https://www.ebi.ac.uk/research/enright/software/chimira) after submitting all the relevant sequencing data files.

3. miRNA_editing_pipeline.sh is a master bash script to process the raw count data table (i.e. subtotal_combined_count_table) and call all relevant Python (sumup_modification_counts.py) and R (identify_editing_sites.R, identify_differentially_edited_sites.R, decide_differentially_edited_sites.R) scripts to process the data and run all relevant statistical tests.

In order to run the pipeline, open a Unix/Linux terminal, go to this folder and type ./miRNA_editing_pipeline.sh
