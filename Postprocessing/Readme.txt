The following scripts are used to generate all analysis following processing of the raw sequencing data.

Note:
* The scripts should be run in the order specified by the filename prefixes. 
* The file 'config.txt' should be updated with the appropriate paths.

The input of the 1st script is the output of the processing pipeline performed on each replicate.
These files are found in the directory 'data/reconstruct_CpG_states' and are named:
* mESC-7cycles_S1_lanes_merged_CpG_state_quants.tsv (Rep 1, 'JSH9')
* CEG1654-MT01-I1528-001_S1_lanes_merged_CpG_state_quants.tsv (Rep 2, 'JSH11')

Also note that the script '02_alt_filtering_reps_separate_reps_for_GEO.Rmd' applies the filtering to each separate replicate. 
The output is not used in downstream steps described in the main manuscript (this used the dataset of the merged replicates).
However, it was used to generate the comparisons between replicates described in the Supplementary Information - and the datasets of individual datasets deposited in the GEO repository.

This downstream analysis can be performed on a laptop. 
Note, however, that when binning the data (03_binning_bed_bedgraph_bigwig_reps_merged.Rmd), this is done iteratively with subsets of chromosomes to avoid exhausting the memory.
