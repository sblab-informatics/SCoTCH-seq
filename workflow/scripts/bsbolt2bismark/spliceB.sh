#!/bin/bash

samtools view -d YS:W_C2T $1 -b > OT_SpliceB.bam

samtools view OT_SpliceB.bam | gawk -F'XB:Z:' '{print $1}' | gawk -F'\t' '{OFS="\t";print $1,$3,$4,$6}' >OT_SpliceB_1

samtools view OT_SpliceB.bam | gawk -F'XB:Z:' '{print $2}' | gawk -F'\t' '{OFS="\t";print $1,$2}' >OT_SpliceB_2

paste OT_SpliceB_1 OT_SpliceB_2 > SpliceB_OT_tmp.txt

# Changing 'rm OT_SpliceB*' to rm OT_SpliceB_*, as otherwise this seems to be deleting the bam file...
rm OT_SpliceB_*

python workflow/scripts/bsbolt2bismark/make_bismark.py SpliceB_OT_tmp.txt
