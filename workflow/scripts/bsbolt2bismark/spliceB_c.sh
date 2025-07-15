#!/bin/bash

samtools view -d YS:C_C2T $1 -b > OB_SpliceB.bam

samtools view OB_SpliceB.bam | gawk -F'XB:Z:' '{print $1}' | gawk -F'\t' '{OFS="\t";print $1,$3,$4,$6}' >OB_SpliceB_1

samtools view OB_SpliceB.bam | gawk -F'XB:Z:' '{print $2}' | gawk -F'\t' '{OFS="\t";print $1,$2}' >OB_SpliceB_2

paste OB_SpliceB_1 OB_SpliceB_2 > SpliceB_OB_tmp.txt

rm OB_SpliceB_*

python workflow/scripts/bsbolt2bismark/make_bismark.py SpliceB_OB_tmp.txt
