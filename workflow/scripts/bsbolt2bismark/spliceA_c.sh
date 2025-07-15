#!/bin/bash

samtools view -d YS:C_C2T $1 -b > OB_SpliceA.bam

samtools view OB_SpliceA.bam | gawk -F'XB:Z:' '{print $1}' | gawk -F'\t' '{OFS="\t";print $1,$3,$4,$6}' >OB_SpliceA_1

samtools view OB_SpliceA.bam | gawk -F'XB:Z:' '{print $2}' | gawk -F'\t' '{OFS="\t";print $1,$2}' >OB_SpliceA_2

paste OB_SpliceA_1 OB_SpliceA_2 >SpliceA_OB_tmp.txt

rm OB_SpliceA_1
rm OB_SpliceA_2

python workflow/scripts/bsbolt2bismark/make_bismark.py SpliceA_OB_tmp.txt


