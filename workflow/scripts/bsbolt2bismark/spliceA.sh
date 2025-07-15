#!/bin/bash

samtools view -d YS:W_C2T $1 -b > OT_SpliceA.bam

samtools view OT_SpliceA.bam | gawk -F'XB:Z:' '{print $1}' | gawk -F'\t' '{OFS="\t";print $1,$3,$4,$6}' >OT_SpliceA_1

samtools view OT_SpliceA.bam | gawk -F'XB:Z:' '{print $2}' | gawk -F'\t' '{OFS="\t";print $1,$2}' >OT_SpliceA_2

paste OT_SpliceA_1 OT_SpliceA_2 >SpliceA_OT_tmp.txt

rm OT_SpliceA_1
rm OT_SpliceA_2

python workflow/scripts/bsbolt2bismark/make_bismark.py SpliceA_OT_tmp.txt
