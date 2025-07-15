#!/bin/bash

samtools view -d YS:W_G2A $1 -b >CTOT_SpliceC.bam

samtools view CTOT_SpliceC.bam | gawk -F'XB:Z:' '{print $1}' | gawk -F'\t' '{OFS="\t";print $1,$3,$4,$6}' >CTOT_SpliceC_1

samtools view CTOT_SpliceC.bam | gawk -F'XB:Z:' '{print $2}' | gawk -F'\t' '{OFS="\t";print $1,$2}' >CTOT_SpliceC_2

paste CTOT_SpliceC_1 CTOT_SpliceC_2 >SpliceC_CTOT_tmp.txt

rm CTOT_SpliceC_*

python workflow/scripts/bsbolt2bismark/make_bismark.py SpliceC_CTOT_tmp.txt

