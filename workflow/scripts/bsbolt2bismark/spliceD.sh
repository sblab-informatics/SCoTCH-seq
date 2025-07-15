#!/bin/bash

samtools view -d YS:W_G2A $1 -b >CTOT_SpliceD.bam

samtools view CTOT_SpliceD.bam | gawk -F'XB:Z:' '{print $1}' | gawk -F'\t' '{OFS="\t";print $1,$3,$4,$6}' >CTOT_SpliceD_1

samtools view CTOT_SpliceD.bam | gawk -F'XB:Z:' '{print $2}' | gawk -F'\t' '{OFS="\t";print $1,$2}' >CTOT_SpliceD_2

paste CTOT_SpliceD_1 CTOT_SpliceD_2 >SpliceD_CTOT_tmp.txt


rm CTOT_SpliceD_*

python workflow/scripts/bsbolt2bismark/make_bismark.py SpliceD_CTOT_tmp.txt

