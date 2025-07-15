#!/bin/bash

samtools view -d YS:C_G2A $1 -b >CTOB_SpliceC.bam

samtools view CTOB_SpliceC.bam | gawk -F'XB:Z:' '{print $1}' | gawk -F'\t' '{OFS="\t";print $1,$3,$4,$6}' >CTOB_SpliceC_1

samtools view CTOB_SpliceC.bam | gawk -F'XB:Z:' '{print $2}' | gawk -F'\t' '{OFS="\t";print $1,$2}' >CTOB_SpliceC_2

paste CTOB_SpliceC_1 CTOB_SpliceC_2 >SpliceC_CTOB_tmp.txt

rm CTOB_SpliceC_*
python workflow/scripts/bsbolt2bismark/make_bismark.py SpliceC_CTOB_tmp.txt
