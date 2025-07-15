#!/bin/bash

samtools view -d YS:C_G2A $1 -b >CTOB_SpliceD.bam

samtools view CTOB_SpliceD.bam | gawk -F'XB:Z:' '{print $1}' | gawk -F'\t' '{OFS="\t";print $1,$3,$4,$6}' >CTOB_SpliceD_1

samtools view CTOB_SpliceD.bam | gawk -F'XB:Z:' '{print $2}' | gawk -F'\t' '{OFS="\t";print $1,$2}' >CTOB_SpliceD_2

paste CTOB_SpliceD_1 CTOB_SpliceD_2 >SpliceD_CTOB_tmp.txt

rm CTOB_SpliceD_*

python workflow/scripts/bsbolt2bismark/make_bismark.py SpliceD_CTOB_tmp.txt
