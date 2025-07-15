#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate snakemake-slurm

snakemake --unlock 

snakemake --profile slurm/
