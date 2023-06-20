#!/bin/bash
##### Load all necessary modules ####
DIROUT=${1}
SCRIPTS=${2}

module load statistics/R/4.2.2
module load bioinfo/Bcftools/1.9
module load bioinfo/samtools/1.14
module load devel/python/Python-3.7.9

Rscript ${SCRIPTS}/lib.r



