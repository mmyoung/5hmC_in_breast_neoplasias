#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4

conda activate rtracklayer
cd /data/nym/20220921-WSL-5hmc/src
Rscript s05_ChIPpeakAnno.R


