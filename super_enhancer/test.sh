#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4

cd /data/nym/tools/rose
conda activate py27
python ROSE_main.py -g hg19 -i /data/nym/20220921-WSL-5hmc/analysis/enhancer_identify/MCF7_enhancer.gff -r /data/nym/20220921-WSL-5hmc/data/ENCODE_MCF7_ChIPseq/MCF7_H3K27ac_align_rep1_sorted.bam -o /data/nym/20220921-WSL-5hmc/analysis/enhancer_identify

