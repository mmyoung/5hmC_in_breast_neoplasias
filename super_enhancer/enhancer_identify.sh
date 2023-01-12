#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


cd /data/nym/20220921-WSL-5hmc/analysis/enhancer_identify
data_dir=/data/nym/20220921-WSL-5hmc/data/ENCODE_MCF7_ChIPseq
../../data/ENCODE_MCF7_ChIPseq/H3K4me3.bed.gz
bedtools subtract -a <(zcat $data_dir/H3K4me1.bed.gz) -b <(zcat $data_dir/H3K4me3.bed.gz) >MCF7_enhancer.bed
bedtools subtract -a <(zcat $data_dir/H3K4me3.bed.gz) -b <(zcat $data_dir/H3K4me1.bed.gz) >MCF7_promoter.bed


