#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


## hypomethylation with histone marks

conda activate base
peak_dir=/data/nym/20220921-WSL-5hmc/data/diff_peak
histone_dir=/data/nym/20220921-WSL-5hmc/data/ENCODE_MCF7_ChIPseq
out_dir=/data/nym/20220921-WSL-5hmc/analysis/MCF7_TFs_correlation

cd $peak_dir
computeMatrix reference-point -S $histone_dir/FOXA1_ChIPseq_MCF7.bigWig \
                               $histone_dir/GATA3_ChIPseq_MCF7.bigWig \
                               $histone_dir/FOS_ChIPseq_MCF7.bigWig \
                               $histone_dir/FOSL2_ChIPseq_MCF7.bigWig \
                               $histone_dir/ESR1_ChIPseq_MCF7_hg19.bigWig.bw \
                             -R UDH_ADH_H_fold1.txt ADH_DCIS_H_fold1.txt DCIS_PT_H_fold1.txt \
                            --referencePoint "center" \
                            --beforeRegionStartLength 3000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $out_dir/hypo_TFs_enrich.mat.gz
cd $out_dir
plotProfile -m hypo_TFs_enrich.mat.gz \
      -out hypo_methy_TFs_enrich.pdf

cd $peak_dir
computeMatrix reference-point -S $histone_dir/FOXA1_ChIPseq_MCF7.bigWig \
                               $histone_dir/GATA3_ChIPseq_MCF7.bigWig \
                               $histone_dir/FOS_ChIPseq_MCF7.bigWig \
                               $histone_dir/FOSL2_ChIPseq_MCF7.bigWig \
                               $histone_dir/ESR1_ChIPseq_MCF7_hg19.bigWig.bw \
                             -R UDH_ADH_H_foldM1.txt ADH_DCIS_H_foldM1.txt DCIS_PT_H_foldM1.txt \
                            --referencePoint "center" \
                            --beforeRegionStartLength 3000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $out_dir/hyper_TFs_enrich.mat.gz
cd $out_dir
plotProfile -m hyper_TFs_enrich.mat.gz \
      -out hyper_methy_TFs_enrich.pdf



