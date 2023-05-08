#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


## hypomethylation with histone marks

conda activate base
peak_dir=/data/nym/20220921-WSL-5hmc/data/ENCODE_MCF7_ChIPseq
histone_dir=/data/nym/20220921-WSL-5hmc/data/ENCODE_MCF7_ChIPseq
out_dir=/data/nym/20220921-WSL-5hmc/analysis/MCF7_TFs_correlation


cd $peak_dir
computeMatrix reference-point -S $histone_dir/H3K27ac.bigWig \
                               $histone_dir/H3K27me3.bigWig \
                               $histone_dir/H3K36me3.bigWig \
                               $histone_dir/H3K4me1.bigWig \
                               $histone_dir/H3K4me3.bigWig \
                              $histone_dir/H3K9me3.bigWig \
                             -R FOXA1_ChIPseq_MCF7.bed \
                            --referencePoint "center" \
                            --beforeRegionStartLength 3000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $out_dir/FOXA1_histone_cov.mat.gz
cd $out_dir
plotHeatmap -m FOXA1_histone_cov.mat.gz \
	    --colorMap Oranges \
      -out FOXA1_histone_cov.pdf


cd $peak_dir
computeMatrix reference-point -S $histone_dir/H3K27ac.bigWig \
                               $histone_dir/H3K27me3.bigWig \
                               $histone_dir/H3K36me3.bigWig \
                               $histone_dir/H3K4me1.bigWig \
                               $histone_dir/H3K4me3.bigWig \
                              $histone_dir/H3K9me3.bigWig \
                             -R GATA3_ChIPseq_MCF7.bed \
                            --referencePoint "center" \
                            --beforeRegionStartLength 3000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $out_dir/GATA3_histone_cov.mat.gz
cd $out_dir
plotHeatmap -m GATA3_histone_cov.mat.gz \
            --colorMap Oranges \
      -out GATA3_histone_cov.pdf
