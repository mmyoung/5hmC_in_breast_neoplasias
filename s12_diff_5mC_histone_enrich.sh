#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


## hypomethylation with histone marks

conda activate base
peak_dir=/data/nym/20220921-WSL-5hmc/data/diff_peak
histone_dir=/data/nym/20220921-WSL-5hmc/data/ENCODE_MCF7_ChIPseq
out_dir=/data/nym/20220921-WSL-5hmc/analysis/diff_peak_with_MCF7_histone
#mkdir $out_dir

cd $peak_dir
computeMatrix reference-point -S $histone_dir/H3K27ac.bigWig \
                               $histone_dir/H3K27me3.bigWig \
                               $histone_dir/H3K36me3.bigWig \
                               $histone_dir/H3K4me1.bigWig \
                               $histone_dir/H3K4me3.bigWig \
                               $histone_dir/H3K9me3.bigWig \
                             -R ADH_DCIS_M_fold1.txt UDH_ADH_M_fold1.txt \
                            --referencePoint "center" \
                            --beforeRegionStartLength 3000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $out_dir/hypo_5mC_matrix.mat.gz
cd $out_dir
plotHeatmap -m hypo_5mC_matrix.mat.gz \
	    --colorMap Oranges \
      -out hypo_5mC_histone_enrich.pdf

cd $peak_dir
computeMatrix reference-point -S $histone_dir/H3K27ac.bigWig \
                               $histone_dir/H3K27me3.bigWig \
                               $histone_dir/H3K36me3.bigWig \
                               $histone_dir/H3K4me1.bigWig \
                               $histone_dir/H3K4me3.bigWig \
                              $histone_dir/H3K9me3.bigWig \
                             -R ADH_DCIS_M_foldM1.txt UDH_ADH_M_foldM1.txt \
                            --referencePoint "center" \
                            --beforeRegionStartLength 3000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $out_dir/hyper_5mC_matrix.mat.gz
cd $out_dir
plotHeatmap -m hyper_5mC_matrix.mat.gz \
	    --colorMap Oranges \
      -out hyper_5mC_histone_enrich.pdf



