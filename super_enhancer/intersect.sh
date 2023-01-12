#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


#cd /data/nym/20220921-WSL-5hmc/analysis/enhancer_identify

#bedtools intersect -wa -wb -a ./enhancer_related_diff_5hmC/DCIS_PT_H_foldM1.txt -b <(cut -f 2- MCF7_enhancer_sort_anno.txt|sed "1d") >DCIS_PT_H_foldM1_enhancer.txt


#bedtools intersect -wa -wb -a ./enhancer_related_diff_5hmC/UDH_ADH_H_fold1.txt -b <(cut -f 2- MCF7_enhancer_sort_anno.txt|sed "1d") >UDH_ADH_H_fold1_enhancer.txt


## get SE and non-SE bed filed for plotHeatmap

#awk '!/^#/&&$9==1' MCF7_enhancer_AllEnhancers.table.txt |cut -f 2- >SE_all.bed

# 1165 SEs
#awk '!/^#/&&$9==0' MCF7_enhancer_AllEnhancers.table.txt |cut -f 2- |shuf -n 1165 >non-SE_enhancer_shuf1165.bed

out_dir=/data/nym/20220921-WSL-5hmc/analysis/enhancer_identify
histone_dir=/data/nym/20220921-WSL-5hmc/data/ENCODE_MCF7_ChIPseq
hmC_bw_dir=/data/nym/20220921-WSL-5hmc/data/sample_peak_bigwig
computeMatrix reference-point -S $histone_dir/H3K27ac.bigWig \
                                $hmC_bw_dir/A10H_unique_sort.bigwig \
                                $hmC_bw_dir/D10H_unique_sort.bigwig \
				$hmC_bw_dir/JH_unique_sort.bigwig \
				$hmC_bw_dir/PT7H_unique_sort.bigwig \
                             -R SE_all.bed non-SE_enhancer_shuf1165.bed \
                            --referencePoint "center" \
                            --beforeRegionStartLength 25000 \
                            --afterRegionStartLength 25000 \
                            --skipZeros -o $out_dir/SE_nonSE_computeMat.mat.gz
cd $out_dir
plotHeatmap -m SE_nonSE_computeMat.mat.gz \
            --colorMap Oranges \
      -out SE_nonSE_computeMat.pdf


