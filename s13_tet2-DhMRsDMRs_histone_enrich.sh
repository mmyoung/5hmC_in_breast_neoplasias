#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


## hypomethylation with histone marks

conda activate base
peak_dir=/data/nym/20220921-WSL-5hmc/analysis/GSE153251_tet2_ChIPseq/DhMRs_tet2_OL
histone_dir=/data/nym/20220921-WSL-5hmc/data/ENCODE_MCF7_ChIPseq
out_dir=/data/nym/20220921-WSL-5hmc/analysis/GSE153251_tet2_ChIPseq/DhMRs_tet2_OL
#mkdir $out_dir

cd $peak_dir
for i in `ls *_tet2_ol.txt`
do
	j=`echo ${i%_tet2_ol.txt}`
	grep "yes" ${i}|awk -F"\t" '$13>1 && $14<0.05' >${j}_fold1.txt
	grep "yes" ${i}|awk -F"\t" '$13<-1 && $14<0.05' >${j}_foldM1.txt
done

computeMatrix reference-point -S $histone_dir/H3K27ac.bigWig \
                               $histone_dir/H3K27me3.bigWig \
                               $histone_dir/H3K36me3.bigWig \
                               $histone_dir/H3K4me1.bigWig \
                               $histone_dir/H3K4me3.bigWig \
                               $histone_dir/H3K9me3.bigWig \
                             -R UDH_ADH_H_fold1.txt ADH_DCIS_H_fold1.txt DCIS_PT_H_fold1.txt \
                            --referencePoint "center" \
                            --beforeRegionStartLength 3000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o H_fold1_matrix.mat.gz
cd $out_dir
plotHeatmap -m  H_fold1_matrix.mat.gz\
	    --colorMap Oranges \
            -out H_fold1_histone_enrich.pdf

computeMatrix reference-point -S $histone_dir/H3K27ac.bigWig \
                               $histone_dir/H3K27me3.bigWig \
                               $histone_dir/H3K36me3.bigWig \
                               $histone_dir/H3K4me1.bigWig \
                               $histone_dir/H3K4me3.bigWig \
                               $histone_dir/H3K9me3.bigWig \
                             -R UDH_ADH_H_foldM1.txt ADH_DCIS_H_foldM1.txt DCIS_PT_H_foldM1.txt \
			     --referencePoint "center" \
			     --beforeRegionStartLength 3000 \
			     --afterRegionStartLength 3000 \
			     --skipZeros -o H_foldM1_matrix.mat.gz
plotHeatmap -m  H_foldM1_matrix.mat.gz\
            --colorMap Oranges \
            -out H_foldM1_histone_enrich.pdf


computeMatrix reference-point -S $histone_dir/H3K27ac.bigWig \
                               $histone_dir/H3K27me3.bigWig \
                               $histone_dir/H3K36me3.bigWig \
                               $histone_dir/H3K4me1.bigWig \
                               $histone_dir/H3K4me3.bigWig \
                               $histone_dir/H3K9me3.bigWig \
                             -R UDH_ADH_M_fold1.txt ADH_DCIS_M_fold1.txt \
			     --referencePoint "center" \
                            --beforeRegionStartLength 3000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o M_fold1_matrix.mat.gz
plotHeatmap -m  M_fold1_matrix.mat.gz\
            --colorMap Oranges \
      	    -out M_fold1_histone_enrich.pdf

computeMatrix reference-point -S $histone_dir/H3K27ac.bigWig \
                               $histone_dir/H3K27me3.bigWig \
                               $histone_dir/H3K36me3.bigWig \
                               $histone_dir/H3K4me1.bigWig \
                               $histone_dir/H3K4me3.bigWig \
                               $histone_dir/H3K9me3.bigWig \
                             -R UDH_ADH_M_foldM1.txt ADH_DCIS_M_foldM1.txt \
			     --referencePoint "center" \
                             --beforeRegionStartLength 3000 \
                             --afterRegionStartLength 3000 \
                             --skipZeros -o M_foldM1_matrix.mat.gz
plotHeatmap -m  M_foldM1_matrix.mat.gz\
            --colorMap Oranges \
      	    -out M_foldM1_histone_enrich.pdf





