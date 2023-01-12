#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4

conda activate base
peak_dir=/data/nym/20220921-WSL-5hmc/data/stage_merged_peak
histone_dir=/data/nym/20220921-WSL-5hmc/data/ENCODE_BE_ChIPseq
out_dir=/data/nym/20220921-WSL-5hmc/analysis/sample_peak_with_histone_mark_BE
mkdir $out_dir
for i in `ls $peak_dir/*.bed`
do
	BaseName=`basename ${i} .bed`
	mkdir $out_dir/$BaseName
	computeMatrix scale-regions -S $histone_dir/H3K27ac.bigWig \
				       $histone_dir/H3K27me3.bigWig \
				       $histone_dir/H3K36me3.bigWig \
				       $histone_dir/H3K4me1.bigWig \
				       $histone_dir/H3K4me3.bigWig \
				       $histone_dir/H3K9me3.bigWig \
				    -R $i \
				    --beforeRegionStartLength 3000 \
	                            --regionBodyLength 5000 \
	                            --afterRegionStartLength 3000 \
	                            --skipZeros -o $out_dir/$BaseName/matrix.mat.gz
	cd $out_dir/$BaseName
	plotHeatmap -m matrix.mat.gz \
	      -out $BaseName.png
done


