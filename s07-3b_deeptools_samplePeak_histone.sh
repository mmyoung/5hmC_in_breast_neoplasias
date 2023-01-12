#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4

conda activate base
peak_dir=/data/nym/20220921-WSL-5hmc/data/sample_peak
MCF7_histone_dir=/data/nym/20220921-WSL-5hmc/data/ENCODE_MCF7_ChIPseq
BE_histone_dir=/data/nym/20220921-WSL-5hmc/data/ENCODE_BE_ChIPseq
out_dir=/data/nym/20220921-WSL-5hmc/analysis/single_sample_histone_enrich

mkdir $out_dir

sample=A10
BaseName="A10"
mkdir $out_dir/$BaseName
computeMatrix scale-regions -S $BE_histone_dir/H3K27ac.bigWig \
			       $BE_histone_dir/H3K27me3.bigWig \
			       $BE_histone_dir/H3K36me3.bigWig \
			       $BE_histone_dir/H3K4me1.bigWig \
			       $BE_histone_dir/H3K4me3.bigWig \
			       $BE_histone_dir/H3K9me3.bigWig \
			    -R $peak_dir/${BaseName}H_peaks.narrowPeak \
			    --beforeRegionStartLength 3000 \
                            --regionBodyLength 5000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $out_dir/$BaseName/matrix.mat.gz
cd $out_dir/$BaseName
plotHeatmap -m matrix.mat.gz \
      -out $BaseName.png

BaseName="D10"
mkdir $out_dir/$BaseName
computeMatrix scale-regions -S $MCF7_histone_dir/H3K27ac.bigWig \
                               $MCF7_histone_dir/H3K27me3.bigWig \
                               $MCF7_histone_dir/H3K36me3.bigWig \
                               $MCF7_histone_dir/H3K4me1.bigWig \
                               $MCF7_histone_dir/H3K4me3.bigWig \
                               $MCF7_histone_dir/H3K9me3.bigWig \
                            -R $peak_dir/${BaseName}H_peaks.narrowPeak \
                            --beforeRegionStartLength 3000 \
                            --regionBodyLength 5000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $out_dir/$BaseName/matrix.mat.gz
cd $out_dir/$BaseName
plotHeatmap -m matrix.mat.gz \
      -out $BaseName.png


BaseName="J"
mkdir $out_dir/$BaseName
computeMatrix scale-regions -S $BE_histone_dir/H3K27ac.bigWig \
                               $BE_histone_dir/H3K27me3.bigWig \
                               $BE_histone_dir/H3K36me3.bigWig \
                               $BE_histone_dir/H3K4me1.bigWig \
                               $BE_histone_dir/H3K4me3.bigWig \
                               $BE_histone_dir/H3K9me3.bigWig \
                            -R $peak_dir/${BaseName}H_peaks.narrowPeak \
                            --beforeRegionStartLength 3000 \
                            --regionBodyLength 5000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $out_dir/$BaseName/matrix.mat.gz
cd $out_dir/$BaseName
plotHeatmap -m matrix.mat.gz \
      -out $BaseName.png


BaseName="PT7"
mkdir $out_dir/$BaseName
computeMatrix scale-regions -S $MCF7_histone_dir/H3K27ac.bigWig \
                               $MCF7_histone_dir/H3K27me3.bigWig \
                               $MCF7_histone_dir/H3K36me3.bigWig \
                               $MCF7_histone_dir/H3K4me1.bigWig \
                               $MCF7_histone_dir/H3K4me3.bigWig \
                               $MCF7_histone_dir/H3K9me3.bigWig \
                            -R $peak_dir/${BaseName}H_peaks.narrowPeak \
                            --beforeRegionStartLength 3000 \
                            --regionBodyLength 5000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $out_dir/$BaseName/matrix.mat.gz
cd $out_dir/$BaseName
plotHeatmap -m matrix.mat.gz \
      -out $BaseName.png

