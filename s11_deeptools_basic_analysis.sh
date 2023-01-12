#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


# coverage over genome of each sample

conda activate base
data_dir=/data/nym/20220921-WSL-5hmc/data/sample_peak_bigwig
out_dir=/data/nym/20220921-WSL-5hmc/analysis/deeptools_basic_analysis
#mkdir $out_dir
file_ls=`ls $data_dir`
#cd $data_dir
#multiBigwigSummary bins -b $file_ls -l A10 A2 AW D10 D20 D3 J PT3 PT5 PT6 PT7 U3 U5 -o ${out_dir}/all_sample.npz

plotPCA --corData ${out_dir}/all_sample.npz --plotFile ${out_dir}/all_sample_5hmC_PCA.pdf --labels A10 A2 AW D10 D20 D3 J PT3 PT5 PT6 PT7 U3 U5

