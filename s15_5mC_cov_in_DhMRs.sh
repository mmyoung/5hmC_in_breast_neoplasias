#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


## hypomethylation with histone marks

peak_dir=/data/nym/20220921-WSL-5hmc/data/diff_peak
bam_dir=/data/nym/20220921-WSL-5hmc/data/sample_5hmC_bam
out_dir=/data/nym/20220921-WSL-5hmc/analysis/5hmC_5mC_correlation

conda activate samtools
cd $bam_dir
#cp /media1/huangchangcai/stu_wu/Project_s348g01027/unique_bam/*M_unique.bam ./
for i in `ls *M_unique.bam`
do
	j=`echo ${i%%_*}`
	samtools reheader -c "grep -v "^@PG"" ${i} >${j}_2.bam
	rm ${i};rm ${i}.bai
	samtools sort -o ${j}_unique_sort.bam -@ 4 ${j}_2.bam
	samtools index ${j}_unique_sort.bam
done


conda activate base

multiBamSummary BED-file --BED $peak_dir/UDH_ADH_H_fold1.txt \
			 --bamfiles $bam_dir/U3M_unique_sort.bam \
			            $bam_dir/U5M_unique_sort.bam \
				    $bam_dir/JM_unique_sort.bam \
				    $bam_dir/A10M_unique_sort.bam \
				    $bam_dir/A2M_unique_sort.bam \
				    $bam_dir/AWM_unique_sort.bam \
			 --labels U3M U5M JM A10M A2M AWM \
			 -o $out_dir/UDH_ADH_H_fold1_cov.npz \
			 --outRawCounts $out_dir/UDH_ADH_H_fold1_readCounts.tab\
			 --scalingFactors $out_dir/UDH_ADH_H_fold1_ScalingFactor.tab


multiBamSummary BED-file --BED $peak_dir/UDH_ADH_H_foldM1.txt \
                         --bamfiles $bam_dir/U3M_unique_sort.bam \
                                    $bam_dir/U5M_unique_sort.bam \
                                    $bam_dir/JM_unique_sort.bam \
                                    $bam_dir/A10M_unique_sort.bam \
                                    $bam_dir/A2M_unique_sort.bam \
                                    $bam_dir/AWM_unique_sort.bam \
                         --labels U3M U5M JM A10M A2M AWM \
                         -o $out_dir/UDH_ADH_H_foldM1_cov.npz \
                         --outRawCounts $out_dir/UDH_ADH_H_foldM1_readCounts.tab \
                         --scalingFactors $out_dir/UDH_ADH_H_foldM1_ScalingFactor.tab


