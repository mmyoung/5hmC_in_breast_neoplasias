#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


## hypomethylation with histone marks

bam_dir=/data/nym/20220921-WSL-5hmc/data/sample_5hmC_bam
peak_dir=/data/nym/20220921-WSL-5hmc/analysis/make_diff_heatmap
bw_dir=/data/nym/20220921-WSL-5hmc/data/sample_peak_bigwig
out_dir=/data/nym/20220921-WSL-5hmc/analysis/DhMRs_5hmC_in_sample


#conda activate samtools
#cd $bam_dir
#for i in `ls *I_unique.bam`
#do
#	j=`echo ${i%%_*}`

#	samtools reheader -c "grep -v "^@PG"" ${i} >${j}_2.bam
#	rm ${i};rm ${i}.bai
#	samtools sort -o ${j}_unique_sort.bam -@ 4 ${j}_2.bam
#	samtools index ${j}_unique_sort.bam
#done

conda activate base
cd $bam_dir
for i in `ls *I_unique_sort.bam`
do
	j=`echo ${i%%_*}`
	bamCoverage -b ${j}_unique_sort.bam -o $bw_dir/${j}_unique_sort.bigwig
done

for i in `ls $bw_dir`
do
        j=`echo ${i%%_*}`
        computeMatrix scale-regions -R $peak_dir/continuously_high_peak.txt \
                                    -S   $bw_dir/${i} \
                                    -b 3000 -a 3000 \
                                    --skipZeros \
                                    -o $out_dir/${j}_high_mat.gz

        computeMatrix scale-regions -R $peak_dir/continuously_low_peak.txt \
                                    -S   $bw_dir/${i} \
                                    -b 3000 -a 3000 \
                                    --skipZeros \
                                    -o $out_dir/${j}_low_mat.gz
done



