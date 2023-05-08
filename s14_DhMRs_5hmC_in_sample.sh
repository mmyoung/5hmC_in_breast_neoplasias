#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


## hypomethylation with histone marks

conda activate base
peak_dir=/data/nym/20220921-WSL-5hmc/analysis/make_diff_heatmap
bw_dir=/data/nym/20220921-WSL-5hmc/data/sample_peak_bigwig
out_dir=/data/nym/20220921-WSL-5hmc/analysis/DhMRs_5hmC_in_sample

cd $peak_dir
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


