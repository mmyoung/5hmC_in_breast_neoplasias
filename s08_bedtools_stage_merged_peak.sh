#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4

raw_peak_dir=/data/nym/20220921-WSL-5hmc/data/sample_peak
out_dir=/data/nym/20220921-WSL-5hmc/data/stage_merged_peak
mkdir /data/nym/20220921-WSL-5hmc/data/stage_merged_peak
sample_ls=`ls $raw_peak_dir/A*H*.narrowPeak`
bedtools multiinter -i $sample_ls -header |awk 'BEGIN{OFS="\t"}gsub(",",",",$5)==2' >$out_dir/ADH_intersect_peak.bed

sample_ls=`ls $raw_peak_dir/D*H*.narrowPeak`
bedtools multiinter -i $sample_ls -header |awk 'BEGIN{OFS="\t"}gsub(",",",",$5)==2' >$out_dir/DCIS_intersect_peak.bed

sample_ls=`ls $raw_peak_dir/PT*.narrowPeak`
bedtools multiinter -i $sample_ls -header |awk 'BEGIN{OFS="\t"}gsub(",",",",$5)==2' >$out_dir/PT_intersect_peak.bed

sample_ls=`ls $raw_peak_dir/[UJ]*H*.narrowPeak`
bedtools multiinter -i $sample_ls -header |awk 'BEGIN{OFS="\t"}gsub(",",",",$5)==2' >$out_dir/UDH_intersect_peak.bed


