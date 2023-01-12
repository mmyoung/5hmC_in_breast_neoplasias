#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


cd /data/nym/20220921-WSL-5hmc/analysis
touch /data/nym/20220921-WSL-5hmc/analysis/CGI_region_distri/diff_peak_CGI_flank_distri.txt
out_dir=/data/nym/20220921-WSL-5hmc/analysis/diff_peak_CGI_region_distri

ucsc_cgi_bed=/data/nym/20220921-WSL-5hmc/data/UCSC_hg19_CGI.bed

mkdir $out_dir
printf "sample\tpeak_counts\tpeak_change\tCGI_region\n" >${out_dir}/diff_peak_CGI_flank_distri.txt
cd /data/nym/20220921-WSL-5hmc/data/diff_peak

for i in `ls *diff_peak_deseq2.txt`
do

awk '$10>1 && $11<0.05' ${i}| awk 'BEGIN{OFS="\t"}{$1=""}1'|awk 'BEGIN{OFS="\t"}{$1=$1}1'|sed s/\"//g |sed "1d" >${out_dir}/${i}_up
awk '$10<-1 && $11<0.05' ${i}|awk 'BEGIN{OFS="\t"}{$1=""}1'|awk 'BEGIN{OFS="\t"}{$1=$1}1'|sed s/\"//g |sed "1d" >${out_dir}/${i}_down
bedtools closest -d -a <(sortBed -i ${out_dir}/${i}_up) -b <(sortBed -i $ucsc_cgi_bed) >$out_dir/${i}_up.cgi_distance
bedtools closest -d -a <(sortBed -i ${out_dir}/${i}_down) -b <(sortBed -i $ucsc_cgi_bed) >$out_dir/${i}_down.cgi_distance
awk '{if($18==0){print "CGI"}if($18>0 && $18<2000){print "shore"}if($18>2000 && $18<4000){print "shelf"}else{print "sea"}}' $out_dir/${i}_up.cgi_distance|sort|uniq -c|awk 'BEGIN{OFS="\t"}{print out_base_2,diff_change,$0}' out_base_2=${i} diff_change="up" >>${out_dir}/diff_peak_CGI_flank_distri.txt
awk '{if($18==0){print "CGI"}if($18>0 && $18<2000){print "shore"}if($18>2000 && $18<4000){print "shelf"}else{print "sea"}}' $out_dir/${i}_down.cgi_distance|sort|uniq -c|awk 'BEGIN{OFS="\t"}{print out_base_2,diff_change,$0}' out_base_2=${i} diff_change="down" >>${out_dir}/diff_peak_CGI_flank_distri.txt

done

for i in `ls ${out_dir}/*.cgi_distance`;do awk 'BEGIN{OFS="\t"}{if($18==0){print $0,"CGI"}else if($18>0 && $18<2000){print $0,"shore"}else if($18>2000 && $18<4000){print $0,"shelf"}else{print $0,"sea"}}' ${i} >tmp && mv tmp ${i};done

## calculate the average modification level of 5hmC and 5mC in each region
touch ${out_dir}/HM_diff_avg_level_in_CGI_regions.txt
printf "sample\ttotal_Fold\ttotal_num\tavg_Fold\n" >>${out_dir}/HM_diff_avg_level_in_CGI_regions.txt
for i in `ls ${out_dir}/*_up.cgi_distance`
do
        out_base=$(basename $i)
        awk -F"\t" '{a[$19] += $9;b[$19]++}END{for (i in a) {print out_base_2, i, "\t",a[i],"\t",b[i],"\t",a[i]/b[i]}}' out_base_2=$out_base ${i} >>${out_dir}/HM_diff_avg_level_in_CGI_regions.txt
done

for i in `ls ${out_dir}/*_down.cgi_distance`
do
        out_base=$(basename $i)
        awk -F"\t" '{a[$19] += $9;b[$19]++}END{for (i in a) {print out_base_2, i, "\t",a[i],"\t",b[i],"\t",a[i]/b[i]}}' out_base_2=$out_base ${i} >>${out_dir}/HM_diff_avg_level_in_CGI_regions.txt
done
