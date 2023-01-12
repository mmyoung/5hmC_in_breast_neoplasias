#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4


cd /data/nym/20220921-WSL-5hmc/analysis
#mkdir CGI_region_distri
touch /data/nym/20220921-WSL-5hmc/analysis/CGI_region_distri/CGI_flank_distri.txt
out_dir=/data/nym/20220921-WSL-5hmc/analysis/CGI_region_distri
ucsc_cgi_bed=/data/nym/20220921-WSL-5hmc/data/UCSC_hg19_CGI.bed

for macs2_bed in `ls /data/nym/20220921-WSL-5hmc/data/sample_peak/*.narrowPeak`
do

	out_base=$(basename $macs2_bed)
#	bedtools closest -d -a <(sortBed -i $macs2_bed) -b <(sortBed -i $ucsc_cgi_bed) >$out_dir/$out_base.cgi_distance
	awk '{if($17==0){print "CGI"}if($17>0 && $17<2000){print "shore"}if($17>2000 && $17<4000){print "shelf"}else{print "sea"}}' $out_dir/$out_base.cgi_distance|sort|uniq -c|awk 'BEGIN{OFS="\t"}{print out_base_2,$0}' out_base_2=$out_base  >>/data/nym/20220921-WSL-5hmc/analysis/CGI_region_distri/CGI_flank_distri.txt

done

for i in `ls *.narrowPeak.cgi_distance`;do awk 'BEGIN{OFS="\t"}{if($17==0){print $0,"CGI"}else if($17>0 && $17<2000){print $0,"shore"}else if($17>2000 && $17<4000){print $0,"shelf"}else{print $0,"sea"}}' ${i} >tmp && mv tmp ${i};done

## calculate the average modification level of 5hmC and 5mC in each region
touch /data/nym/20220921-WSL-5hmc/analysis/CGI_region_distri/HM_avg_level_in_CGI_regions.txt
printf "sample\ttotal_level\ttotal_num\tavg_level\n" >>/data/nym/20220921-WSL-5hmc/analysis/CGI_region_distri/HM_avg_level_in_CGI_regions.txt
for i in `ls /data/nym/20220921-WSL-5hmc/analysis/CGI_region_distri/*.narrowPeak.cgi_distance`
do
	out_base=$(basename $i)
	out_base=`echo ${out_base%_peaks*}`
	awk -F"\t" '{a[$18] += $7;b[$18]++}END{for (i in a) {print out_base_2, i, "\t",a[i],"\t",b[i],"\t",a[i]/b[i]}}' out_base_2=$out_base ${i} >>/data/nym/20220921-WSL-5hmc/analysis/CGI_region_distri/HM_avg_level_in_CGI_regions.txt
done

