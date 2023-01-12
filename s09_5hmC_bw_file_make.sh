#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l p=4

conda activate samtools
out_dir=/data/nym/20220921-WSL-5hmc/data/sample_peak_bigwig
cd /data/nym/20220921-WSL-5hmc/data/sample_5hmC_bam
#for i in `ls /data/nym/20220921-WSL-5hmc/data/sample_5hmC_bam/*H*_unique_sort.bam|grep -v "A10H"|grep -v "^PT"`
#do
#	base_name=`basename ${i} .bam`
#	samtools reheader -c "grep -v "^@PG"" ${i} >${base_name}_2.bam
#	rm ${i};rm ${i}.bai
#	mv ${base_name}_2.bam ${i}
	#samtools sort -o ${base_name}_sort.bam ${i}
#	samtools index ${i}
	#bamCoverage -b ${base_name}_sort.bam -o $out_dir/$base_name.bigwig
#done


conda activate base
for i in `ls /data/nym/20220921-WSL-5hmc/data/sample_5hmC_bam/*H*_unique_sort.bam|grep -v "A10H"|grep -v "^PT"`
do
	base_name=`basename ${i} .bam`
	bamCoverage -b ${i} -o $out_dir/$base_name.bigwig
done



