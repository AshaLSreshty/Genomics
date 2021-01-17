
### Script to post process the alignment file from mapping -  the WGS datasets #####s


#### Input files from mapping directory #####



#!/bin/bash

cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/MAPPING

#sra=(SRR4009215 SRR4009218 SRR4009233 SRR4009237 SRR4009298 SRR4009227 SRR4009228 SRR4009229 SRR4009230 SRR4009301 SRR4009305)

sra=(SRR4009227 SRR4009228 SRR4009229 SRR4009230)


for i in "${sra[@]}";
	do
	echo $i;
	cd $i;
	picard SortSam I=${i}_RG.bam O=${i}_namesorted.bam SO=queryname

	picard MarkDuplicates I=${i}_namesorted.bam O=${i}_namesorted_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=${i}_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT

	picard SortSam I=${i}_namesorted_mrkdup.bam O=${i}_sorted_mrkdup.bam SO=coordinate

	picard BuildBamIndex I=${i}_sorted_mrkdup.bam

	gatk BaseRecalibrator -R ../../HGINDEX/BWA/hg38.fa -I ${i}_sorted_mrkdup.bam -O ${i}_sorted_mrkdup_bqsr.table --known-sites ../../HGINDEX/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites ../../HGINDEX/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites ../../HGINDEX/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching

	gatk ApplyBQSR -R ../../HGINDEX/BWA/hg38.fa -I ${i}_sorted_mrkdup.bam -O ${i}_mrkdup_bqsr.bam --bqsr-recal-file ${i}_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30

	cd ..
	done



	  
