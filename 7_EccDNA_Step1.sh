
### Script to map reads on to the reference genome #####s

#### Input files from TRIMMING directory #####

#!/bin/bash

#mkdir ../eccDNA

cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/eccDNA

#input = "../List.txt"

#sras = ('SRR4009202', 'SRR4009221', 'SRR4009232', 'SRR4009233', 'SRR4009237', 'SRR4009255', 'SRR4009288', 'SRR4009299', 'SRR4009301', 'SRR4009307', 'SRR4009309', 'SRR4009312', 'SRR4009201', 'SRR4009220', 'SRR4009231', 'SRR4009235', 'SRR4009244', 'SRR4009266', 'SRR4009298', 'SRR4009300', 'SRR4009306', 'SRR4009308', 'SRR4009311', 'SRR4009328');

#sra=(SRR4009202 SRR4009221 SRR4009232 SRR4009237 SRR4009255 SRR4009288 SRR4009299 SRR4009301 SRR4009307 SRR4009309 SRR4009312 SRR4009201 SRR4009220 SRR4009231 SRR4009233 SRR4009235 SRR4009244 SRR4009266 SRR4009298 SRR4009300 SRR4009306 SRR4009308 SRR4009311 SRR4009328)

sra=(SRR4009212 SRR4009213 SRR4009215 SRR4009216 SRR4009217 SRR4009218 SRR4009227 SRR4009228 SRR4009229 SRR4009230 SRR4009303 SRR4009304 SRR4009305 SRR2009288 SRR4009298 SRR4009299)

for i in "${sra[@]}";
	do
	echo $i;
	mkdir $i;
	cd $i;
	samtools sort -@ 16 -n ../../MAPPING/$i/${i}_map.sam -o qname_${i}_map.bam
	cd ..
	done



	  
