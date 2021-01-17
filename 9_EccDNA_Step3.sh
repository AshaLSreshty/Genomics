
### Script to identify eccDNA from the WGS datasets #####s

#### Input files from eccDNA directory #####

##Step 3###

##Preparing the files for Circle-Map##


#!/bin/bash

#mkdir ../eccDNA

cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/eccDNA

sra=(SRR4009220 SRR4009233 SRR4009234 SRR4009235 SRR4009237 SRR4009244 SRR4009288 SRR4009298 SRR4009299 SRR4009300 SRR4009301 SRR4009328 SRR4009213 SRR4009215 SRR4009216 SRR4009217 SRR4009218 SRR4009227 SRR4009228 SRR4009229 SRR4009230 SRR4009303 SRR4009304 SRR4009305)

#SRR4009201 SRR4009292 SRR4009212 SRR4009221 SRR4009231 SRR4009232

for i in "${sra[@]}";
	do
	echo $i;
	mkdir $i;
	cd $i;
	Circle-Map Realign -t 16 -i sort_circular_read_candidates_${i}.bam -qbam qname_${i}_map.bam -sbam sorted_${i}_circle.bam -fasta ../../HGINDEX/BWA/hg38.fa -o circle_${i}.bed
	cd ..
	done



	  
