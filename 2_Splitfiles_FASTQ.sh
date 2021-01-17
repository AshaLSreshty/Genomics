### Script to download SRA files from DRA database #####

# GBM - WGS SRA files

#!/bin/bash

cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/FASTQ_Files

while IFS= read -r line;
	do
	echo $line
        mkdir $line
        nohup fastq-dump --gzip --outdir $line --split-files ../Datasets/${line}.1 &
        done < "../List.txt"
	
	  
