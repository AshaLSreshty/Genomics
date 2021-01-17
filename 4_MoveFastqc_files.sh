
### Script to quality control of fastq.gz files #####s

#### Input files from FASTQ_Files directory #####

#!/bin/bash


while IFS= read -r line
	do
        mv ../FASTQ_Files/$line/*.html /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/FASTQC/$line/
        mv ../FASTQ_Files/$line/*.zip /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/FASTQC/$line/
        done < ../List-2.txt


