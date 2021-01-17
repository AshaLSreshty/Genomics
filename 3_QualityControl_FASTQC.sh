
### Script to quality control of fastq.gz files #####s

#### Input files from FASTQ_Files directory #####

#!/bin/bash

#mkdir ../FASTQC

cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/FASTQC

#input = "../List.txt"

while IFS= read -r line
	do
	echo $line;
        mkdir $line;
   
        nohup fastqc ../FASTQ_Files/$line/${line}.1_1.fastq.gz ../FASTQ_Files/$line/${line}.1_2.fastq.gz -o $line &
#        mv ../FASTQ_Files/$line/*.html /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/FASTQC/$line/
#        mv ../FASTQ_Files/$line/*.zip /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/FASTQC/$line/
       
        done < ../List-2.txt

cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/SCRIPTs





#for i in $sranum;
#	do
#	echo $sranum;
#	mkdir $sranum;
#	cd $sranum;
#	nohup fastqc ../../FASTQ_Files/${i}.1_1.fastq.gz ../../FASTQ_Files/${i}.1_2.fastq.gz &
#	cd ..
#	done

#cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/SCRIPTs
	
	  
