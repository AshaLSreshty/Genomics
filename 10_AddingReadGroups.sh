
### Script to add Read Groups to the bam file from mapping -  the WGS datasets #####s


#### Input files from mapping directory #####


#!/bin/bash

cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/MAPPING

while read a b c d e;
	do
	cd $a 
	picard AddOrReplaceReadGroups I=$a"_map.bam" O=$a"_RG.bam" RGLB=$b RGPL=$c RGPU=$d RGSM=$e
	cd ..	
	done < ReadGroups_NSC-3.txt
