
### Script to download SRA files from DRA database #####

# GBM - WGS SRA files

#!/bin/bash

cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/Datasets

#touch Downloading_status.txt

for i in $(cat ../List_SRRLinks-2.txt)
	do
	echo $i
	nohup wget -c $i &
#	echo "Downloading $i Complete" >> Downloading_status.txt
	done

cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/SCRIPTs




























