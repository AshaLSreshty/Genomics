#####Script to Somatic Coy Number alterations using GATK with -  the WGS datasets #####


#### Input files from mapping directory #####


##Preparing the files for GATK Processing##

#!/bin/bash

#mkdir /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/CNAs 

cd /home/group_srivatsan02/WORK/ASHA/datasets/WGS/GBM/NSC/CNAs/CNA

sra=(SRR4009220 SRR4009231 SRR4009234 SRR4009244) 

#sra=(SRR4009299 SRR4009303 SRR4009328 SRR4009213 SRR4009217 SRR4009221 SRR4009232 SRR4009235 SRR4009288 SRR4009300 SRR4009304 SRR4009202 SRR4009215 SRR4009218 SRR4009233 SRR4009237 SRR4009298 SRR4009301 SRR4009305 SRR4009227 SRR4009228 SRR4009229 SRR4009230)

#sra=(SRR4009201 SRR4009212 SRR4009216)

for i in "${sra[@]}";
        do
        echo $i;
#	mkdir $i;
        cd $i;

#	gatk CollectReadCounts -I ../../MAPPING/$i/${i}_mrkdup_bqsr.bam -L ../preprocessed_intervals.interval_list --interval-merging-rule OVERLAPPING_ONLY -O ${i}.counts.hdf5

#	gatk DenoiseReadCounts -I ${i}.counts.hdf5 --standardized-copy-ratios ${i}_standardized_CR.tsv --denoised-copy-ratios ${i}_denoised_CR.tsv --count-panel-of-normals ../PanelOfNormals/cnv.pon.hdf5 --number-of-eigensamples 20 --annotated-intervals ../hg38_annotated.tsv

#	gatk CollectAllelicCounts -R ../../HGINDEX/BWA/hg38.fa -L ../../HGINDEX/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -I ../../MAPPING/$i/${i}_sorted_mrkdup.bam -O ${i}_GBM_allelicCounts.tsv

       gatk ModelSegments --denoised-copy-ratios ${i}_denoised_CR.tsv --allelic-counts ${i}_GBM_allelicCounts.tsv --normal-allelic-counts ../SRR4009227/SRR4009227_GBM_allelicCounts.tsv --number-of-changepoints-penalty-factor 1.0 -O ModSegs27 --output-prefix ${i}_GBM

       gatk CallCopyRatioSegments -I ModSegs27/${i}_GBM.cr.seg -O ${i}_GBM_27.called.seg










       gatk ModelSegments --denoised-copy-ratios ${i}_denoised_CR.tsv --allelic-counts ${i}_GBM_allelicCounts.tsv --normal-allelic-counts ../SRR4009228/SRR4009228_GBM_allelicCounts.tsv --number-of-changepoints-penalty-factor 1.0 -O ModSegs28 --output-prefix ${i}_GBM

        gatk CallCopyRatioSegments -I ModSegs28/${i}_GBM.cr.seg -O ${i}_GBM_28.called.seg

       gatk ModelSegments --denoised-copy-ratios ${i}_denoised_CR.tsv --allelic-counts ${i}_GBM_allelicCounts.tsv --normal-allelic-counts ../SRR4009229/SRR4009229_GBM_allelicCounts.tsv --number-of-changepoints-penalty-factor 1.0 -O ModSegs29 --output-prefix ${i}_GBM

        gatk CallCopyRatioSegments -I ModSegs29/${i}_GBM.cr.seg -O ${i}_GBM_29.called.seg

	gatk ModelSegments --denoised-copy-ratios ${i}_denoised_CR.tsv --allelic-counts ${i}_GBM_allelicCounts.tsv --normal-allelic-counts ../SRR4009230/SRR4009230_GBM_allelicCounts.tsv --number-of-changepoints-penalty-factor 1.0 -O ModSegs30 --output-prefix ${i}_GBM

	gatk CallCopyRatioSegments -I ModSegs30/${i}_GBM.cr.seg -O ${i}_GBM_30.called.seg	

	cd ..
	done








