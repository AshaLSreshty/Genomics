
####### DOWNLOADING Human Genome Reference files #####################

# 1) This script will download the hg38 reference genome from the UCSC browser  
# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/
#
# 2) The files will be then indexed using BWA 
#
# 3) The downloaded files will be saved and processed in the "eccDNA/4_REF_FILES" directory
#
# 4) Authors: Dr. Ashalatha, Khushboo & Dr. Srivatsan

############################################################################

# Install following packages

# 1) install.packages("downloader")
# 2) install.packages("curl")
# 3) install.packages("devtools")
# 4) install.packages("pdftools")

########################################################################

# Libraries to load in R
#!/usr/bin/env Rscript

library(devtools)
library(RCurl)
require(downloader)
library(grid)
library(gridExtra)
library(ggplot2)
library(pdftools)


#############################

# Setting the path

basedir <- ("/home/ibab/WORK/PROJECT/R/eccDNA_DNA/")

setwd(basedir)

# Informing the user that hg38 reference genome from UCSC bowser will be dowbloaded

print("The human reference genome - hg38.fullAnalysisSet.chroms.tar.gz will be downloaded from the UCSC browser")

############################

wd <- "4_REF_FILES"

urllink1 = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.fullAnalysisSet.chroms.tar.gz"

dest_file2 <- (paste0(wd, "/", "hg38.fullAnalysisSet.chroms.tar.gz", sep = ""))

# Downloading the ref genome - hg38 file

download.file(urllink1, destfile = dest_file2, mode = "wb")

##################Preparing index files using BWA tool################

# Prompt to specify the path of the bwa executable

print("BWA alignment tool must be installed. bwa will be used to generate the index for the downloaded hg38 human genome reference file")

bwa_exec <- readline(prompt = "Provide the path for BWA package's - bwa executable: ");

print(paste("The path for bwa is", bwa_exec))

##################Extract the Human genome files from tar compressed file###################

untar("hg38.fullAnalysisSet.chroms.tar.gz")

# Concatenate all files

setwd("hg38.fullAnalysisSet.chroms")

system2(command = "cat", args = c("*.fa > ../hg38.fa"))

setwd(basedir)

##################Generate the index files######################

#bwa_exec <- ("/home/ibab/PROGRAMS/GENOMICS/bwa-0.7.17/bwa")

setwd(wd)

system2(command = bwa_exec, args = c("index -a bwtsw hg38.fa")) # invoke command

#################################

############ Downloading Annotation files for computing copy number alteration #########3

## All the required files will be downloaded form http://genomedata.org/pmbio-workshop/references/gatk/

wd <- "4_REF_FILES"

urllink2 = "http://genomedata.org/pmbio-workshop/references/gatk/1000G_omni2.5.hg38.vcf.gz"
urllink3 = "http://genomedata.org/pmbio-workshop/references/gatk/1000G_omni2.5.hg38.vcf.gz.tbi"
urllink4 = "http://genomedata.org/pmbio-workshop/references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
urllink5 = "http://genomedata.org/pmbio-workshop/references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
urllink6 = "http://genomedata.org/pmbio-workshop/references/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz"
urllink7 = "http://genomedata.org/pmbio-workshop/references/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
urllink8 = "http://genomedata.org/pmbio-workshop/references/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
urllink9 = "http://genomedata.org/pmbio-workshop/references/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"

dest_file2 <- (paste0(wd, "/", "1000G_omni2.5.hg38.vcf.gz", sep = ""))
dest_file2 <- (paste0(wd, "/", "1000G_omni2.5.hg38.vcf.gz.tbi", sep = ""))
dest_file2 <- (paste0(wd, "/", "Homo_sapiens_assembly38.dbsnp138.vcf.gz", sep = ""))
dest_file2 <- (paste0(wd, "/", "Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi", sep = ""))
dest_file2 <- (paste0(wd, "/", "Homo_sapiens_assembly38.known_indels.vcf.gz", sep = ""))
dest_file2 <- (paste0(wd, "/", "Homo_sapiens_assembly38.known_indels.vcf.gz.tbi", sep = ""))
dest_file2 <- (paste0(wd, "/", "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", sep = ""))
dest_file2 <- (paste0(wd, "/", "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi", sep = ""))

# Downloading the annotation files for hg38 genome

download.file(urllink2, destfile = dest_file2, mode = "wb")
download.file(urllink3, destfile = dest_file3, mode = "wb")
download.file(urllink4, destfile = dest_file4, mode = "wb")
download.file(urllink5, destfile = dest_file5, mode = "wb")
download.file(urllink6, destfile = dest_file6, mode = "wb")
download.file(urllink7, destfile = dest_file7, mode = "wb")
download.file(urllink8, destfile = dest_file8, mode = "wb")
download.file(urllink9, destfile = dest_file9, mode = "wb")

######################################




