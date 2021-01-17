
####### POST ALIGNMNET QUALITY CONTROL OF REFERERENCE ALIGNED FILES #####################

# 1) This script will take fastq.gz files as input for the SRA_IDs provided by the user in the INPUT_FILES directory for trimming the fastq files (linux built) 
#
# 2) The user has to provide a csv file with SraIds. By default it will take the file in the INPUT_FILES Directory
#
# 3) The downloaded files will be saved in the "eccDNA/3_TRIM_FILTER" directory
#
# 4) Authors: Dr. Ashalatha Sreshty, Khushboo & Dr. Srivatsan

############################################################################

# Install following packages

# 1) install.packages("downloader")
# 2) install.packages("curl")
# 3) install.packages("devtools")
# 4) install.packages("pdftools")
# 5) install.packages("filesstrings")
# 6) BiocManager::install("Rsamtools") # To convert sam to bam

########################################################################

# Libraries to load in R
#!/usr/bin/env Rscript

#library(devtools)
library(RCurl)
require(downloader)
library(grid)
library(gridExtra)
library(ggplot2)
library(pdftools)
library(Rsamtools)
library(filesstrings)

#############################

print("Setting the path if you are not in the path of main working directory")

wd = ("/home/ibab/WORK/PROJECT/R/eccDNA_DNA/")

setwd(wd)

###### Calling the input csv file from INPUT_Files directory #########

# The input sra files will be found in the 1_DATASETS

# The sra files to be split will be taken from the csvfile in the input directory -  INPUT_FILES (previously used)

srasplit_df <- read.csv(file.path("INPUT_FILES", "sraids_readgroups.csv"), header=TRUE, sep=",")

print("The sra files to be used for splitting: ")

print(srasplit_df["SRA_Ids"])

print("Taking the readgroup information")

##########################

### Prompt to specify the path of GATK and samtools executable

gatk_exec <- readline(prompt = "Provide the path for GATK package's - gatk executable: ");

print(paste("The path for GATK is", gatk_exec))

samtools_exec <- readline(prompt = "Provide the path for samtools executable: ");

print(paste("The path for Samtools is", samtools_exec))

##########Creating a Dataframe to store the information of this task ############

post_alignment_tbl <- data.frame("Sno" = "", "SRA_Ids" = "", "SeqTech" = "", "Input_Bam_file" = "", "Output_BQSR_file" = "", "Folder" = "", stringsAsFactors = FALSE, row.names = NULL)

### Accessing the bam from 1_DATASETS folder #########

input_dir <- "5_ALIGNMENT"

output_dir <- "6_POST_ALIGNMENT"

### Taking the values in csv files into vectors ######

sra_Ids = as.vector(srasplit_df[,'SRA_Ids'])
class(sra_Ids)

seq_tech = as.vector(srasplit_df[,'Seq_Tech'])
class(seq_tech)

RGLB = as.vector(srasplit_df[,'RGLB'])
class(seq_tech)

RGPL = as.vector(srasplit_df[,'RGPL'])
class(seq_tech)

RGPU = as.vector(srasplit_df[,'RGPU'])
class(seq_tech)

RGSM = as.vector(srasplit_df[,'RGSM'])
class(seq_tech)

############# Function #############

### Adding Read Groups ##########

#gatk_exec <- ("/home/ibab/PROGRAMS/GENOMICS/gatk-4.1.9.0/gatk")

for (i in 1:length(sra_Ids)) {

	sraId <- (sra_Ids[i])
	
	input_subdir <- paste0(input_dir, "/", sraId, sep = "")
	output_subdir <- paste0(output_dir, "/", sraId, sep = "")
	dir.create(output_subdir, showWarnings = TRUE)

	infile <- paste0(input_subdir, "/", sraId, ".bam", sep = "")
	outfile <- paste0(output_subdir, "/", sraId, "_RG.bam", sep = "")
	
	# Adding readgroups to the bam files
	# gatk4.1.9.0 has picard toolkit embeded so can call with gatk
	# defining the readgroups for each sra_Id
	
	RGLB <- (RGLB[i])
	RGPL <- (RGPL[i])
	RGPU <- (RGPU[i])
	RGSM <- (RGSM[i])
	
	system2(command = gatk_exec, args = c("AddOrReplaceReadGroups", "I=", infile, "O=", outfile, "RGLB=", RGLB, "RGPL=", RGPL, "RGPU=", RGPU, "RGSM=", RGSM)) # invoke command
	
}

##############

#gatk_exec <- ("/home/ibab/PROGRAMS/GENOMICS/gatk-4.1.9.0/gatk")

#samtools_exec <- ("/home/ibab/PROGRAMS/GENOMICS/samtools-1.11/samtools")

### Generating .fa.fai index and .dict dictionary files for the reference genome to be used by GATK

system2(command = samtools_exec, args = c("faidx", ref_index))
	
system2(command = gatk_exec, args = c("CreateSequenceDictionary", "-R", ref_index))

for (i in 1:length(sra_Ids)) {

	sraId <- (sra_Ids[i])
	
	wd_dir <- paste0(output_dir, "/", sraId, sep = "")
		
	infile1 <- paste0(wd_dir, "/", sraId, "_RG.bam", sep = "")
	outfile1 <- paste0(wd_dir, "/", sraId, "_namesorted.bam", sep = "")
	
	system2(command = gatk_exec, args = c("SortSam", "I=", infile1, "O=", outfile1, "SO=queryname")) # invoke command
	system2(command = samtools_exec, args = c("index", "I=", infile1)) # invoke command
	
	outfile21 <- paste0(wd_dir, "/", sraId, "_namesorted_mrkdup.bam", sep = "")
	outfile22 <- paste0(wd_dir, "/", sraId, "_metrics.txt", sep = "")
	
	system2(command = gatk_exec, args = c("MarkDuplicates", "I=", outfile1, "O=", outfile21, "ASSUME_SORT_ORDER=queryname", "METRICS_FILE=", outfile22, "QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT"))
	
	outfile3 <- paste0(wd_dir, "/", sraId, "_sorted_mrkdup.bam", sep = "")

	system2(command = gatk_exec, args = c("SortSam", "I=", outfile21, "O=", outfile3, "SO=coordinate"))
	
	system2(command = gatk_exec, args = c("BuildBamIndex", "I=", outfile3))
	
#	ref_index <- paste0("4_REF_FILES", "/", "hg38.fa", sep = "")
	ref_index <- paste0("4_REF_FILES", "/", "chr7", "/", "chr7.fa", sep = "")
	annotation_file1 <- paste0("4_REF_FILES", "/", "Homo_sapiens_assembly38.dbsnp138.vcf.gz", sep = "")
	annotation_file2 <- paste0("4_REF_FILES", "/", "Homo_sapiens_assembly38.known_indels.vcf.gz", sep = "")
	annotation_file3 <- paste0("4_REF_FILES", "/", "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", sep = "")
	
	outfile4 <- paste0(wd_dir, "/", sraId, "_sorted_mrkdup_bqsr.table", sep = "")
	
	system2(command = gatk_exec, args = c("BaseRecalibrator", "-R", ref_index, "-I", outfile3, "-O", outfile4, "--known-sites", annotation_file1, "--known-sites", annotation_file2, "--known-sites", annotation_file3, "--preserve-qscores-less-than 6 --disable-bam-index-caching" ))
	
	outfile5 <- paste0(wd_dir, "/", sraId, "_mrkdup_bqsr.bam", sep = "")

	system2(command = gatk_exec, args = c("ApplyBQSR", "-R", ref_index, "-I", outfile3, "-O", outfile5, "--bqsr-recal-file", outfile4, "--preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30"))

#######################################

# Adding the information to the dataframe

	Sno <- (i)
	Seqtech <- (seq_tech[i])
	
	bam_file <- infile
	outfile <- paste0(sraId, "_", "mrkdup_bqsr.bam", sep = "")
	
	(exists("srafiles_df")) && is.data.frame(get("srafiles_df"))
	new_row <- c("Sno" = Sno, "SRA_Ids" = sraId, "SeqTech" = Seqtech, "Input_Bam_file" = bam_file, "Output_BQSR_file" = outfile, "Folder" = wd_dir)
	post_alignment_tbl <- rbind(post_alignment_tbl, new_row)
}

### Printing the dataframe ######

print(post_alignment_tbl)

optimal.scale <- function(w,h, wanted.w, wanted.h) max(c(w/wanted.w,  h/wanted.h))

tg = gridExtra::tableGrob(post_alignment_tbl)
h = grid::convertHeight(sum(tg$heights), "mm", TRUE)
w = grid::convertWidth(sum(tg$widths), "mm", TRUE)
scale = optimal.scale(w,h, 279, 210) + 0.1 #A4 = 279 x 210 in landscape 
ggplot2::ggsave("12_FINAL_REPORT/Table_Details_of_Post_Alignment_Quality_Control.pdf", tg, width = 279, height = 210, units = 'mm' , scale = scale)

##################################









