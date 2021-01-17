
####### ALIGNMENT OF FASTQ files TO REFERENCE GENOME #####################

# 1) This script will take fastq.gz files as input for the SRA_IDs provided by the user in the INPUT_FILES directory for trimming the fastq files (linux built) 
#
# 2) The user has to provide a csv file with SraIds. By default it will take the file in the INPUT_FILES Directory
#
# 3) The downloaded files will be saved in the "eccDNA/3_TRIM_FILTER" directory
#
# 4) Authors: Dr. Ashalatha, Khushboo & Dr. Srivatsan

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

srasplit_df <- read.csv(file.path("INPUT_FILES", "sraids_split.csv"), header=TRUE, sep=",")

print("The sra files to be used for splitting: ")

print(srasplit_df["SRA_Ids"])

##########################

# Prompt to specify the path of the bwa executable

print("BWA alignment tool must be installed. bwa will be used to generate the index for the downloaded hg38 human genome reference file")

bwa_exec <- readline(prompt = "Provide the path for BWA package's - bwa executable: ");

print(paste("The path for bwa is", bwa_exec))

# Prompt to user to select the number of threads to use to execute bwa

print("Provide the number of threads to use for BWA alignment")

num_threads <- readline(prompt = "Input the number of threads for BWA: ");

print(paste("Threads to take are", num_threads))


##########Creating a Dataframe to store the information of this task ############

alignment_tbl <- data.frame("Sno" = "", "SRA_Ids" = "", "SeqTech" = "", "fastq_1" = "", "fastq_2" = "", "Folder" = "", stringsAsFactors = FALSE, row.names = NULL)

### Accessing the sra files from 1_DATASETS folder #########

input_dir <- "3_TRIM_FILTER"

output_dir <- "5_ALIGNMENT"

sra_Ids = as.vector(srasplit_df[,'SRA_Ids'])
class(sra_Ids)

seq_tech = as.vector(srasplit_df[,'Seq_Tech'])
class(seq_tech)

############# Function #############

#bwa_exec <- ("/home/ibab/PROGRAMS/GENOMICS/bwa-0.7.17/bwa")

for (i in 1:length(sra_Ids)) {

	sraId <- (sra_Ids[i])
	
	input_subdir <- paste0(input_dir, "/", sraId, sep = "")
	output_subdir <- paste0(output_dir, "/", sraId, sep = "")
	dir.create(output_subdir, showWarnings = TRUE)

	inputfile1 <- paste0(input_dir, "/", sraId, "/", sraId, "_1_val_1.fq.gz", sep = "")
	inputfile2 <- paste0(input_dir, "/", sraId, "/", sraId, "_2_val_2.fq.gz", sep = "")
#	ref_index <- paste0("4_REF_FILES", "/", "hg38.fa", sep = "")
	ref_index <- paste0("4_REF_FILES", "/", "chr7", "/", "chr7.fa", sep = "")
	sam_file <- paste0(output_subdir, "/", sraId, "_map.sam", sep = "")
	bam_file <- paste0(output_subdir, "/", sraId, "_sorted.bam", sep = "")
	threads <- paste0("-t", num_threads, sep = " ")
	
	system2(command = bwa_exec, args = c("mem", threads, ref_index, inputfile1, inputfile2, "-o", sam_file)) # invoke command
	
	asBam(sam_file, output_subdir, overwrite=TRUE, indexDestination=TRUE)
	
	# moving the bam and indexed bai files to respective directories 
	 
	file.move((paste0(output_dir, sraId, ".bam", sep = "")), (paste0(output_dir, "/", sraId, sep = "")))
	file.move((paste0(output_dir, sraId, ".bam.bai", sep = "")), (paste0(output_dir, "/", sraId, sep = "")))
	
#######################################

# Adding the information to the dataframe

	Sno <- (i)
	Seqtech <- (seq_tech[i])
	fastq_1 = paste0(sraId, "_", "_1_val_1.fq.gz", sep = "")
	fastq_2 = paste0(sraId, "_", "_2_val_2.fq.gz", sep = "")

	(exists("srafiles_df")) && is.data.frame(get("srafiles_df"))
	new_row <- c("Sno" = Sno, "SRA_Ids" = sraId, "SeqTech" = Seqtech, "fastq_1" = fastq_1, "fastq_2" = fastq_2, "Folder" = output_subdir)
	alignment_tbl <- rbind(alignment_tbl, new_row)
}


### Printing the dataframe ######

print(alignment_tbl)

optimal.scale <- function(w,h, wanted.w, wanted.h) max(c(w/wanted.w,  h/wanted.h))

tg = gridExtra::tableGrob(alignment_tbl)
h = grid::convertHeight(sum(tg$heights), "mm", TRUE)
w = grid::convertWidth(sum(tg$widths), "mm", TRUE)
scale = optimal.scale(w,h, 279, 210) + 0.1 #A4 = 279 x 210 in landscape 
ggplot2::ggsave("12_FINAL_REPORT/Table_Details_of_Reference_Alignment of fastq files.pdf", tg, width = 279, height = 210, units = 'mm' , scale = scale)

##################################









