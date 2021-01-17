
####### SPLITTING FASTQ files #####################

# 1) This script will take fastq.gz files as input for the SRA_IDs provided by the user in the INPUT_FILES directory for quality control using fastqc (linux built) 
#
# 2) The user has to provide a csv file with SraIds. By default it will take the file in the INPUT_FILES Directory

#
# 3) The downloaded files will be saved in the "eccDNA/2_QUALITY_CONTROL" directory
#
# 4) Authors: Dr. Ashalatha, Khushboo & Dr. Srivatsan

############################################################################

# Install following packages

# 1) install.packages("downloader")
# 2) install.packages("curl")
# 3) install.packages("fastqcr") # for making consolidated report

# Install following packages in linux

# 1) install fastqcr

########################################################################

# Libraries to load in R
#!/usr/bin/env Rscript

library(dplyr)
library(RCurl)
library(fastqcr)
library(grid)
library(gridExtra)
library(ggplot2)

#############################

print("Setting the path if you are not in the path of main working directory")

wd = ("/home/ibab/WORK/PROJECT/R/eccDNA_DNA/")

setwd(wd)

###### Calling the input csv file from INPUT_Files directory #########

# The input sra files will be found in the 1_DATASETS

# The sra files to be split will be taken from the csvfile in the input directory -  INPUT_FILES (previously used)

srasplit_df <- read.csv(file.path("INPUT_FILES", "sraids_split.csv"), header=TRUE, sep=",")

print("The sra files to be used for splitting: \n")

print(srasplit_df["SRA_Ids"])

##########################

# Prompt to specify the path of the fastqc executable

print("FASTQC must be installed. fastqc will be used for performing quality control of the gzipped fastq files")

fastqc_exec <- readline(prompt = "Provide the path for FASTQC package's - fastqc executable: ");

print(paste("The path for fastqc is", fastqc_exec))

##########Creating a Dataframe to store the information of this task ############

fastqc_tbl <- data.frame("Sno" = "", "SRA_Ids" = "", "SeqTech" = "", "fastq_1" = "", "fastq_2" = "", "Folder" = "", stringsAsFactors = FALSE, row.names = NULL)

### Accessing the sra files from 1_DATASETS folder #########

input_dir <- "1_DATASETS"

output_dir <- "2_QUALITY_CONTROL"

sra_Ids = as.vector(srasplit_df[,'SRA_Ids'])
class(sra_Ids)

seq_tech = as.vector(srasplit_df[,'Seq_Tech'])
class(seq_tech)

############# Function #############

fastqc_exec <- ("/home/ibab/PROGRAMS/GENOMICS/FastQC/fastqc")

for (i in 1:length(sra_Ids)) {

	sraId <- (sra_Ids[i])
	
	input_subdir <- paste0(input_dir, "/", sraId, sep = "")
	output_subdir <- paste0(output_dir, "/", sraId, sep = "")
	dir.create(output_subdir, showWarnings = TRUE)
	
	inputfile1 <- paste0(input_dir, "/", sraId, "/", sraId, "_1.fastq.gz", sep = "")
	inputfile2 <- paste0(input_dir, "/", sraId, "/", sraId, "_2.fastq.gz", sep = "")

	system2(command = fastqc_exec, args = c(inputfile1, inputfile2, "-o", output_subdir))  # invoke command

#######################################

# Adding the information to the dataframe

	Sno <- (i)
	Seqtech <- (seq_tech[i])
	fastq_1 = paste0(sraId, "_", "1.fastq.gz", sep = "")
	fastq_2 = paste0(sraId, "_", "2.fastq.gz", sep = "")
	
	(exists("srafiles_df")) && is.data.frame(get("srafiles_df"))
	new_row <- c("Sno" = Sno, "SRA_Ids" = sraId, "SeqTech" = Seqtech, "fastq_1" = fastq_1, "fastq_2" = fastq_2, "Folder" = output_subdir)
	fastqc_tbl <- rbind(fastqc_tbl, new_row)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# making a copy of files to temp_dir

temp_dir <- paste0(output_dir, "/", "temp_fastq", sep = "")
dir.create(temp_dir, showWarnings = TRUE)

for (i in 1:length(sra_Ids)) {

	sraId <- (sra_Ids[i])
	
	output_subdir <- (paste0(output_dir, "/", sraId, sep = ""))

# Transfering all files to a temp directory to make a multi-QC report

	# find the files that you want
	list.of.files <- list.files(output_subdir, full.names = TRUE, recursive = TRUE)

	# copy the files to the new folder
	file.copy(list.of.files, temp_dir)	

}

# Inspecting QC Problems
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Building Multi QC Reports
fastqc_dir <- (paste0(output_dir, "/", "temp_fastq", sep = ""))

qc <- qc_aggregate(fastqc_dir)

report <- (paste0("12_FINAL_REPORT", "/", "multi-qc-report", sep = ""))
	
qc_report(fastqc_dir, result.file = report, experiment = "WGS Sequencing of Glioblastoma patients")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


### Printing the dataframe ######

print(fastqc_tbl)

optimal.scale <- function(w,h, wanted.w, wanted.h) max(c(w/wanted.w,  h/wanted.h))

tg = gridExtra::tableGrob(fastqc_tbl)
h = grid::convertHeight(sum(tg$heights), "mm", TRUE)
w = grid::convertWidth(sum(tg$widths), "mm", TRUE)
scale = optimal.scale(w,h, 279, 210) + 0.1 #A4 = 279 x 210 in landscape 
ggplot2::ggsave("12_FINAL_REPORT/Table_Details_of_Quality_Control.pdf", tg, width = 279, height = 210, units = 'mm' , scale = scale)

############






