
####### SPLITTING FASTQ files #####################

# 1) This script will split the sra file into gzipped fastq files provided by the user  
#
# 2) The user has to provide a csv file with SRAids and the links
#
# 3) The downloaded files will be saved in the "eccDNA/1_DATASETS" directory
#
# 4) Authors: Dr. Ashalatha, Khushboo & Dr. Srivatsan

############################################################################

# Install following packages

# 1) install.packages("downloader")
# 2) install.packages("curl")

# Install following packages in linux

# 1) install sratoolkit in linux

########################################################################

# Libraries to load in R
#!/usr/bin/env Rscript

library(RCurl)
require(downloader)
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

# Prompt to specify the path of the fastq-dump executable

print("SRA Toolkit must be installed. fastq-dump will be used for converting sra file to fastq")

fastqdump_exec <- readline(prompt = "Provide the path for sratoolkit's fastq-dump executable: ");

print(paste("The path for fastq-dump is", fastqdump_exec))

##########Creating a Dataframe to store the information of this task ############

sra2fastq_tbl <- data.frame("Sno" = "", "SRA_Ids" = "", "SeqTech" = "", "fastq_1" = "", "fastq_2" = "", "Folder" = "", stringsAsFactors = FALSE, row.names = NULL)

### Accessing the sra files from 1_DATASETS folder #########

sub_dir <- "1_DATASETS"

sra_Ids = as.vector(srasplit_df[,'SRA_Ids'])
class(sra_Ids)

seq_tech = as.vector(srasplit_df[,'Seq_Tech'])
class(seq_tech)

############# Function #############

#fastqdump_exec <- ("/home/ibab/PROGRAMS/GENOMICS/sratoolkit.2.10.8-ubuntu64/bin/fasterq-dump")
#fastqdump_exec = ("/home/ibab/PROGRAMS/GENOMICS/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump")

for (i in 1:length(sra_Ids)) {
	sraId <- (sra_Ids[i])
	sra_dir <- paste0(sub_dir, "/", sraId, sep = "")
	srafile <- paste0(sra_Ids[i], ".sra", sep = "")
	sraAccession <- (paste0(sub_dir, "/", sraId, "/", srafile, sep = ""))
	dest <- paste0(sub_dir, "/", sraId, sep = "")
    	setwd(sra_dir)	
	system2(command = fastqdump_exec, args = c("--split-files", srafile, "--gzip", "--outdir", "."))  # invoke command
#	system2(command = fastqdump_exec, args = c("--split-files", srafile, "--outdir", ".", "--threads", 4))  # invoke command
	setwd(wd)
	
# Adding the information to the dataframe
	Sno <- (i)
	Seqtech <- (seq_tech[i])
	fastq_1 = paste0(sraId, "_", "1.fastq.gz", sep = "")
	fastq_2 = paste0(sraId, "_", "2.fastq.gz", sep = "")
	(exists("srafiles_df")) && is.data.frame(get("srafiles_df"))
	new_row <- c("Sno" = Sno, "SRA_Ids" = sraId, "SeqTech" = Seqtech, "fastq_1" = fastq_1, "fastq_2" = fastq_2, "Folder" = dest)
	sra2fastq_tbl <- rbind(sra2fastq_tbl, new_row)
}

### Printing the dataframe ######

print(sra2fastq_tbl)

optimal.scale <- function(w,h, wanted.w, wanted.h) max(c(w/wanted.w,  h/wanted.h))

tg = gridExtra::tableGrob(sra2fastq_tbl)
h = grid::convertHeight(sum(tg$heights), "mm", TRUE)
w = grid::convertWidth(sum(tg$widths), "mm", TRUE)
scale = optimal.scale(w,h, 279, 210) + 0.1 #A4 = 279 x 210 in landscape 
ggplot2::ggsave("12_FINAL_REPORT/Table_Details_of_Split_FASTQ_Files.pdf", tg, width = 279, height = 210, units = 'mm' , scale = scale)

############






