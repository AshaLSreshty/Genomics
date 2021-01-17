
####### TRIMMING FASTQ files #####################

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
# 3) install.packages(pdftools)

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

# Prompt to specify the path of the trim_galore executable

print("Trimgalore must be installed. trim_galore will be used to trim the adapters and low quality reads of the gzipped fastq files")

trimgalore_exec <- readline(prompt = "Provide the path for TRIMGALORE package's - trim_galore executable: ");

print(paste("The path for trimgalore is", trimgalore_exec))

##########Creating a Dataframe to store the information of this task ############

trimgalore_tbl <- data.frame("Sno" = "", "SRA_Ids" = "", "SeqTech" = "", "fastq_1" = "", "fastq_2" = "", "Folder" = "", stringsAsFactors = FALSE, row.names = NULL)

### Accessing the sra files from 1_DATASETS folder #########

input_dir <- "1_DATASETS"

output_dir <- "3_TRIM_FILTER"

sra_Ids = as.vector(srasplit_df[,'SRA_Ids'])
class(sra_Ids)

seq_tech = as.vector(srasplit_df[,'Seq_Tech'])
class(seq_tech)

############# Function #############

#trimgalore_exec <- ("/home/ibab/PROGRAMS/GENOMICS/TrimGalore-0.6.6/trim_galore")

for (i in 1:length(sra_Ids)) {

	sraId <- (sra_Ids[i])
	
	input_subdir <- paste0(input_dir, "/", sraId, sep = "")
	output_subdir <- paste0(output_dir, "/", sraId, sep = "")
	dir.create(output_subdir, showWarnings = TRUE)

	inputfile1 <- paste0(input_dir, "/", sraId, "/", sraId, "_1.fastq.gz", sep = "")
	inputfile2 <- paste0(input_dir, "/", sraId, "/", sraId, "_2.fastq.gz", sep = "")
	system2(command = trimgalore_exec, args = c("-j 6 --fastqc --paired --trim-n --length 20", inputfile1, inputfile2, "-o", output_subdir)) # invoke command
	
#######################################

# Adding the information to the dataframe

	Sno <- (i)
	Seqtech <- (seq_tech[i])
	fastq_1 = paste0(sraId, "_", "1.fastq.gz", sep = "")
	fastq_2 = paste0(sraId, "_", "2.fastq.gz", sep = "")

	(exists("srafiles_df")) && is.data.frame(get("srafiles_df"))
	new_row <- c("Sno" = Sno, "SRA_Ids" = sraId, "SeqTech" = Seqtech, "fastq_1" = fastq_1, "fastq_2" = fastq_2, "Folder" = output_subdir)
	trimgalore_tbl <- rbind(trimgalore_tbl, new_row)
}


#######################################

# Combinig the output summary files to a single summary text file

listfile <- list.files(output_dir, pattern = "txt", full.names = T, recursive = TRUE)

out.file<-""
for(i in 1:length(listfile)){
	file <- read.table(listfile[i],header=TRUE, sep=";", stringsAsFactors=FALSE)
	out.file <- rbind(out.file, file)
}

write.table(out.file, file = "12_FINAL_REPORT/Summary_Statistics_Trimming_QualityControl.txt", sep="", row.names = FALSE, qmethod = "double",fileEncoding="UTF-8")

################################

# Converting the Summary text file to pdf

require(rmarkdown)
	
file_path <- ("12_FINAL_REPORT/Summary_Statistics_Trimming_QualityControl.txt")
file_out <- ("12_FINAL_REPORT/Summary_Statistics_Trimming_QualityControl.pdf")
my_text <- readLines(file_path) 
cat(my_text, sep="  \n", file = "my_text.Rmd")
render("my_text.Rmd", pdf_document(), file_out)
file.remove("my_text.Rmd") #cleanup

### Printing the dataframe ######

print(trimgalore_tbl)

optimal.scale <- function(w,h, wanted.w, wanted.h) max(c(w/wanted.w,  h/wanted.h))

tg = gridExtra::tableGrob(trimgalore_tbl)
h = grid::convertHeight(sum(tg$heights), "mm", TRUE)
w = grid::convertWidth(sum(tg$widths), "mm", TRUE)
scale = optimal.scale(w,h, 279, 210) + 0.1 #A4 = 279 x 210 in landscape 
ggplot2::ggsave("12_FINAL_REPORT/Table_Details_of_Trimming of the Fastq files.pdf", tg, width = 279, height = 210, units = 'mm' , scale = scale)

##################################

# Appending the two pdf file into one report and cleaning

library(pdftools)

pdf1 = ("12_FINAL_REPORT/Table_Details_of_Trimming of the Fastq files.pdf")
pdf2 = ("12_FINAL_REPORT/Summary_Statistics_Trimming_QualityControl.pdf")

final_report = ("12_FINAL_REPORT/Summary Statistics of Trimming for Quality Control of Fastq files.pdf")

pdf_combine(c(pdf1, pdf2), output = final_report)

####################################3

# Cleaning of files

file.remove("12_FINAL_REPORT/Table_Details_of_Trimming of the Fastq files.pdf") #cleanup
file.remove("12_FINAL_REPORT/Summary_Statistics_Trimming_QualityControl.pdf")








