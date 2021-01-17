
####### PREDICT SOMATIC COPY NUMBER ALTERATIONS (CNA) ####################

# 1) This script will take post quality control - base quality sequence recalibrated and duplicates marked bam files from the Post Alignment folder. The Files will be processed based on the SRA_IDs provided by the user in the INPUT_FILES directory for computing the somatic CNA (linux built) 
#
# 2) The user has to provide a csv file with SraIds. By default it will take the file in the INPUT_FILES Directory
#
# 3) The downloaded files will be saved in the "eccDNA/8_SOMATIC_CNAs" directory
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

##########################

### Prompt to specify the path of GATK and samtools executable

gatk_exec <- readline(prompt = "Provide the path for GATK package's - gatk executable: ");

print(paste("The path for GATK is", gatk_exec))

##########Creating a Dataframe to store the information of this task ############

predict_cna_tbl <- data.frame("Sno" = "", "SRA_Ids" = "", "SeqTech" = "", "Input_BQSR_bamfile" = "", "Output_Predict_CNAs" = "", "Folder" = "", stringsAsFactors = FALSE, row.names = NULL)

### Accessing the bam from 1_DATASETS folder #########

input_dir <- "6_POST_ALIGNMENT"

output_dir <- "8_SOMATIC_CNAs"

### Taking the values in csv files into vectors ######

sra_Ids = as.vector(srasplit_df[,'SRA_Ids'])
class(sra_Ids)

seq_tech = as.vector(srasplit_df[,'Seq_Tech'])
class(seq_tech)

condition = as.vector(srasplit_df[,'Condition'])
class(condition)

############# Function #############

### Adding Read Groups ##########

#gatk_exec <- ("/home/ibab/PROGRAMS/GENOMICS/gatk-4.1.9.0/gatk")

###########################################

### Generating the interval list file for the reference genome to be used for read counts and in further processing

#index_file <- paste0("4_REF_FILES", "/", "hg38.fa", sep = "")
index_file <- paste0("4_REF_FILES", "/", "chr7", "/", "chr7.fa", sep = "")

outref_dir <- paste0("8_SOMATIC_CNAs", "/", "Reference", sep = "")
dir.create(outref_dir, showWarnings = TRUE)

outfile_intervals <- paste0(outref_dir, "/", "preprocessed_intervals.interval_list", sep = "")

system2(command = gatk_exec, args = c("PreprocessIntervals", "-R", index_file, "--bin-length 1000 --padding 0", "-O", outfile_intervals))

print("preprocessed_intervals.interval_list for the index is generated")

##############################

### Generating Readcounts for the aligned, dorted and duplicate marked, BQSR bamfile

for (i in 1:length(sra_Ids)) {

	sraId <- (sra_Ids[i])
	
	input_subdir <- paste0(input_dir, "/", sraId, sep = "")
	output_subdir <- paste0(output_dir, "/", sraId, sep = "")
	dir.create(output_subdir, showWarnings = TRUE)
	
	infile <- paste0(input_subdir, "/", sraId, "_mrkdup_bqsr.bam", sep = "")
	outfile <- paste0(output_subdir, "/", sraId, ".counts.hdf5", sep = "")
	preprocess_intervals <- paste0(output_dir, "/", "Reference", "/", "preprocessed_intervals.interval_list", sep = "")
	
	system2(command = gatk_exec, args = c("CollectReadCounts", "-I", infile, "-L", preprocess_intervals, "--interval-merging-rule OVERLAPPING_ONLY", "-O", outfile)) # invoke command

}
####################Creating Panel of Normals#################

normals <- c() # creating a vector containing normal sample Ids

for (i in 1:length(sra_Ids)) {
	
	if (condition[i] == 'Normal') {
		normals <- c(normals, sra_Ids[i])
	} else {
		print(paste0(sra_Ids[i], " is not a normal"))
	}
}


#### Generating a command to input for the create panel of normals step ######

dir.create(paste0(output_dir, "/", "PanelofNormals"))

pon_path <- paste0(output_dir, "/", "PanelofNormals")

files.to.copy <- paste0(output_dir, "/", normals, "/", normals, ".counts.hdf5", sep = "")

file.copy(files.to.copy, pon_path)

### Create Panel of Normals ######

in_pon <- NULL; # creating a null string

list_ponfiles <- list.files(pon_path)

in_pon <- paste0("-I ", pon_path, "/", list_ponfiles, collapse = " ") # making a string with the ids in the normal vector generated from the above command

#print(in_pon)

annot_ref <- paste0(output_dir, "/", "hg38_annotated.tsv", sep="")

outfile1 <- paste0(output_dir, "/", "PanelofNormals", "/", "cnv.pon.hdf5", sep = "")

system2(command = gatk_exec, args=c("CreateReadCountPanelOfNormals", in_pon, "-O", outfile1))

################################

for (i in 1:length(sra_Ids)) {

	sraId <- (sra_Ids[i])
	
	input_subdir <- paste0(input_dir, "/", sraId, sep = "")
	output_subdir <- paste0(output_dir, "/", sraId, sep = "")
	
	infile2 <- paste0(input_subdir, "/", sraId, "_mrkdup_bqsr.bam", sep = "")
	outfile2 <- paste0(output_subdir, "/", sraId, ".counts.hdf5", sep = "")
	interval_list <- paste0(output_dir, "/", "Reference", "/", "preprocessed_intervals.interval_list", sep = "")
	
	system2(command = gatk_exec, args = c("CollectReadCounts", "-I", infile2, "-L", interval_list, "--interval-merging-rule OVERLAPPING_ONLY", "-O", outfile2)) # invoke command
	
	outfile3 <- paste0(output_subdir, "/", sraId, "_standardized_CR.tsv", sep = "")
	outfile4 <- paste0(output_subdir, "/", sraId, "_denoised_CR.tsv", sep = "")

	system2(command = gatk_exec, args = c("DenoiseReadCounts", "-I", outfile2, "--standardized-copy-ratios", outfile3, "--denoised-copy-ratios", outfile4, "--count-panel-of-normals", outfile1, "--number-of-eigensamples 20")) # invoke command

	site_file <- paste0("4_REF_FILES", "/", "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", sep = "")
	infile3 <- paste0(input_subdir, "/", sraId, "_sorted_mrkdup.bam", sep = "")
	outfile5 <- paste0(output_subdir, "/", sraId, "_GBM_allelicCounts.tsv", sep = "")

	system2(command = gatk_exec, args = c("CollectAllelicCounts", "-R", index_file, "-L", site_file, "-I", infile3, "-O", outfile5)) # invoke command
	
	infile4 <- paste0(output_subdir, "/", sraId, "_GBM_allelicCounts.tsv", sep = "")
	
	outfile6 <- paste0(sraId, "_GBM", sep = "")
	
	system2(command = gatk_exec, args = c("ModelSegments", "--denoised-copy-ratios", outfile4, "--allelic-counts", outfile5, "--normal-allelic-counts", infile4, "--number-of-changepoints-penalty-factor 1.0", "-O", "ModSegs", "--output-prefix", outfile6 )) # invoke command
	
	infile5 <- paste0(output_subdir, "/", "ModSegs", "/", sraId, "_GBM.cr.seg", sep = "")
	outfile7 <- paste0(output_subdir, "/", sraId, "_GBM.called.seg", sep = "")
	
	system2(command = gatk_exec, args = c("CallCopyRatioSegments", "-I", infile5, "-O", outfile7))
	
#######################################

# Adding the information to the dataframe

	Sno <- (i)
	Seqtech <- (seq_tech[i])
	
	BQSR_bamfile <- infile2
	CNA_Pred <- outfile7
		
	(exists("srafiles_df")) && is.data.frame(get("srafiles_df"))
	new_row <- c("Sno" = Sno, "SRA_Ids" = sraId, "SeqTech" = Seqtech, "Input_BQSR_bamfile" = BQSR_bamfile, "Output_Predict_CNAs" = CNA_Pred, "Folder" = output_dir)
	predict_cna_tbl <- rbind(predict_cna_tbl, new_row)
}

### Printing the dataframe ######

print(predict_cna_tbl)

optimal.scale <- function(w,h, wanted.w, wanted.h) max(c(w/wanted.w,  h/wanted.h))

tg = gridExtra::tableGrob(predict_cna_tbl)
h = grid::convertHeight(sum(tg$heights), "mm", TRUE)
w = grid::convertWidth(sum(tg$widths), "mm", TRUE)
scale = optimal.scale(w,h, 279, 210) + 0.1 #A4 = 279 x 210 in landscape 
ggplot2::ggsave("12_FINAL_REPORT/Table_Details_of_Predicted_Somatic_Copy_Number_Alterations.pdf", tg, width = 279, height = 210, units = 'mm' , scale = scale)

##################################











