
####### PREDICT SOMATIC COPY NUMBER ALTERATIONS (CNA) ####################

# 1) This script will take post quality control - base quality sequence recalibrated and duplicates marked bam files from the Post Alignment folder. The Files will be processed based on the SRA_IDs provided by the user in the INPUT_FILES directory for computing the somatic CNA (linux built) 
#
# 2) The user has to provide a csv file with SraIds. By default it will take the file in the INPUT_FILES Directory
#
# 3) The downloaded files will be saved in the "eccDNA/7_eccDNA" directory
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

# install Circle-Map package wth Biopython v1.77, samtools and bedtools

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

### Prompt to specify the path of Circle-Map and samtools executable

circlemap_exec <- readline(prompt = "Provide the path for Circle-Map package's - Circle-Map executable: ");

print(paste("The path for Circle-Map is", circlemap_exec))

samtools_exec <- readline(prompt = "Provide the path for samtools executable: ");

print(paste("The path for Samtools is", samtools_exec))

# Prompt to user to select the number of threads to use to execute samtools

print("Provide the number of threads to use for samtools")

num_threads <- readline(prompt = "Input the number of threads for samtools: ");

print(paste("Threads to take are", num_threads))


##########Creating a Dataframe to store the information of this task ############

predict_eccDNA_tbl <- data.frame("Sno" = "", "SRA_Ids" = "", "SeqTech" = "", "Input_Sorted_bamfile" = "", "Output_Predict_eccDNAs" = "", "Folder" = "", stringsAsFactors = FALSE, row.names = NULL)

### Accessing the bam from 1_DATASETS folder #########

input_dir <- "6_POST_ALIGNMENT"

output_dir <- "7_eccDNA"

### Taking the values in csv files into vectors ######

sra_Ids = as.vector(srasplit_df[,'SRA_Ids'])
class(sra_Ids)

seq_tech = as.vector(srasplit_df[,'Seq_Tech'])
class(seq_tech)

condition = as.vector(srasplit_df[,'Condition'])
class(condition)

############# Function #############

### Adding Read Groups ##########

#circlemap_exec <- ("/home/ibab/anaconda3/bin/Circle-Map")

#samtools_exec <- ("/home/ibab/PROGRAMS/GENOMICS/samtools-1.11/samtools")

###########################################

## Extracting the reads and sorting using samtools

for (i in 1:length(sra_Ids)) {

	sraId <- (sra_Ids[i])
	
	input_subdir <- paste0(input_dir, "/", sraId, sep = "")
	output_subdir <- paste0(output_dir, "/", sraId, sep = "")
	dir.create(output_subdir, showWarnings = TRUE)
	
	## The input file will be the name sorted file obtained during the Post alignment quality control process
	
	infile1 <- paste0(input_subdir, "/", sraId, "_namesorted", ".bam", sep = "")
	outfile1 <- paste0(output_subdir, "/", sraId, "_sorted_unknown_circle", ".bam", sep = "")
	
	system2(command = samtools_exec, args = c("sort", "-@", num_threads, infile1, "-o", outfile1))
	
	outfile2 <- paste0(output_subdir, "/", sraId, "_circread_candidates", ".bam", sep = "")
	
	system2(command = circlemap_exec, args = c("ReadExtractor", "-i", infile, "-o", outfile2)) # invoke command
	
	outfile3 <- paste0(output_subdir, "/", sraId, "_sort_circread_candidates", ".bam", sep = "")
	
	system2(command = samtools_exec, args = c("sort", "-@", num_threads, outfile2, "-o", outfile3))
	
	# Indexing the BAM files #
	
	system2(command = samtools_exec, args = c("index", "-@", num_threads, outfile1)) # invoke command
	
	system2(command = samtools_exec, args = c("index", "-@", num_threads, outfile3)) # invoke command
	
	# 
	
	# ref_index <- paste0("4_REF_FILES", "/", "hg38.fa", sep = "")
	ref_index <- paste0("4_REF_FILES", "/", "chr7", "/", "chr7.fa", sep = "")
	
	dir.create(

	outfile4 <- paste0(output_subdir, "/", sraId, "_circle", ".bed", sep = "")
	
	
	system2(command = circlemap_exec, args = c("Realign", "-t", num_threads, "-i", outfile3, "-qbam", infile1, "-sbam", outfile1, "-fasta", ref_index, "-o", outfile4)) # invoke command

Circle-Map Realign -t 16 -i sort_circular_read_candidates_${i}.bam -qbam qname_${i}_map.bam -sbam sorted_${i}_circle.bam -fasta ../../HGINDEX/BWA/hg38.fa -o circle_${i}.bed


	
	
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











