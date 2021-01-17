
####### DOWNLOADING SRA files #####################

# 1) This script will download the sra files provided by the user  
#
# 2) The user has to provide a csv file with SRAids and the links
#
# 3) The downloaded files will be saed in the "eccDNA/1_DATASETS" directory
#
# 4) Authors: Dr. Ashalatha, Khushboo & Dr. Srivatsan

############################################################################

# Install follwoing packages

# 1) install.packages("downloader")
# 2) install.packages("curl")
# 3) install.packages("devtools")

########################################################################

# Libraries to load in R
#!/usr/bin/env Rscript

library(RCurl)
require(downloader)
library(grid)
library(gridExtra)
library(ggplot2)

#############################


# Setting the path

#setwd("/home/ibab/WORK/PROJECT/R/eccDNA_DNA/")

# Prompting the user to prepare a CSV file with S.No, SRA_Ids and the url links to download

print("Prepare a sraids_links.csv file with S.No, SRA_Ids and the urls and save it in the directory - 'INPUT_Files'")

csvfile <- readline(prompt = "Provide the name of the csv file: ");

print(paste("The file name entered is", csvfile))

###### Calling the input csv file from INPUT_Files directory #########3

sra_df <- read.csv(file.path("INPUT_FILES", csvfile), header=TRUE, sep=",")

print(sra_df)

### Creating the directories based on sra_ids #########

for (sra in sra_df$SRA_Ids) {

	if (file.exists(paste("1_DATASETS", sra))) {
		print(sra, "Directory already exists")
	} else {
		dir.create(file.path("1_DATASETS", sra))
		print("Directory created")
	}
}
	
##########Creating a Dataframe to store the information of this task ############

srafiles_df <- data.frame("Sno" = "", "SRA_Ids" = "", "SeqTech" = "", "urls" = "", "Folder" = "", stringsAsFactors = FALSE, row.names=NULL)

### Downloading the sra files based on the links provided #########

wd <- "1_DATASETS"

sra_Ids = as.vector(sra_df[,'SRA_Ids'])
class(sra_Ids)

url_vector = as.vector(sra_df[,'urls'])
class(url_vector)

seq_tech = as.vector(sra_df[,'Seq_Tech'])
class(seq_tech)

for(i in 1:length(url_vector)) {
	print(sra_Ids[i])
	print(url_vector[i])
	sraId <- (sra_Ids[i])
	sra <- paste0(sra_Ids[i], ".sra", sep = "")
	urllink <-(url_vector[i])
	dest <- (paste0(wd, "/", sraId, "/", sra, sep = ""))
	print(dest)
	try(download.file(urllink, destfile=dest, mode = "wb"))

# Adding the information to the dataframe
	Sno <- (i)
	Seqtech <- (seq_tech[i])
	(exists("srafiles_df")) && is.data.frame(get("srafiles_df"))
	new_row <- c("Sno" = Sno, "SRA_Ids" = sraId, "SeqTech" = Seqtech, "urls" = urllink, "Folder" = dest)
	srafiles_df <- rbind(srafiles_df, new_row)
}

### Printing the dataframe ######

print(srafiles_df)

optimal.scale <- function(w,h, wanted.w, wanted.h) max(c(w/wanted.w,  h/wanted.h))

tg = gridExtra::tableGrob(srafiles_df)
h = grid::convertHeight(sum(tg$heights), "mm", TRUE)
w = grid::convertWidth(sum(tg$widths), "mm", TRUE)
scale = optimal.scale(w,h, 279, 210) + 0.1 #A4 = 279 x 210 in landscape 
ggplot2::ggsave("12_FINAL_REPORT/Table_Details_of_SRA_Files.pdf", tg, width = 279, height = 210, units = 'mm' , scale = scale)

############

# download.file("https://sra-download.ncbi.nlm.nih.gov/traces/sra61/SRR/010453/SRR10704195", destfile = "1_DATASETS/SRR10704195/SRR10704195.sra", mode = "wb")

# download.file("https://sra-download.ncbi.nlm.nih.gov/traces/sra73/SRR/010453/SRR10704205", destfile = "1_DATASETS/SRR10704205/SRR10704205.sra", mode = "wb")

# download.file("https://sra-download.ncbi.nlm.nih.gov/traces/sra32/SRR/010453/SRR10704211", destfile = "1_DATASETS/SRR10704211/SRR10704211.sra", mode = "wb")

#### Prepare a dataframe of the output ############

	




