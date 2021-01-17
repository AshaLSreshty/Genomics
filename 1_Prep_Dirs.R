
################ STRUCTURE OF DIRECTORIES ################################

# 1) This script will establish the architecture for the directories in the 
# user defined path.
#
# 2) This is an user interactive mode of script
#
# 3) Authors: Dr. Ashalatha, Khushboo & Dr. Srivatsan

############################################################################

#!/usr/bin/env Rscript

# Ask the user to input the path where the directories has to be generated

print("Running this script will generate the required directories inthe user defined path")

input_path <- readline(prompt = "Give the path to location: ");

print(paste("The path entered is", input_path))

# Creating the directories

mainDir <- "eccDNA_DNA"

# Check if the main working directory exists

if (file.exists(paste(input_path, mainDir, "/", sep = "/", collapse = "/"))) {
    cat(mainDir, " exists in", input_path, "and is a directory. \n The path is set to the directory \n")
    setwd(file.path(input_path, mainDir))
    print(paste(getwd()))

# Check if any file name with the main working directory exists

} else if (file.exists(paste(input_path, mainDir, sep = "/", collapse = "/"))) {
    cat(mainDir, " exists in", input_path, "but is a file \n")

# if the main working directory does not exist, creating one
   
} else {
    cat(mainDir, " does not exist in ", input_path, "\nCreating the main working directory \n")
    dir.create(file.path(input_path, mainDir))
    setwd(file.path(input_path, mainDir))
    print(paste(getwd()))
}

### Creating subdirectories ##########

## 1. For Storing datasets ###

## if (file.exists("1_DATASETS")){
##     print("1_DATASETS exists")
## } else {
##    dir.create(("1_DATASETS"), showWarnings = FALSE)
## }

## 1. For Storing datasets ###
dir.create(("1_DATASETS"), showWarnings = TRUE)

## 2. For Storing Fastqc files ###
dir.create(("2_QUALITY_CONTROL"), showWarnings = TRUE)

## 3. For Storing Trimmed AND Filtered files ###
dir.create(("3_TRIM_FILTER"), showWarnings = TRUE)

## 4. For Storing Genome Reference files ###
dir.create(("4_REF_FILES"), showWarnings = TRUE)

## 5. For Storing Ref Mapped files ###
dir.create(("5_ALIGNMENT"), showWarnings = TRUE)

## 6. For Storing Post Alignment files ###
dir.create(("6_POST_ALIGNMENT"), showWarnings = TRUE)

## 7. For Storing Processed eccDNA output files ###
dir.create(("7_eccDNA"), showWarnings = TRUE)

## 8. For Storing Processed Somatic Copy Number Alterations output files ##
dir.create(("8_SOMATIC_CNAs"), showWarnings = TRUE)

## 9. For Storing Processed Data Visualization files ##
dir.create(("9_DATA_VISUALS"), showWarnings = TRUE)

## 10. For Storing Annotation of eccDNA ##
dir.create(("10_ANNOTATIONS"), showWarnings = TRUE)

## 11. For Integrated Analysis of eccDNA ##
dir.create(("11_ANNOTATIONS"), showWarnings = TRUE)

## 12. For Final Report on eccDNA analysis ##
dir.create(("12_FINAL_REPORT"), showWarnings = TRUE)

## 13. For input raw files ##
dir.create(("INPUT_FILES"), showWarnings = TRUE)

#############Printing the directories generated ############

print("The structure of the working directories is as follows: \n");

print("1. For Storing datasets: 1_DATASETS \n");

print("2. For Storing Fastqc files: 2_QUALITY_CONTROL \n");

print("3. For Storing Trimmed AND Filtered files: 3_TRIM_FILTER \n");

print("4. For Storing Genome Reference files: 4_REF_FILES \n");

print("5. For Storing Ref Mapped files: 5_ALIGNMENT \n");

print("6. For Storing Post Alignment files: 6_POST_ALIGNMENT \n");

print("7. For Storing Processed eccDNA output files: 7_eccDNA \n");

print("8. For Storing Processed Somatic Copy Number Alterations output files: 8_SOMATIC_CNAs \n");

print("9. For Storing Processed Data Visualization files: 9_DATA_VISUALS \n");

print("10. For Storing Annotation of eccDNA: 10_ANNOTATIONS \n");

print("11. For Integrated Analysis of eccDNA: 11_ANNOTATIONS \n");

print("12. For Final Report on eccDNA analysis: 12_FINAL_REPORT \n");

print("13. For Raw input files to be provided by user: INPUT_FILES \n");

##### The output of this script is creating the main directory and the subdirectories ########## 
