#!/usr/bin/env Rscript

### Script for mutation matrix validation and memory parsing

message(paste("Checking validating input mutation matrix format and generating resource requirements. \n"))

### Load in required packages 
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(argparse))

options(stringsAsFactors = F )

### Setting up CLI argument parsing 
# Create parser
parser = ArgumentParser(prog = 'Input matrix matrix check and mSigHdp memory requirements generation', description = 'Mutation matrix validation and memory requirements generation.')
#Command line arguments
parser$add_argument("matrix_path", nargs = 1, help = "Specify path to input mutational matrix.") 

#Parse arguments
args <- parser$parse_args()

mutation_matrix_path <- args$matrix_path

if(!exists("mutation_matrix_path")) {
  stop(sprintf("Mutation matrix not provided. Please specify by providing path at end of command; Use -h for further information."))
}

message("Importing user datasets and conducting necessary data wrangling.")

### Data wrangling
##Mutation matrix 
#change context to rownames
mutation_matrix <- read.csv(file = mutation_matrix_path, 
                           header = TRUE, sep="\t")

if (ncol(mutation_matrix) == 1 ) {
  stop(sprintf("Error: Incorrect format of input mutation matrix. Please input tab delimited matrix. Stopping mSigHdp pipeline."))
}

### Check that the input mutation matrix is in SigProfilerMatrixGenerator output format

if ("MutationType" %in% colnames(mutation_matrix)) {
  mutation_matrix <- tibble::column_to_rownames(mutation_matrix, "MutationType")
} else {
  stop(sprintf("Error: Input mutation matrix does not provide mutations under a column labelled as 'MutationType'. Please conduct the necessary data wrangling to ensure your mutation matrix is compatible with the pipeline. Stopping mSigHdp pipeline."))
}

### Count number of samples
#Samples will be by row
#calculate number of rows 
sample_number <- ncol(mutation_matrix)
message(paste0("Input mutation matrix has ", sample_number," samples. \n"))

### Identify mutation burden
#calculate row sums
#sum row sums
#divide by number of rows
mutation_burden <- mean(colSums(mutation_matrix))
message(paste0("Calculated mutation burden: ", round(mutation_burden, digits = 2),". \n"))

### Input into linear equations
mutation_threshold <- 10000

if (mutation_burden > mutation_threshold) {
  message(paste0("Based on mutation burden of ", round(mutation_burden, digits = 2),", input sample detected to be tumour. Calculating memory requirements based on tumour samples. \n"))
  a1 <- 0.023
  b1 <- 0.6136
  memory_required <- (a1 + b1 * sample_number)*1.10
  message(paste0("Required memory for mSigHdp run calculated to be ", round(memory_required, digits = 2)," GB, including 10% leeway. \n"))
} else {
  message(paste0("Based on mutation burden of ", round(mutation_burden, digits = 2),", input sample detected to be normal. Calculating memory requirements based on normal samples. \n"))
  a1 <- 6.27
  b1 <- 0.0832
  memory_required <- (a1 + b1 * sample_number)*1.10
  message(paste0("Required memory for mSigHdp run calculated to be ", round(memory_required, digits = 2)," GB, including 10% leeway. \n"))
}
### Make output directory
message(paste0("Creating output directory for memory requirements"))  
main_dir <- getwd()
sub_dir <- paste0("memory_requirements")
if (!file.exists(sub_dir)) {
  dir.create(file.path(main_dir, sub_dir))
  u.work.dir <- file.path(main_dir,sub_dir)
  } else {
  u.work.dir <- file.path(main_dir,sub_dir)
  message(paste0("Work directory is ",u.work.dir))
  }

setwd(u.work.dir)

message(paste0("Output directory is ",u.work.dir))

### Generate dataframe of results
memory_requirements_df <- data.frame(
  Sample_number = round(as.integer(sample_number), digits = 2),
  Mutation_burden = round(as.integer(mutation_burden),digits = 2),
  Memory_required = round(as.integer(memory_required), digits = 2)
)

### Export dataframe as CSV, saving into directory
write.table(memory_requirements_df, file = paste0(u.work.dir,"/memory_requirements.csv"), sep = ",",
                quote = FALSE, row.names = FALSE, col.names = TRUE)