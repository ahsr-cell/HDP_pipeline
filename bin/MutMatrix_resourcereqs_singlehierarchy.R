#!/usr/bin/env Rscript

### Script for mutation matrix validation and memory parsing

message(paste("Checking validating input mutation matrix format and generating resource requirements. \n"))

### Load in required packages 
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(argparse))

options(stringsAsFactors = F )

### Setting up CLI argument parsing 
# Create parser
parser = ArgumentParser(prog = 'Input matrix matrix check and HDP memory requirements generation', description = 'Mutation matrix validation and memory requirements generation.')
#Command line arguments
parser$add_argument("matrix_path", nargs = 1, help = "Specify path to input mutational matrix.") 
parser$add_argument("-hierarchy","--hierarchy_type", type = 'character', help = "If available, specify hierarchy type.", required=FALSE)
parser$add_argument("-hmat","--hierarchy_matrix", type = 'character', help = "If available, specify path to hierarchy matrix.", required=FALSE)
parser$add_argument("-hp1","--hierarchy_parameter1", type = 'character', help = "Specify primary hierarchy parameter as listed in input hierarchy matrix (e.g., column name). Used to identify 1st hierarchy column.", required=FALSE)
parser$add_argument("-umem","--user_defmemory", type = 'double', help = "If known and optimised, specify memory required for mSigHdp run.", required=FALSE)
#Parse arguments
args <- parser$parse_args()

mutation_matrix_path <- args$matrix_path

if(is.null(mutation_matrix_path)) {
  stop(sprintf("Mutation matrix not provided. Please specify by providing path at end of command; Use -h for further information."))
}

if (!is.null(args$hierarchy_matrix)) {
  hierarchy_matrix <- args$hierarchy_matrix
}

if (!is.null("args$hierarchy_parameter1")) {
  hp1 <- args$hierarchy_parameter1
}


u_mem <- args$user_defmemory

message("Importing user datasets and validating required formatting.")

##Mutation matrix 
#change context to rownames
mutation_matrix <- read.csv(file = mutation_matrix_path, 
                           header = TRUE, sep="\t")

if (ncol(mutation_matrix) == 1 ) {
  stop(sprintf("Error: Incorrect format of input mutation matrix. Please input tab delimited matrix. Stopping mSigHdp pipeline."))
} else {
  message("Mutation matrix ")
}

### Check that the input mutation matrix is in SigProfilerMatrixGenerator output format

if ("MutationType" %in% colnames(mutation_matrix)) {
  mutation_matrix <- tibble::column_to_rownames(mutation_matrix, "MutationType")
  message("Mutation matrix correctly formatted. Assessing hierarchy matrix.")
} else {
  stop(sprintf("Error: Input mutation matrix does not provide mutations under a column labelled as 'MutationType'. Please conduct the necessary data wrangling to ensure your mutation matrix is compatible with the pipeline. Stopping mSigHdp pipeline."))
}

### Check input hierarchy table is correct format, listed hierarchy parameter is correct, and data wrangling occurs
if (!is.null(hierarchy_matrix)) {
  hierarchy_key <- read.csv(file = hierarchy_matrix, 
                           header = TRUE, sep="\t")
  if (ncol(hierarchy_key) == 1 ) {
    stop(sprintf("Incorrect format of input hierarchy matrix. Please provide tab delimited matrix. Stopping mSigHdp pipeline."))
  }
  if (hp1 %in% colnames(hierarchy_key)) {
    hp1_i <- which(colnames(hierarchy_key)==hp1)
    message(paste0("Hierarchy parameter ", hp1, " detected in column ", hp1_i,". Hierarchy matrix correctly formatted. Assessing integration into mutation matrix."))
  } else {
  stop(sprintf("Error: Input hierarchy parameter not found in hierarchy matrix. Please correct provided hierarchy parameter. Stopping mSigHdp pipeline."))
  }
} else {
  message(paste("No hierarchy matrix provided, please note that mSigHdp will run assuming no hierarchy."))
}

if (is.null(u_mem) | (u_mem == 0)) {
  message(paste0("Calculating memory requirements according to input mutation matrix."))

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
} else if (!is.null(u_mem) | (u_mem != 0)) {
  sample_number <- ncol(mutation_matrix)
  mutation_burden <- mean(colSums(mutation_matrix))

  message(paste0("Memory requirements for mSigHdp run detected. Bypassing calculations and using user inputs for resource requirements"))
  memory_required <- u_mem
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

message(paste0("Mutation matrix and hierarchy key successfully validated. Memory requirements successfully calculated. Proceeding with HDP pipeline."))