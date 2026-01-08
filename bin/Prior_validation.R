#!/usr/bin/env Rscript

### Script for prior matrix validation

message(paste("Checking validating input prior matrix format.\n"))

### Load in required packages 
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(argparse))

options(stringsAsFactors = F )

### Setting up CLI argument parsing 

# Create parser
parser = ArgumentParser(prog = 'Input prior matrix check', description = 'Prior matrix validation and memory requirements generation.')
#Command line arguments
parser$add_argument("-prior_mat","--prior_matrix", type = 'character', help = "If available, specify path to prior matrix.", required=FALSE)
#parser$add_argument("-pseudo","--prior_pseudocounts", type = 'integer', nargs='+', default = "1000", help = "Specify pseudocounts weighitng for prior signatures.", required=FALSE)
message("parsing pseudocounts")
parser$add_argument("-pseudo","--prior_pseudocounts", type='character', default = "1000", help = "Specify pseudocounts weighitng for prior signatures. If mutiple, specify a list separated by commas (e.g., 1,2,3)", required=FALSE)

# Function to parse multiple pseudocounts values
pseudocount_list <- function(arg) {
  return(as.integer(unlist(strsplit(arg,","))))
}

#Parse arguments
args <- parser$parse_args()

if (!is.null(args$prior_matrix)) {
  
  prior_matrix <- args$prior_matrix
  
  if (!is.null(args$prior_pseudocounts)) {
    u_pseudocounts <- pseudocount_list(args$prior_pseudocounts)  
    #u_pseudocounts <- args$prior_pseudocounts
  }
} 

message("Importing user datasets and validating required formatting.")

##Prior matrix 
ref = read.table(prior_matrix, header = T, stringsAsFactors = F, sep = '\t')
if (ncol(ref) == 1) {
  ref <- read.table(prior_matrix, header=T, sep = ",")
  prior_signatures <- colnames(ref)
}

if ("MutationType" %in% colnames(ref)) {
  ref <- tibble::column_to_rownames(ref, "MutationType")
  prior_signatures <- colnames(ref)
} else {
  stop(sprintf("Error: Input prior matrix is incorrectly formatted (e.g., does not provide mutations under a column labelled as 'MutationType'. Please conduct the necessary data wrangling to ensure input prior matrix is compatible with the pipeline. Seek GitHub example files for reference. Stopping HDP pipeline."))
}

prior_sigs = as.matrix(ref)

prior_signatures
u_pseudocounts
message(paste0("Number of prior signatures: ",length(prior_signatures)))
message(paste0("Pseudocounts provided: ",length(u_pseudocounts)))

if (length(unique(u_pseudocounts))!=1) {
  if (length(unique(prior_signatures))==length(u_pseudocounts)) {
    message(paste0("Multiple pseudocounts values detected. Assigning ",u_pseudocounts," to prior signatures in corresponding order."))
    
  } else {
    stop(sprintf("Error: The number of pseudocounts values does not match the number of prior signatures. Please provide an equal number of pseudocounts values and prior signatures."))
  }   
} else if (length(unique(u_pseudocounts))==1) {
  message(paste0("Single pseudocounts value,", u_pseudocounts, ", provided. Assigning ",u_pseudocounts," pseudocounts to each prior signature."))
  
}
### Generate dataframe of results
message(paste0("Generating matrix of prior signatures-pseudocounts breakdown."))

priorsig_pseudocount_df <- data.frame(
  Prior_Signatures = prior_signatures,
  Pseudocounts = u_pseudocounts
)

priorsig_pseudocount_df

### Make output directory
message(paste0("Creating output directory for prior matrix validation"))  
main_dir <- getwd()
sub_dir <- paste0("validation")
if (!file.exists(sub_dir)) {
  dir.create(file.path(main_dir, sub_dir))
  u.work.dir <- file.path(main_dir,sub_dir)
  } else {
  u.work.dir <- file.path(main_dir,sub_dir)
  message(paste0("Work directory is ",u.work.dir))
  }

setwd(u.work.dir)

message(paste0("Output directory is ",u.work.dir))

### Export dataframe as CSV, saving into directory
write.table(priorsig_pseudocount_df, file = paste0(u.work.dir,"/prior_sig_pseudocounts.csv"), sep = ",",
                quote = FALSE, row.names = FALSE, col.names = TRUE)

message(paste0("Prior matrix successfully validated. Proceeding with HDP pipeline."))