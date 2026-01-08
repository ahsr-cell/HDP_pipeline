#!/usr/bin/env Rscript

### Script for HDP run with no hierarchy

application <- "HDP mutational signature extraction pipeline."

message(paste("Welcome to", application, "For any questions, please consort the original GitHub/vignette: (https://github.com/nicolaroberts/hdp/blob/master/vignettes/mutation_signatures.Rmd), or the Stratton group. \n"))

### Load in required packages
suppressPackageStartupMessages(require(hdp))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(argparse))

options(stringsAsFactors = FALSE)

### Setting up CLI argument parsing
# Create parser
parser <- ArgumentParser(prog = 'HDP', description='Hdp pipeline')
#Command line arguments
parser$add_argument("mutation_matrix", nargs = 1, help = "Specify path to input mutational matrix.") 

parser$add_argument("-prior_mat","--prior_matrix", type = 'character', help = "If available, specify path to prior matrix.", required=FALSE)

parser$add_argument("-pseudo","--prior_pseudocounts", type='character', default = "1000", help = "Specify pseudocounts weighitng for prior signatures. If mutiple, specify a list separated by commas (e.g., 1,2,3)", required=FALSE)
#parser$add_argument("-pseudo","--prior_pseudocounts", type = 'integer', nargs='+', default = "1000", help = "Specify pseudocounts weighitng for prior signatures.", required=FALSE)

parser$add_argument("-a", "--analysis_type", type = "character", default = "Testing", help = "Specify type of analysis run. Options are [testing] or [analysis].", required=TRUE)

parser$add_argument("-b", "--burnin_iterations", type = 'double', default = "30000", help = "Specify number of burn-in iterations. Default set to 30000.", required=FALSE) 

parser$add_argument("-o", "--posterior", type = 'double', default = "100", help = "Specify number of posterior samples to collect. Default set to 100.", required=FALSE) 

parser$add_argument("-i", "--posterior_iterations", type = 'double', default = "200", help = "Specify number of iterations collected between posterior samples. Default set to 1000.", required=FALSE) 

parser$add_argument("-n", "--chain_index", type = 'double', help = "Chain index")

parser$add_argument("-t", "--threshold", type = 'double', default = "0", help = "Specify threshold for minimum mutations required. Default set to 0.")

# Function to parse multiple pseudocounts values
pseudocount_list <- function(arg) {
  return(as.integer(unlist(strsplit(arg,","))))
}

#Parse arguments
args <- parser$parse_args()

mutation_matrix <- args$mutation_matrix
if(!exists("mutation_matrix")) {
  stop(sprintf("Mutation matrix not provided. Please specify by providing path at end of command; Use -h for further information."))
}

if (!is.null(args$prior_matrix)) {
  prior_matrix <- args$prior_matrix
  if (!is.null(args$prior_pseudocounts)) {
    u_pseudocounts <- pseudocount_list(args$prior_pseudocounts)
  }
}

if(!is.null(args$chain_index)) {
    chain_index <- args$chain_index
}

if(!is.null(args$threshold)) {
    threshold <- args$threshold
}

lower_threshold <- threshold

u_analysis_type <- args$analysis_type

if (u_analysis_type == 'analysis' | u_analysis_type == 'Analysis') {
  message(paste0("Analysis run selected. Please note that this is intended to be run on a HPC as it requires 20 threads. \n"))

  if (!is.null(args$burnin_iterations)) {
    u_burnin <- args$burnin_iterations
    }
  if (!is.null(args$posterior)) {
    u_post <- args$posterior
  }
  if (!is.null(args$posterior_iterations)) {
    u_post_space <- args$posterior_iterations
  }
}

n <- as.numeric(chain_index)

##### Setting up HDP
message(paste("Setting up HDP posterior sampling chain ", n, " of 20. \n"))

#message(paste0("Creating output subdirectory for run"))  
#  main_dir <- getwd()
#  sub_dir <- paste0("HDP_chains")
#  if (!file.exists(sub_dir)){
#    dir.create(file.path(main_dir, sub_dir))
#    u.work.dir <- file.path(main_dir,sub_dir)
#    u.work.dir
#  } else {
#    u.work.dir <- file.path(main_dir,sub_dir)
#    message(paste0("Work directory is ",u.work.dir))
#  }

message("Chain ", n, ": Importing user datasets. \n")

### Import mutation matrix and conduct necessary data wrangling
mutations=read.table(mutation_matrix, header = TRUE, check.names = FALSE, sep = "\t", quote = "", row.names = 1)

if (ncol(mutations) == 1) {
  mutations <- read.table(mutation_matrix, header = TRUE, sep = ",", row.names = 1)
}

#if ("MutationType" %in% colnames(mutation_types)) {
#  mutation_types <- tibble::column_to_rownames(mutation_types, "MutationType")
#} else {
#  stop(sprintf("Error: Input mutation matrix does not provide mutations under a column labelled as 'MutationType'. Please conduct the necessary data wrangling to ensure your mutation matrix is compatible with the pipeline. Stopping mSigHdp pipeline."))
#}

#if (ncol(mutations) > 96 | nrow(mutations) == 96) {
#  message("Input mutation matrix detected with rows as mutation type. Conducting data wrangling to make compatible with HDP pipeline.")
#  mutations <- tibble::rownames_to_column(mutations, "MutationType")
#  mutations <- t(mutations)
#  colnames(mutations) <- as.character(mutations[1, ])
#  mutations <- as.data.frame(mutations[-1,])
#}
mutations <- as.data.frame(t(mutations))

message("Chain ",n,": successfully imported input mutation matrix. \n")

message(paste0("Chain ", n,": prior matrix provided. Extracting prior signatures to incorporate into HDP structure. \n"))

ref = read.table(prior_matrix, header = TRUE, stringsAsFactors = FALSE, sep = '\t', row.names = 1)
if (ncol(ref) == 1 ) {
  ref <- read.table(prior_matrix, header=T, sep = ",", row.names = 1)
}

if ("MutationType" %in% colnames(ref)) {
  ref <- tibble::column_to_rownames(ref, "MutationType")
  prior_signatures <- colnames(ref)
} else {
  stop(sprintf("Error: Input prior matrix does not provide mutations under a column labelled as 'MutationType'. Please conduct the necessary data wrangling to ensure your prior matrix is compatible with the pipeline. Stopping HDP pipeline."))
}

#rownames(ref) <- ref[,1]
#ref <- ref[,-1]
#ref <- ref[tinuc_sort,]

prior_sigs <- as.matrix(ref)
      
message(paste0("Chain ", n,": prior matrix imported and signatures extracted. Adjusting and initialising HDP structure. \n"))

nps <- ncol(prior_sigs)

ppindex <- c(1, rep(1+nps+1, nrow(mutations)))
cpindex <- c(3, rep(4, nrow(mutations)))

if (length(unique(u_pseudocounts))!=1) {
  message(paste0("Multiple pseudocounts provided, assigning ",u_pseudocounts," to prior signatures", prior_signatures," in corresponding order."))
  hdp_prior <- hdp_prior_init(prior_distn = prior_sigs,
                            prior_pseudoc = as.integer(u_pseudocounts),
                            hh = rep(1, 96), # prior is uniform over 96 categories
                            alphaa = rep(1, 2), # shape hyperparameters for 2 CPs
                            alphab = rep(1, 2))  # rate hyperparameters for 2 CPs
}
if (length(unique(u_pseudocounts))==1) {
  message(paste0("Single pseudocount provided, assigning ",u_pseudocounts," pseudocounts to all prior signatures."))
  hdp_prior <- hdp_prior_init(prior_distn = prior_sigs,
                            prior_pseudoc = rep(as.integer(u_pseudocounts), nps),
                            hh = rep(1, 96), # prior is uniform over 96 categories
                            alphaa = rep(1, 2), # shape hyperparameters for 2 CPs
                            alphab = rep(1, 2))  # rate hyperparameters for 2 CPs
}

hdp_prior <- hdp_addconparam(hdp_prior,
                              alphaa = rep(1,length(unique(cpindex))), # shape hyperparams for 2 new CPs
                              alphab = rep(1,length(unique(cpindex)))) # rate hyperparams for 2 new CPs

hdp_prior <- hdp_adddp(hdp_prior,
                        numdp = 1 + nrow(mutations),
                        ppindex = ppindex,
                        cpindex = cpindex)

hdp_prior <- hdp_setdata(hdp_prior,
                        dpindex = (1+nps+1)+1:nrow(mutations),
                        mutations)

hdp_activated <- dp_activate(hdp_prior,
                              dpindex = (1+nps+1):numdp(hdp_prior), initcc=nps+5, seed=n*300)

message(paste0("Chain ", n,": HDP structure initialised with priors and no hierarchy. \n"))

if (u_analysis_type == 'analysis' | u_analysis_type == 'Analysis') {

  message(paste0("Chain ",n,": Executing posterior sampling chain ", n, " with analysis run settings: ",u_burnin," burn-in iterations, collecting ",u_post," posterior samples off each chain with ",u_post_space," iterations between each. \n"))
  
  chain=hdp_posterior(hdp_activated,
                    burnin=u_burnin,
                    n=u_post,
                    seed=n*1000,
                    space=u_post_space,
                    cpiter=3)
}

if (u_analysis_type == 'testing' | u_analysis_type == 'Testing' | u_analysis_type == 'test' | u_analysis_type == 'Test') {
  message(paste0("Executing posterior sampling chain number, ", n, ". Running with test run settings: 100 burn-in iterations, collecting 10 posterior samples off each chain with 10 iterations between each. \n"))

  chain=hdp_posterior(hdp_activated,
                    burnin=100,
                    n=10,
                    seed=n*1000,
                    space=10,
                    cpiter=3)
}

saveRDS(chain,paste0("hdp_chain_",n,".Rdata"))

message(paste0("Posterior sampling chain number", n, " completed. Successfully saved chain .Rdata in /HDP_chains for subsequent step. \n"))