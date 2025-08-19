#!/usr/bin/env Rscript

### Script for HDP run with single hierarchy

application <- "HDP mutational signature extraction pipeline."

message(paste("Welcome to", application, "For any questions, please consort the original GitHub/vignette: (https://github.com/nicolaroberts/hdp/blob/master/vignettes/mutation_signatures.Rmd), or the Stratton group. \n"))

### Load in required packages 
suppressPackageStartupMessages(require(hdp))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(argparse))

options(stringsAsFactors = F)

### Setting up CLI argument parsing 
# Create parser
parser = ArgumentParser(prog = 'HDP', description='Hdp pipeline')
#Command line arguments
parser$add_argument("mutation_matrix", nargs = 1, help = "Specify path to input mutational matrix.") 

parser$add_argument("-hierarchy","--hierarchy_matrix", type = 'character', help = "If available, specify path to hierarchy matrix.", required=FALSE) 

parser$add_argument("-hp1","--hierarchy_parameter1", type = 'character', help = "Specify primary hierarchy parameter as listed in input hierarchy matrix (e.g., column name). Used to identify column.", required=FALSE)

parser$add_argument("-prior","--prior_matrix", type = 'character', help = "If available, specify path to prior matrix.", required=FALSE)

parser$add_argument("-a", "--analysis_type", type = "character", default = "Testing", help = "Specify type of analysis run. Options are [testing] or [analysis].", required=TRUE)

parser$add_argument("-b", "--burnin_iterations", type = 'double', default = "30000", help = "Specify number of burn-in iterations. Default set to 30000.", required=FALSE) 

parser$add_argument("-o", "--posterior", type = 'double', default = "100", help = "Specify number of posterior samples to collect. Default set to 100.", required=FALSE) 

parser$add_argument("-i", "--posterior_iterations", type = 'double', default = "200", help = "Specify number of iterations collected between posterior samples. Default set to 1000.", required=FALSE) 

parser$add_argument("-n", "--chain_index", type = 'character', help = "Chain index")

parser$add_argument("-t", "--threshold", type = 'character', default = "0", help = "Specify threshold for minimum mutations required. Default set to 0.")

#Parse arguments
args <- parser$parse_args()

mutation_matrix <- args$mutation_matrix
if (!exists("mutation_matrix")) {
  stop(sprintf("Mutation matrix not provided. Please specify by providing path at end of command; Use -h for further information. \n"))
}

if (!is.null("args$hierarchy_matrix")) {
  hierarchy_matrix <- args$hierarchy_matrix
}

if (!is.null("args$hierarchy_parameter1")) {
  hp1 <- args$hierarchy_parameter1
}

if (!is.null("args$prior_matrix")) {
  prior_matrix <- args$prior_matrix
} 

if (!exists("chain_index")) {
  chain_index <- args$chain_index
}

if (!exists("threshold")) {
  threshold <- args$threshold
}

lower_threshold=threshold

u_analysis_type <- args$analysis_type

if (u_analysis_type == 'analysis' | u_analysis_type == 'Analysis') {
  message(paste0("Analysis run selected. Please note that this is intended to be run on a HPC as it requires 20 threads. \n"))

  if (!is.null("args$burnin_iterations")) {
    u_burnin <- args$burnin_iterations
    }
  if (!is.null("args$posterior")) {
    u_post <- args$posterior
  }
  if (!is.null("args$posterior_iterations")) {
    u_post_space <- args$posterior_iterations
  }
} else {
  message(paste0("Testing run selected. Executing run with minimal HDP settings. \n"))
}

n <- as.numeric(chain_index)

##### Setting up HDP
message(paste("Setting up HDP posterior sampling chain ", n, " of 20. \n"))

message(paste0("Creating output subdirectory for run \n"))  
main_dir <- getwd()
sub_dir <- paste0("HDP_chains")
if (!file.exists(sub_dir)){
  dir.create(file.path(main_dir, sub_dir))
  u.work.dir <- file.path(main_dir,sub_dir)
  u.work.dir
} else {
  u.work.dir <- file.path(main_dir,sub_dir)
  message(paste0("Work directory is ",u.work.dir))
}

message("Chain ", n, ": Importing user datasets. \n")

### Import mutation matrix and conduct necessary data wrangling 
mutations=read.table(mutation_matrix, header = TRUE, check.names = FALSE, sep = "\t", quote = "", row.names = 1)

if (ncol(mutations) == 1) {
  mutations <- read.table(mutation_matrix, header = TRUE, sep = ",")
}

if ("MutationType" %in% colnames(mutation_types)) {
  mutation_types <- tibble::column_to_rownames(mutation_types, "MutationType")
} else {
  stop(sprintf("Error: Input mutation matrix does not provide mutations under a column labelled as 'MutationType'. Please conduct the necessary data wrangling to ensure your mutation matrix is compatible with the pipeline. Stopping mSigHdp pipeline."))
}

#if (ncol(mutations) > 96 | nrow(mutations) == 96) {
#  message("Input mutation matrix format detected with mutation type as rownames/in rows. Conducting data wrangling to make input matrix compatible with HDP pipeline. \n")
#  mutations <- tibble::rownames_to_column(mutations, "MutationType")
#  mutations <- t(mutations)
#  colnames(mutations) <- as.character(mutations[1, ])
#  mutations <- as.data.frame(mutations[-1,])
#}

tinuc_sort <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

mutations <- mutations[tinuc_sort]
#mutations <- as.data.frame(mutations)[, unlist(tinuc_sort)]

message(paste0("Chain ",n,": mutation matrix successfully imported. \n"))

key_table=read.table(hierarchy_matrix, header=T, check.names=FALSE, sep="\t",quote = "")
if (ncol(key_table) == 1 ) {
    key_table <- read.table(hierarchy_matrix, header=T, sep = ",")
  }
message(paste0("Chain ",n,": hierarchy matrix successfully imported. \n"))

hp1i <- which(colnames(key_table)==hp1)

freq <- table(key_table[,hp1i])

message(paste0("Chain ", n,": prior matrix provided. Extracting prior signatures to incorporate into HDP structure. \n"))

ref = read.table(prior_matrix, header = T, stringsAsFactors = F, sep = '\t')
if (ncol(ref) == 1) {
  ref <- read.table(prior_matrix, header=T, sep = ",")
}

if ("MutationType" %in% colnames(ref)) {
  ref <- tibble::column_to_rownames(ref, "MutationType")
} else {
  stop(sprintf("Error: Input prior matrix does not provide mutations under a column labelled as 'MutationType'. Please conduct the necessary data wrangling to ensure your prior matrix is compatible with the pipeline. Stopping HDP pipeline."))
}

#rownames(ref) <- ref[,1]
#ref <- ref[,-1]
#ref <- ref[tinuc_sort,]

prior_sigs = as.matrix(ref)

message(paste0("Chain ", n,": prior matrix imported and signatures extracted. Adjusting and intialising HDP structure. \n"))  

# number of prior signatures to condition on (8)
nps <- ncol(prior_sigs)

#If requiring a minimum number of mutations:
sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]

#Set Hierarchy

hdp_prior <- hdp_prior_init(prior_distn = prior_sigs, # matrix of prior sigs
                            prior_pseudoc = rep(1000, nps), # pseudocount weights
                            hh=rep(1, 96), # uniform prior over 96 categories
                            alphaa=c(1, 1), # shape hyperparams for 2 CPs
                            alphab=c(1, 1)) # rate hyperparams for 2 CPs
#hdp_prior
numdp(hdp_prior)
pseudoDP(hdp_prior)
conparam(hdp_prior)
ppindex(hdp_prior)
cpindex(hdp_prior)
dpstate(hdp_prior) # 2 for active node, 1 for frozen, 0 for heldout

# make two more CPs available for the data we will add
hdp_prior <- hdp_addconparam(hdp_prior,
                              alphaa = rep(1,length(freq)+2), # shape hyperparams for x new CPs
                              alphab = rep(1,length(freq)+2)) # rate hyperparams for x new CPs

hdp_prior <- hdp_adddp(hdp_prior,
                        numdp = nrow(mutations) + 1,
                        ppindex = c(1, rep(1+nps+1:(length(freq)), times=freq)),
                        cpindex = c(3, rep(4:(length(freq)+3), times=freq)))

# assign the data to the relevant DP nodes
hdp_prior <- hdp_setdata(hdp_prior,
                          dpindex = (1 + nps + 1) + 1:nrow(mutations), 
                          mutations) # mutation counts in all GCTs


hdp_activated <- dp_activate(hdp_prior, 
                              dpindex = (1+nps+1)+0:nrow(mutations), 
                              initcc = nps+5,
                              seed = n * 1000)

message(paste0("Chain ", n,": HDP structure initialised with priors and single hierarchy. \n"))


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