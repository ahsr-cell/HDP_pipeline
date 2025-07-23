#!/usr/bin/env Rscript

### Script for HDP run with double hierarchy

application <- "HDP mutational signature extraction pipeline."

message(paste("Welcome to", application, "For any questions, please consort the original GitHub/vignette: (https://github.com/nicolaroberts/hdp/blob/master/vignettes/mutation_signatures.Rmd), or the Stratton group. \n"))

### Load in required packages
suppressPackageStartupMessages(require(hdp))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(argparse))

options(stringsAsFactors = FALSE)

### Setting up CLI argument parsing
# Create parser
parser = ArgumentParser(prog = 'HDP', description='Hdp pipeline')
#Command line arguments
parser$add_argument("mutation_matrix", nargs = 1, help = "Specify path to input mutational matrix.") 

parser$add_argument("-hierarch","--hierarchy_matrix", type = 'character', help = "If available, specify path to hierarchy matrix.", required=FALSE) 

parser$add_argument("-prior","--prior_matrix", type = 'character', help = "If available, specify path to prior matrix.", required=FALSE)

parser$add_argument("-a", "--analysis_type", type = "character", default = "Testing", help = "Specify type of analysis run. Options are [testing] or [analysis].", required=TRUE)

parser$add_argument("-b", "--burnin_iterations", type = 'double', default = "30000", help = "Specify number of burn-in iterations. Default set to 30000.", required=FALSE) 

parser$add_argument("-o", "--posterior", type = 'double', default = "100", help = "Specify number of posterior samples to collect. Default set to 100.", required=FALSE) 

parser$add_argument("-i", "--posterior_iterations", type = 'double', default = "1000", help = "Specify number of iterations collected between posterior samples. Default set to 1000.", required=FALSE) 

parser$add_argument("-c", "--mutational_context", type = 'character', default = "SBS96", help = "Specify context of mutational matrix; options are SBS96 (default), SBS288, SBS1536, DBS78, or ID83.", required = TRUE)

parser$add_argument("-n", "--n_iter", type = 'character', default = "20", help = "n iteration, provided by for loop")

parser$add_argument("-t", "--threshold", type = 'character', default = "0", help = "Specify threshold for minimum mutations required. Default set to 0.")

#Parse arguments
args <- parser$parse_args()

mutation_matrix <- args$mutation_matrix
if (!exists("mutation_matrix")) {
  stop(sprintf("Mutation matrix not provided. Please specify by providing path at end of command; Use -h for further information."))
}

if (!is.null("args$hierarchy_matrix")) {
  hierarchy_matrix <- args$hierarchy_matrix
}

if (!is.null("args$prior_matrix")) {
  prior_matrix <- args$prior_matrix
}

if (!is.null("args$mutational_context")) {
  mut_context <- args$mutational_context
}

if (!exists("mut_context")) {
  stop(sprintf("Mutational signature context not specified. Please specify using -c or --mutational_context; Use -h for further information."))
}

if (!exists("n_iter")) {
  n_iter <- args$n_iter
}

if (!exists("threshold")) {
  threshold <- args$threshold
}

lower_threshold <- threshold

u_analysis_type <- args$analysis_type

if (u_analysis_type == 'analysis' | u_analysis_type == 'Analysis') {
  if (!is.null("args$burnin_iterations")) {
    u_burnin <- args$burnin_iterations
    }
  if (!is.null("args$posterior")) {
    u_post <- args$posterior
  }
  if (!is.null("args$posterior_iterations")) {
    u_post_space <- args$posterior_iterations
  }
}

if (mut_context == 'SBS96' | mut_context == 'SBS288' | mut_context == 'SBS1536') {
  u.mc <- 'SBS'
}
if (mut_context == 'DBS78') {
  u.mc <- 'DBS'
}
if (mut_context == 'ID83') {
    u.mc <- 'ID'
}

n <- as.numeric(n_iter)

##### Setting up HDP
message("Importing user datasets and conducting necessary data wrangling. \n")

### Import mutation matrix and conduct necessary data wrangling
if (ncol(mutations) == 1) {
  mutations <- read.table(mutation_matrix, header = TRUE, sep = ",")
}

if (ncol(mutations) > 96 | nrow(mutations) == 96) {
  message("Input mutation matrix provided detected to be formatted with rows as mutation type. Conducting data wrangling to make compatible with HDP pipeline.")
  mutations <- mutations <- tibble::rownames_to_column(mutations, "MutationType")
  mutations <- t(mutations)
  colnames(mutations) <- as.character(mutations[1, ])
  mutations <- as.data.frame(mutations[-1,])
}

tinuc_sort <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

mutations <- mutations[tinuc_sort]
#mutations <- as.data.frame(mutations)[, unlist(tinuc_sort)]

key_table=read.table(hierarchy_matrix, header=T, check.names=FALSE, sep="\t",quote = "")
if (ncol(key_table) == 1 ) {
  key_table <- read.table(hierarchy_matrix, header=T, sep = ",")
}
#If requiring a minimum number of mutations:
sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]

key_table$Type <- factor(key_table$Patient)

key_table$Type = paste0(key_table$Patient, key_table$Tissue)

key_table$GP <- as.numeric(factor(key_table$Tissue, levels = unique(key_table$Tissue)))
table(key_table$GP)

key_table$CD <- as.numeric(factor(key_table$Type, levels = unique(key_table$Type)))
table(key_table$CD)

gp <- key_table$GP
pp <- key_table %>% select(Type, GP) %>% distinct() %>% .$GP
cd <- key_table$CD

dps_to_add <- c(0, 
                rep(1, max(gp)), 
                pp + 1, 
                cd + max(gp) + 1)

#loading priors list
if (exists("prior_matrix")) {
  ref = read.table(prior_matrix, header = T, stringsAsFactors = F, sep = '\t')
  rownames(ref) <- ref[,1]
  ref <- ref[,-1]
  ref <- ref[tinuc_sort,]

  prior_sigs = as.matrix(ref)
  # number of prior signatures to condition on (8)
  nps <- ncol(prior_sigs)

  #with PID and Method as parents (2 hierarchy)
  hdp_PD_tissue_prior <- hdp_prior_init(prior_distn = prior_sigs, # matrix of prior sigs
                                          prior_pseudoc = rep(1000, nps), # pseudocount weights
                                          hh=rep(1, 96), # uniform prior over 96 categories
                                          alphaa=c(1, 1), # shape hyperparams for 2 CPs
                                          alphab=c(1, 1)) # rate hyperparams for 2 CPs

  hdp_PD_tissue_prior <- hdp_addconparam(hdp_PD_tissue_prior,
                                           alphaa = rep(1,length(unique(dps_to_add))), # shape hyperparams for 2 new CPs
                                           alphab = rep(1,length(unique(dps_to_add)))) # rate hyperparams for 2 new CPs

  #Adjustment for the priors
  pd <- c(1, dps_to_add[-1]+nps+1)

  hdp_PD_tissue_prior <- hdp_adddp(hdp_PD_tissue_prior,
                                     numdp = length(dps_to_add),
                                     ppindex = pd,
                                     cpindex = dps_to_add+3)

  # assign data to the relevant DP nodes (samples to appropriate tissues)
  hdp_PD_tissue_prior <- hdp_setdata(hdp_PD_tissue_prior,
                                       dpindex = (((length(dps_to_add) - nrow(mutations))+nps+1)+1):length(dpstate(hdp_PD_tissue_prior)),
                                       mutations)



  hdp_PD_tissue_prior_activated <- dp_activate(hdp_PD_tissue_prior,
                                                 dpindex = ((1+nps+1)+0):length(dpstate(hdp_PD_tissue_prior)),
                                                 initcc = 10,
                                                 seed = n*300)
} else {
  hdp_PD_tissue <- hdp_init(ppindex = dps_to_add, # index of parental node
                          cpindex = dps_to_add +1, # index of the CP to use
                          hh = rep(1, 96), # prior is uniform over 96 categories
                          alphaa = rep(1,length(unique(dps_to_add))), # shape hyperparameters for 2 CPs
                          alphab = rep(1,length(unique(dps_to_add))))  # rate hyperparameters for 2 CPs
  ppindex(hdp_PD_tissue)

  hdp_PD_tissue <- hdp_setdata(hdp_PD_tissue,
                             dpindex = (1+  length(unique(key_table$GP)) + length(unique(key_table$CD)) +1 ):numdp(hdp_PD_tissue), # index of nodes to add data to
                             mutations)

  hdp_act_PD_tissue <- dp_activate(hdp_PD_tissue, 1:numdp(hdp_PD_tissue), initcc=10,seed=n*300)
}

if (u_analysis_type == 'analysis' | u_analysis_type == 'Analysis') {
    chain_PD_tissue=hdp_posterior(hdp_PD_tissue_prior_activated,
                              burnin=u_burnin,
                              n=u_post,
                              seed=n*1000,
                              space=u_post_space,
                              cpiter=3)
}

if (u_analysis_type == 'testing' | u_analysis_type == 'Testing' | u_analysis_type == 'test' | u_analysis_type == 'Test') {
    chain_PD_tissue=hdp_posterior(hdp_PD_tissue_prior_activated,
                              burnin=100,
                              n=10,
                              seed=n*1000,
                              space=10,
                              cpiter=3)
}

saveRDS(chain_PD_tissue,paste0("hdp_chain_",n,"_PD_tissue_prior.Rdata"))  
