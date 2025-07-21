#!/usr/bin/env Rscript

### Script for mSigHdp 

application <- "HDP mutational signature extraction pipeline."

message(paste("Welcome to", application, "For any questions, please consort the original GitHub/vignette: (https://github.com/nicolaroberts/hdp/blob/master/vignettes/mutation_signatures.Rmd), or the Stratton group. \n"))

### Load in required packages 
suppressPackageStartupMessages(require(hdp))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(argparse))

options(stringsAsFactors = F)
lower_threshold=0

### Setting up CLI argument parsing 
# Create parser
parser = ArgumentParser(prog = 'HDP', description='Hdp pipeline')
#Command line arguments
parser$add_argument("mutation_matrix", nargs = 1, help = "Specify path to input mutational matrix.") 

parser$add_argument("-h","--hierarchy_matrix", type = 'character', help = "If available, specify path to hierarchy matrix.", required=FALSE) 

parser$add_argument("-a", "--analysis_type", type = "character", default = "Testing", help = "Specify type of analysis run. Options are [testing] or [analysis].", required=TRUE)

parser$add_argument("-b", "--burnin_iterations", type = 'double', default = "30000", help = "Specify number of burn-in iterations. Default set to 30000.", required=FALSE) 
parser$add_argument("-o", "--posterior", type = 'double', default = "100", help = "Specify number of posterior samples to collect. Default set to 100.", required=FALSE) 
parser$add_argument("-i", "--posterior_iterations", type = 'double', default = "1000", help = "Specify number of iterations collected between posterior samples. Default set to 1000.", required=FALSE) 

parser$add_argument("-c", "--mutational_context", type = 'character', default = "SBS96", help = "Specify context of mutational matrix; options are SBS96 (default), SBS288, SBS1536, DBS78, or ID83.", required = TRUE)

parser$add_argument("-n", "--n_iter", type = 'character', default = "20", help = "n iteration, provided by for loop")

#Parse arguments
args <- parser$parse_args()

mutation_matrix <- args$mutation_matrix
if(!exists("mutation_matrix")) {
  stop(sprintf("Mutation matrix not provided. Please specify by providing path at end of command; Use -h for further information."))
}

if (!is.null("args$hierarchy_matrix")) {
  hierarchy_matrix <- args$hierarchy_matrix
}

if (!is.null("args$mutational_context")) {
  mut_context <- args$mutational_context  
}

if(!exists("mut_context")) {
  stop(sprintf("Mutational signature context not specified. Please specify using -c or --mutational_context; Use -h for further information."))
}

if(!exists("n_iter")) {
    n_iter <- args$n_iter
}

u.analysis.type <- args$analysis_type

if (u.analysis.type == 'analysis') {
  if (!is.null("args$burnin_iterations")) {
    u.burnin <- args$burnin_iterations
    }
  if (!is.null("args$posterior")) {
    u.post <- args$posterior
  }
  if (!is.null("args$posterior_iterations")) {
    u.post.space <- args$posterior_iterations
  }
}

if (mut_context == 'SBS96' | mut_context == 'SBS288' | mut_context == 'SBS1536') {
  u.mc <- 'SBS'
}
if (mut_context == 'DBS78') {
  u.mc = 'DBS'
}
  if (mut_context == 'ID83') {
    u.mc = 'ID'
}

n=as.numeric(commandArgs(T)[1])

##### Setting up HDP
message("Importing user datasets and conducting necessary data wrangling.")

### Import mutation matrix and conduct necessary data wrangling 
mutations=read.table(mutation_matrix, header=T,check.names =F, sep="\t",quote = "", row.names=1)
if (ncol(mutation_matrix) == 1 ) {
  mutations <- read.table(mutation_matrix, header=T, sep = ",")
}

mutations=read.table("trinuc_mut_mat.txt", header=T,check.names =F, sep="\t",quote = "", row.names=1)
tinuc_sort <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
mutation_table <- as.data.frame(mutation_table)[, unlist(tinuc_sort)]


key_table=read.table("key_table.txt", header=T, check.names=FALSE, sep="\t",quote = "")


ref = read.table("sigpro_ref.txt", header = T, stringsAsFactors = F, sep = '\t')

rownames(ref) <- ref[,1]
ref <- ref[,-1]
ref <- ref[tinuc_sort,]


prior_sigs = as.matrix(ref)


# number of prior signatures to condition on (8)
nps <- ncol(prior_sigs)

#If requiring a minimum number of mutations:
sample_remove=rownames(mutations)[rowSums(mutations)<lower_threshold]
mutations=mutations[!rownames(mutations)%in%sample_remove,]
key_table=key_table[!key_table$Sample%in%sample_remove,]

#Hierarchy is set per patient, can change if wanted
freq = table(key_table$Tissue)

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

chain=hdp_posterior(hdp_activated,
                    burnin=30000,
                    n=100,
                    seed=n*1000,
                    space=200,
                    cpiter=3)

saveRDS(chain,paste0("hdp_prior_chain_",n,".Rdata"))

