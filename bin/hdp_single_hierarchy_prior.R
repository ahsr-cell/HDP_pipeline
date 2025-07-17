options(stringsAsFactors = F)
library(hdp)
lower_threshold=0

n=as.numeric(commandArgs(T)[1])

mutations=read.table("trinuc_mut_mat.txt", header=T,check.names =F, sep="\t",quote = "", row.names=1)
tinuc_sort <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
mutations <- mutations[tinuc_sort]

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

