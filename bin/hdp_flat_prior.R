options(stringsAsFactors = F)
library(hdp)
lower_threshold=0

n=as.numeric(commandArgs(T)[1])
mutations=read.table("trinuc_mut_mat_96_combined.txt", header=T,check.names =F, sep="\t",quote = "", row.names=1)
tinuc_sort <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
mutations <- mutations[tinuc_sort]

ref = read.table("sigpro_ref.txt", header = T, stringsAsFactors = F, sep = '\t')

rownames(ref) <- ref[,1]
ref <- ref[,-1]
ref <- ref[tinuc_sort,]


prior_sigs = as.matrix(ref)

nps <- ncol(prior_sigs)


ppindex <- c(1, rep(1+nps+1, nrow(mutations)))
cpindex <- c(3, rep(4, nrow(mutations)))


hdp_prior <- hdp_prior_init(prior_distn = prior_sigs,
                            prior_pseudoc = rep(1000, nps),
                            hh = rep(1, 96), # prior is uniform over 96 categories
                            alphaa = rep(1, 2), # shape hyperparameters for 2 CPs
                            alphab = rep(1, 2))  # rate hyperparameters for 2 CPs

hdp_prior <- hdp_addconparam(hdp_prior,
                              alphaa = rep(1,length(unique(cpindex))), # shape hyperparams for 2 new CPs
                              alphab = rep(1,length(unique(cpindex)))) # rate hyperparams for 2 new CPs


hdp_prior <- hdp_adddp(hdp_prior,
                       numdp = 1 + nrow(mutations),
                       ppindex = ppindex,
                       cpindex = cpindex)

hdp_prior <- hdp_setdata(hdp_prior,
                          dpindex = (1+nps+1)+1:nrow(mutations)
                          mutations)

hdp_activated <- dp_activate(hdp_prior,
                               dpindex = (1+nps+1):numdp(hdp_prior), initcc=nps+5, seed=i*300)

chain=hdp_posterior(hdp_activated,
                    burnin=30000,
                    n=100,
                    seed=n*1000,
                    space=200,
                    cpiter=3)

saveRDS(chain,paste0("hdp_chain_",n,".Rdata"))
