#!/usr/bin/env Rscript

### Script for HDP run with double hierarchy

message(paste("Combining HDP chains \n"))

### Load in required packages
suppressPackageStartupMessages(require(hdp))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require(argparse))

options(stringsAsFactors = FALSE)

### Setting up CLI argument parsing
# Create parser
parser = ArgumentParser(prog = 'HDP', description='Hdp pipeline')
#Command line arguments
parser$add_argument("mutation_matrix", nargs = 1, help = "Specify path to input mutational matrix.") 

#parser$add_argument("-hierarchy","--hierarchy_matrix", type = 'character', help = "If available, specify path to hierarchy matrix.", required=FALSE) 

parser$add_argument("-h_chains", "--HDP_chains", type = 'character', default = "0", help = "Specify path to HDP chains.")

parser$add_argument("-num_chains", "--number_chains", type = 'character', default = "0", help = "Specify total number of HDP chains.")

parser$add_argument("-t", "--threshold", type = 'character', default = "0", help = "Specify threshold for minimum mutations required. Default set to 0.")

parser$add_argument("-c", "--mutation_context", type = 'character', default = "SBS96", help = "Specify context of mutational matrix; options are SBS96 (default), SBS288, SBS1536, DBS78, or ID83.")

#Parse arguments
args <- parser$parse_args()

mutation_matrix <- args$mutation_matrix

#if (!is.null("args$hierarchy_matrix")) {
#  hierarchy_matrix <- args$hierarchy_matrix
#}

if (!is.null("args$HDP_chains")) {
  HDP_chain_path <- args$HDP_chains
}

if (!is.null("args$number_chains")) {
  numofchains <- args$number_chains
}

if (!exists("threshold")) {
  threshold <- args$threshold
}

if (!is.null("args$mutation_context")) {
  mut_context <- args$mutation_context
}

lower_threshold=threshold

if (mut_context == 'SBS96' | mut_context == 'SBS288' | mut_context == 'SBS1536') {
  u.mc <- 'SBS'
}
if (mut_context == 'DBS78') {
  u.mc = 'DBS'
}
if (mut_context == 'ID83') {
    u.mc = 'ID'
}

chlist <- vector("list", as.integer(numofchains))
for (i in 1:as.integer(numofchains)){
  if (file.exists(paste0(HDP_chain_path, "hdp_chain_", i, ".Rdata"))) {
    chlist[[i]] <- readRDS(paste0(HDP_chainp_path, "hdp_chain_", i, ".Rdata"))
  }
}
message(paste("Successfully imported HDP chains. Generating QC plots. \n"))

message(paste0("Creating output subdirectory for run"))  
main_dir <- getwd()
sub_dir <- paste0("HDP_ExtractedSigs")
if (!file.exists(sub_dir)){
  dir.create(file.path(main_dir, sub_dir))
  u.work.dir <- file.path(main_dir,sub_dir)
  u.work.dir
  } else {
    u.work.dir <- file.path(main_dir,sub_dir)
    message(paste0("Work directory is ",u.work.dir))
  }

if(any(unlist(lapply(chlist,is.null)))) chlist=chlist[-which(unlist(lapply(chlist,is.null)))]

mut_example_multi <- hdp_multi_chain(chlist)
pdf("QC_plots_chain.pdf")
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")
dev.off()

mut_example_multi <- hdp_extract_components(mut_example_multi) 
#This step can take a while. If too long, submit R script as job
saveRDS(mut_example_multi, "HDP_multi_chain.Rdata")

pdf("muts_attributed.pdf")
plot_comp_size(mut_example_multi, bty = "L")
dev.off()

trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))
mut_colours = c("dodgerblue", "black", "red", "grey70", "olivedrab3", "plum2")

#dev.new(width=12,height=4)
#par(mfrow=c(3,4))


for (i in 0:mut_example_multi@numcomp){
  pdf(paste0("hdp_component_", i, ".pdf"), width = 12, height = 4)

  plot_comp_distn(mut_example_multi, cat_names = trinuc_context,
                  grouping = group_factor, col = mut_colours, comp = i,
                  col_nonsig = "grey80", show_group_labels = TRUE)
  dev.off()
}

nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

plot_dp_comp_exposure(mut_example_multi,
                      dpindices=2:4, incl_numdata_plot=FALSE,
                      col=mycolors,
                      incl_nonsig=TRUE, cex.names=0.8,
                      ylab_exp = 'Signature exposure', leg.title = 'Signature')

mutations=read.table(mutation_matrix, header = TRUE, check.names = FALSE, sep = "\t", quote = "", row.names=1)
if (ncol(mutation_matrix) == 1) {
  mutations <- read.table(mutation_matrix, header = TRUE, sep = ",")
}
#key_table=read.table(hierarchy_matrix, header = TRUE, check.names = FALSE, sep = "\t", quote = "")
#if (ncol(key_table) == 1) {
#  key_table <- read.table(hierarchy_matrix, header = TRUE, sep = ",")
#}
#If requiring a minimum number of mutations:
sample_remove <- rownames(mutations)[rowSums(mutations) < lower_threshold]
mutations <- mutations[!rownames(mutations) %in% sample_remove, ]
#key_table <- key_table[!key_table$Sample %in% sample_remove, ]

#freq <- table(key_table$Patient)

#pdf("signature_attribution.pdf",width=10,height=8)
#plot_dp_comp_exposure(mut_example_multi, dpindices=(length(freq)+2):length(mut_example_multi@comp_dp_counts), incl_nonsig = T, ylab_exp = 'Signature exposure', leg.title = 'Signature', col=mycolors, incl_numdata_plot=F)
#dev.off()

message(paste("QC plots completed. Generating matrices containing raw profiles and attributions of extracted de novo signatures. \n"))
dp_distn <- comp_dp_distn(mut_example_multi)
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)
#mean_assignment <- t(dp_distn$mean[length(freq)+1+1:nrow(mutations),drop=FALSE])

mean_assignment <- as.data.frame(comp_dp_distn(mut_example_multi)$mean)
write.table(mean_assignment, "mean_assignment_hdp.txt")

mean_sigs <- as.data.frame(t(comp_categ_distn(mut_example_multi)$mean))

tinuc_sort <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G",
                "C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C",
                "T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A",
                "C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T",
                "T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G",
                "A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C",
                "G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A",
                "A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T",
                "G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G",
                "T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C",
                "C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A",
                "T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T",
                "C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G",
                "G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

rownames(mean_sigs) <- tinuc_sort

tri_nuc_cOrder <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")

mean_sigs <- as.data.frame(mean_sigs)[unlist(tri_nuc_cOrder),]

mean_sigs <- rownames_to_column(mean_sigs,"MutationType")

write.table(mean_sigs, file = paste0(u.work.dir,"/HDP_deNovoSignatures.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

mean_sigs_SigPa <- mean_sigs

mut_context <- u_mut

colnames(mean_sigs_SigPa) = c('MutationType', paste0(mut_context, LETTERS[1:ncol(mean_sigs_SigPa) - 1]))

write.table(mean_sigs_SigPa, file = paste0(u.work.dir,"/HDP_deNovoSigs_sigPADecomp.txt"), 
            quote = FALSE, row.names = FALSE, sep = '\t')

message(paste("Extracted de novo signatures matrices generation completed. Output found in /HDP_ExtractedSigs. \n"))