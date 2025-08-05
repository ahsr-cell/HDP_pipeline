#!/bin/bash -ue
Rscript --vanilla /lustre/scratch125/casm/teams/team267/projects/Pipelines/HDP_pipeline/bin/HDP_single.R --hierarchy_matrix EAC_Manuscript_sample_key.csv --hierarchy_parameter1 sample_type  --analysis_type analysis --burnin_iterations 50 --posterior 10 --posterior_iterations 200 --threshold 0 --chain_index 9 EAC_GEC_Manuscript_Subset_v1.SBS96.txt
