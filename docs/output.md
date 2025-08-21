# HDP_pipeline: Output

## Introduction

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline output overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [HDP_Chains](#HDP_Chains) - HDP run
- [HDP_ExtractedSigs](#HDP_ExtractedSigs) - HDP run
- [Signature_Spectra](#Signature_Spectra) - SigProfilerPlotting spectra
- [SigProfilerDecomposition](#SigProfilerDecomposition) - SigProfilerAssignment decomposition

### HDP output
- `HDP_Chains/`
  - `hdp_chain_1-#.Rdata`: RData files of individual HDP chains. Note that the total number of chains (and therefore files) will depend on user specifications (i.e., the value provided to `numchains`). For example, if a user specified `numchains=20`, there would be `hdp_chain_1.Rdata` to `hdp_chain_20.Rdata`.

- `HDP_ExtractedSigs/`
  - `HDP_multi_chain.Rdata`:  
  - `hdp_component_#.pdf`: Initial visualisation of extracted/de novo signatures. Note that there will be multiple files of this nature, depending on the number of signatures extracted (e.g., if 25 signatures are extracted, expect `hdp_component_0.pdf` to `hdp_component_24.pdf`).
  - `signature_attribution.pdf`
  - `QC_plots_chain.pdf`
  - `muts_attributed.pdf`
  - `mean_assignment_hdp.txt`
  - `hdp_sigs.txt`
  - `HDP_deNovoSignatures.txt`: Extracted signatures matrix (mutation types as rows and columns as extracted signatures)
  - `HDP_deNovoSigs_sigPADecomp.txt`: Extracted signatures matrix formatted to serve as input to SigProfilerAssignment

### SigProfilerPlotting output
- `Signature_Spectra/`
  - `DeNovoSignatures/`
    - `pkl/`
      - `SBS96.pkl`
    - `SBS_96_plots_deNovoSignatures.pdf`

### SigProfilerAssignment output
- `SigProfilerDecomposition/`
  - `Decompose_Solution/`
    - `Activities/`
      - `Decompose_Solution_Activities.txt`
      - `Decompose_Solution_Activity_Plots.pdf`
      - `Decompose_Solution_TMB_plot.pdf`
      - `Decomposed_MutationType_Probabilities.txt`
    - `Signatures/`
      - `Decompose_Solution_Activities.txt`
      - `SBS_96_plots_Decompose_Solution.pdf`
    - `Solution_Stats/`
      - `Cosmic_SBS96_Decomposition_Log.txt`
      - `Decompose_Solution_Samples_Stats.txt`
      - `Decompose_Solution_Signature_Assignment_log.txt`
    - `De_Novo_map_to_COSMIC_SBS96.csv`
    - `SBS96_Decomposition_Plots.pdf`     
  - `pkl/`
    - `CosmicTemplates/`
      - `COSMIC_v3.4_SBS_GRCH38.json`
    - `SBS96.pkl`
  - `JOB_METADATA_SPA.txt/`