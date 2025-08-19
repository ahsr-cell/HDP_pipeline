## Introduction

![HDP pipeline overview](/assets/images/HDP_overview.png)

**HDP pipeline** is a bioinformatics pipeline for standardised mutational signature extraction using [HDP](https://github.com/nicolaroberts/hdp).

There are three main processes of the pipeline: HDP, SigProfilerPlotting, and SigProfilerAssignment. 

The pipeline executes all three processes (by default). HDP runs first, and the results it generates (i.e., a matrix containing the *de novo* signatures) are subsequently fed to [SigProfilerPlotting](https://github.com/AlexandrovLab/SigProfilerPlotting) for signature spectra plotting and [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment) for a preliminary decomposition. Depending on user needs, SigProfilerPlotting and SigProfilerAssignment can be turned off. The pipeline is designed to be compatible with the following mutation type classifications: SBS96, SBS288, SBS1536, DBS78, ID83. 


### HDP
HDP uses hierarchical Dirichlet processes to identify mutational signatures present within samples. 

The HDP 

The primary input of HDP is a `mutational matrix` (specified via /path/to/mutation_matrix), with an expected format of one row per mutation type (e.g., the 96 SBS, A[C>A]A) and one column per sample. 

HDP has two run modes (`analysis` and `testing`) set by `analysis_type`. 

**Testing** is intended for initial runs (aka "is this working" scenarios). It is run with minimal settings (1 Gibbs sampling chain using one thread, running 100 burn-in iterations, collecting 10 posterior samples off of each chain with 10 iterations between each, to allow for a short execution time. As this is for testing purposes, these run settings cannot be changed. 

**Analysis** is used for full-analysis runs, utilising 20 Gibbs sampling chains across 20 threads. These chains run 30,000 burn-in iterations and collect 100 posterior samples from each chain, with 200 iterations collected between each sample. These are standardised settings, optimised and conducted routinely. Users can change these settings by specifying `--burnin_iterations`, `--burnin_multiplier`, `--posterior`, and `--posterior_iterations`.

HDP can be run with **hierarchy** (setting `hierarchy = true`) or **no hierarchy** (setting `hierarchy=false`). In hierarchy runs, the hierarchy table is specified via `--hierarchy_matrix /path/to/hierarchy_matrix`. In its current form, the HDP pipeline has functionality for single- and double-tier hierarchy, specified by providing either one or two hierarchy parameters (for double: `--hierarchy_parameter1` and `--hierarchy_parameter2`, for single: `--hierarchy_parameter1`). The default for the pipeline is **no hierarchy**, thus providing these inputs turns the hierarchy functionality on.

HDP can be run with **prior** (setting `prior = true`, specifying the prior key via `--prior_matrix /path/to/prior_matrix`) or **no prior** (setting `prior = false`). The pipeline defaults to **no priors**, thus this functionality is turned on by providing a `prior_matrix`. 

The primary output of HDP, found under the subdirectory `/deNovo_signatures/`) include the extracted *de novo* signature table (`mSigHdp_deNovoSignatures.txt`) and, if detected, the low confidence signature table (`mSigHdp_lowConfSignatures.txt`). Standard QC plots and matrices are also included within `/deNovo_signatures/`.

### SigProfilerPlotting
SigProfilerPlotting is used to plot the extracted *de novo* signature spectra. The pipeline will feed the extracted *de novo* signature table (i.e., `mSigHdp_deNovoSignatures.txt`) and the user-specified mutational context (done via `--mutational_context SBS96`, default is SBS96) to SigProfilerPlotting for plotting. The output of SigProfilerPlotting can be found under the directory `/Signature_Spectra/`, further broken down to subdirectories `/DeNovoSignatures/` and, if detected, `/LowConfidenceSignatures/`, both containing PDFs of their plots. 

Setting `plotting` to `false` will turn off this functionality. 

### SigProfilerAssignment
SigProfilerAssignment is used to decompose the extracted *de novo* signatures. The pipeline feeds the extracted *de novo* signature table and the previous, user-provided mutational matrix (specified via `mutational matrix`) to SigProfilerAssignment, which executes its `decompose_fit()` function. The output of SigProfilerAssignment is located under the directory `/SigProfilerDecomposition/`, containing subdirectories `/Activities/`, `/Signatures/`, and `/Solution_Stats/` and decomposition plots and mappings to COSMIC signatures (e.g., [COSMIC SBS](https://cancer.sanger.ac.uk/signatures/sbs/)).

Setting `decompose` to `false` will turn off this functionality. 

## Dependencies
* Nextflow >= 24.04.2 required
* Python, required packages: [SigProfilerPlotting](), [SigProfilerAssignment](), [argparse](https://docs.python.org/3/library/argparse.html), []()
* R, required packages: [mSigHdp](https://github.com/steverozen/mSigHdp), [hdpx](https://github.com/steverozen/hdpx), [tidyverse](https://www.tidyverse.org/), [argparse](https://cran.r-project.org/web/packages/argparse/index.html)

## Installation
Clone this repository via

 > git clone git@github.com:ahsr-cell/mSigHdp_pipeline.git

## Usage

### Input files

| Input      | Description |
| ----------- | ----------- |
| `mutation_matrix`      | Required input file, provided as a tab-delimited file (.tsv). The expected format is a matrix with one row per mutation type and one column per sample. Include the mutation types under a column labelled as 'MutationType'. It is highly recommended to generate mutation matrices via [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator). For an example. please see [example_mutation_matrix.tsv](https://github.com/ahsr-cell/HDP_pipeline/blob/main/docs/example_data/example_mutation_matrix.tsv) for an example.        |
| `hierarchy`   | Required value, provided as a string. Options are `true` or `false`         |
| `hierarchy_matrix`   | Optional input file, provided if `hierarchy = true`. The expected format is a matrix with one column per sample ID (matching the input mutation matrix) and one column specifying hierarchy groupings. See [example_hierarchy_matrix.tsv](https://github.com/ahsr-cell/HDP_pipeline/blob/main/docs/example_data/example_hierarchy_matrix.tsv) for an example.         |
| `hierarchy_parameter`   | Optional value, provided as a string if `hierarchy = true`. This should be formatted exactly as the column name specifying hierarchy groupings in the input hierarchy_matrix. E.g., if a user provided [example_hierarchy_matrix.tsv](https://github.com/ahsr-cell/HDP_pipeline/blob/main/docs/example_data/example_hierarchy_matrix.tsv), `hierarchy_parameter` would be `hierarchy_parameter=Grouping`             |
| `prior_matrix`   | Optional input file, provided if `prior = true`. The expected format is a matrix with one column per mutational signature (e.g., SBS1, SBS5, etc.) and one column per mutation type (e.g., A[C>A]A). Note that the prior matrix variant class should match the mutation matrix variant class; in other words, if the inputted mutation matrix is SBS96, the prior matrix should be SBS96. See [example_prior_matrix.tsv](https://github.com/ahsr-cell/HDP_pipeline/blob/main/docs/example_data/example_prior_matrix.tsv) for an example.         |
| `analysis_type`   | Required value, provided as a string. Options are `analysis` or `testing`, default is `analysis`         |
| `mutational_context`   | Required value, provided as a string. Options are `SBS96`, `SBS288`, `SBS1536`, `DBS78`, `ID83`, default is `SBS96`         |
| `plotting`   | Required value, provided as a string. Options are `true` or `false`, default is `true`         |
| `decompose`   | Required value, provided as a string. Options are `true` or `false`, default is `true`         |
| `outdir`   | Required path, specifying the location of output files generated by pipeline. See [`output.md`](docs/output.md) for the files/directories that will be contained within this directory.        |

The pipeline can be run using:

```
nextflow run /path/to/HDP_pipeline/main.nf \
     -profile <docker/singularity/.../institute> \
     -c /path/to/config_file \
     --mutational_matrix /path/to/mutation_matrix.tsv \
     --hierarchy <true/false> \
     --hierarchy_matrix /path/to/hierarchy_key.tsv \
     --hierarchy_parameter1 <variable_used_for_primary_hierarchy_as_the_column-name_in_hierarchy_key> \
     --hierarchy_parameter2 <variable_used_for_secondary_hierarchy_as_the_column-name_in_hierarchy_key> \
     --prior_matrix /path/to/prior_matrix.tsv \
     --mutational_context <SBS96/SBS288/SBS1536/DBS78/ID83> \
     --analysis_type <analysis/testing> \
     --numchains <number_of_desired_chains>
     --outdir /path/to/outdir/ \
     --plotting <true/false> \
     --decompose <true/false> 
```

### Sanger Users
Sangers can run the pipeline using the following wrapper script:

```
#! /usr/bin/env bash

#BSUB -J USER_JOB_NAME
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -u USER@sanger.ac.uk
#BSUB -q week #set according to sample type and time requirements
#BSUB -n 20 #if testing (e.g., analysis_type=test), set to 1; otherwise, this should match the value inputted to numchains
#BSUB -M50000 #set according to sample type and memory requirements*
#BSUB -R "select[mem>50000] rusage[mem=50000]" #set according to sample type and memory requirements

#Change these options according to run requirements
mutational_matrix=/path/to/input/mutation_matrix.tsv

hierarchy=false #if running hierarchy (single-tier only), set to true

hierarchy_matrix=/path/to/input/hierarchy_key.tsv #if running no hierarchy, delete/comment this line out 
hierarchy_parameter1=sample_type #set this equal to the column name of input hierarchy_key if running no hierarchy, delete/comment this line out
hierarchy_parameter2=sample_type #set this equal to the column name of input hierarchy_key if running no hierarchy, delete/comment this line out

prior=false #if running prior, set to true

prior_matrix=/path/to/input/prior_matrix.tsv #if running no prior, delete/comment this line out

mutational_context=SBS96 #For SigProfilerPlotting, options include SBS96, SBS288, SBS1536, DBS78, ID83

analysis_type=analysis #'testing' for test run, 'analysis' for full analysis run

outdir=/path/to/desired/output/ #set to location of output files

plotting=true #set to false if you do not want to plot extracted signature spectra

decompose=true #set to false if you do not want to decompose extracted signatures 

module load cellgen/nextflow/25.04.4
module load ISG/singularity/3.11.4
config_file=/lustre/scratch125/casm/teams/team267/projects/Pipelines/HDP_pipeline/conf/sanger_lsf.config
main_script=/lustre/scratch125/casm/teams/team267/projects/Pipelines/HDP_pipeline/main.nf
profile=singularity

nextflow run ${main_script} \
     -profile ${profile} \
     -c ${config_file} \
     --mutational_matrix ${mutational_matrix} \
     --hierarchy ${hierarchy} \
     --hierarchy_matrix ${hierarchy_matrix} \
     --hierarchy_parameter1 ${hierarchy_parameter1} \
     --hierarchy_parameter2 ${hierarchy_parameter2} \
     --prior_matrix ${prior_matrix} \
     --mutational_context ${mutational_context} \
     --analysis_type ${analysis_type} \
     --outdir ${outdir} \
     --plotting ${plotting} \
     --decompose ${decompose} 
```
*Note - The pipeline resource (i.e., memory and time) requirements varies throughout the different processes. mSigHdp is the most memory and time intensive process and therefore the amount of resources requested from HPCs (for Sanger users. LSF) should be requested to ensure HDP successfully completes. In other words, if HDP requires 25Gb memory to successfully run, is recommended to request 27-30 Gb memory (25Gb + a recommended +10% additional resources requested to err on the cautious side). 

The following tables have been generated (through empirical testing, using mSigHdp (50k burn-in iterations, 250 posterior samplings, 100 iterations between each chain) on synthetic datasets of increasing sample size) as a guide for Sanger users (submitting their job to the LSF queue):

#### Normal tissue samples (max mutation burden 10,000)
| Hierarchy      | Sample size      | Max memory required [Gb] | Time required [00 hours : 00 minutes : 00 seconds] |
| ----------- | ----------- | ----------- | ----------- |
| Flat         | 25   | 6.20         | 00:18:09         |
| Flat         | 50   | 7.85         | 00:42:18         |
| Flat         | 100   | 11.28         | 01:19:10         |
| Flat         | 250   | 26.40         | 02:52:53         |
| Flat         | 500   | 47.23         | 05:59:36         |
| Flat         | 1000   | 103.00         | 13:18:21         |

#### Cancer samples (max mutation burden 50,000)
| Hierarchy      | Sample size      | Max memory required [Gb] | Time required [00 hours : 00 minutes : 00 seconds] |
| ----------- | ----------- | ----------- | ----------- |
| Flat         | 25   | 16.44         | 03:02:05         |
| Flat         | 50   | 36.55         | 05:24:18         |
| Flat         | 100   | 58.97         | 09:23:48         |
| Flat         | 250   | 160.56         | 22:21:59         |
| Flat         | 500   | 287.30         | 47:56:14         |
| Flat         | 1000   | 621.58         | 82:49:12         |

## Pipeline output

For more details about the output files and reports, please refer to the [`output.md`](docs/output.md) output documentation.

## Credits

The HDP pipeline was written by Andrew Ramos and Phuong Le. 

We thank the following people and teams for their assistance in the development of this pipeline:

Sarah Moody, Mimy Pham, Yichen Wang, and CASM IT

## Contributions and Support

Please feel free to contribute by either creating a pull request or create a new issue on this GitHub repo.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
