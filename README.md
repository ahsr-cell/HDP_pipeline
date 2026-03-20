## Introduction

![HDP pipeline overview](/docs/images/HDP_overview.png)

**HDP pipeline** is a bioinformatics pipeline for standardised mutational signature extraction using [HDP](https://github.com/nicolaroberts/hdp). This pipeline was developed in conjunction with the [mSigHdp pipeline](https://github.com/ahsr-cell/mSigHdp_pipeline), another mutation signature extraction pipeline that may be of interests to users. 

There are four main processes of the pipeline: Validation, HDP, SigProfilerPlotting, and SigProfilerAssignment. 

The pipeline executes all four processes (by default). Validation runs first, followed by HDP. The results of HDP (i.e., a matrix containing the *de novo* signatures) are subsequently fed to [SigProfilerPlotting](https://github.com/AlexandrovLab/SigProfilerPlotting) for signature spectra plotting and [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment) for a preliminary decomposition. Validation and HDP will always run, but depending on user needs, SigProfilerPlotting and SigProfilerAssignment can be turned off. The pipeline is compatible with the following mutation type classifications: SBS96, DBS78, ID83. 

### Validation
The validation step of the pipeline serves two purposes. First, it checks the input matrices' format(s), (e.g., for the input mutational matrix, is it in the standardised SigProfiler format: mutation types as rows, samples as columns). If they are not in the required format, the pipeline will stop and print out (in error logs), which input is incorrectly formatted. See [Inputs](https://github.com/ahsr-cell/HDP_pipeline?tab=readme-ov-file#input-files) for more information on required formatting.

Second, the pipeline provides functionality to automatically calculate resource requirements for HDP runs. HDP memory and runtime requirements scale with sample number and maximum mutation burden. Within the HDP pipeline, linear mixed-effects models are used to approximate memory requirements, accounting for sample number and mutation burden. The pipeline take the inputted mutation matrix, calculate these two metrics, use the linear mixed-effect model to calculate the memory requirement, and finally pass this value for each HDP chain. 

### HDP
HDP uses hierarchical Dirichlet processes to identify mutational signatures present within samples. 

The primary input of HDP is a `mutational matrix` (specified via /path/to/mutation_matrix), with an expected format of one row per mutation type (e.g., the 96 SBS, A[C>A]A) and one column per sample. 

HDP has two run modes (`analysis` and `testing`) set by `analysis_type`. 

**Testing** is intended for initial runs (aka "is this working" scenarios). It is run with minimal settings (1 Gibbs sampling chain using one thread, running 100 burn-in iterations, collecting 10 posterior samples off of each chain with 10 iterations between each, to allow for a short execution time. As this is for testing purposes, these run settings cannot be changed. 

**Analysis** is used for full-analysis runs, utilising 20 Gibbs sampling chains across 20 threads. These chains run 30,000 burn-in iterations and collect 100 posterior samples from each chain, with 200 iterations collected between each sample. These are standardised settings, optimised and conducted routinely. Users can change these settings by specifying `--burnin_iterations`, `--burnin_multiplier`, `--posterior`, and `--posterior_iterations`.

The HDP pipeline can be run with **hierarchy**, supporting single and double tier hierarchy (setting `hierarchy = double` or `hierarchy = single`) or **no hierarchy** (setting `hierarchy=flat`). In hierarchy runs, the hierarchy table is specified via `--hierarchy_matrix /path/to/hierarchy_matrix`. Hierarchy parameter(s), which define how hierarchy is constructed, must be their own column(s) in the hierarchy matrix and may specified via their respective column name(s), specified to the pipeline by `--hierarchy_parameter1 <column_name_of_hierarchy_parameter1>` and, if double hierarchy, `--hierarchy_parameter2 <column_name_of_hierarchy_parameter2>`. The default for the pipeline is **no hierarchy**. 

HDP can be run with **prior** (setting `prior = true`, specifying the prior matrix via `--prior_matrix /path/to/prior_matrix`) or **no prior** (setting `prior=false`). The pipeline defaults to **no priors**. Pseudocounts can be assigned to prior signatures provided via `--prior_pseudocounts <pseudocounts1,pseudocounts2,pseudocounts3>`. Further guidance is provided below. 

The primary output of HDP, found under the subdirectory `/HDP_ExtractedSigs/`) include the extracted *de novo* signature table (`HDP_deNovoSignatures.txt`) and standard QC plots and matrices.

### SigProfilerPlotting
SigProfilerPlotting is used to plot the extracted *de novo* signature spectra. The pipeline will feed the extracted *de novo* signature table (i.e., `mSigHdp_deNovoSignatures.txt`) and the user-specified mutational context (done via `--mutational_context SBS96`, default is SBS96) to SigProfilerPlotting for plotting. The output of SigProfilerPlotting can be found under the directory `/Signature_Spectra/`, further broken down to subdirectory `/DeNovoSignatures/`, containing PDFs of the spectra plots. 

Setting `plotting` to `false` will turn this functionality off. 

### SigProfilerAssignment
SigProfilerAssignment is used to decompose the extracted *de novo* signatures. The pipeline feeds the extracted *de novo* signature table and the previous, user-provided mutational matrix (specified via `mutational matrix`) to SigProfilerAssignment, which executes its `decompose_fit()` function. The output of SigProfilerAssignment is located under the directory `/SigProfilerDecomposition/`, containing subdirectories `/Activities/`, `/Signatures/`, and `/Solution_Stats/` and decomposition plots and mappings to COSMIC signatures (e.g., [COSMIC SBS](https://cancer.sanger.ac.uk/signatures/sbs/)).

Setting `decompose` to `false` will turn this functionality off. 

## Dependencies
* Nextflow >= 24.04.2 required
* Python, required packages: [SigProfilerPlotting](https://github.com/AlexandrovLab/SigProfilerPlotting), [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment), [argparse](https://docs.python.org/3/library/argparse.html)
* R, required packages: [HDP](https://github.com/nicolaroberts/hdp), [tidyverse](https://www.tidyverse.org/), [argparse](https://cran.r-project.org/web/packages/argparse/index.html), [Tidyverse](https://tidyverse.org/)

### Docker Images
The following Docker images have been created to facilitate the running of the HDP pipeline and are available from DockerHub. The Docker files used to create them are available under [docker_files](https://github.com/ahsr-cell/HDP_pipeline/tree/main/docker_files). 

* HDP: ar39944/hdp_amd64:0.0.2
* SigProfilerPlotting: ar39944/sigprofilerplotting_amd64:0.0.1
* SigProfilerAssignment ar39944/sigprofilerassignment_amd64:0.0.3

## Installation
Clone this repository via

 > git clone git@github.com:ahsr-cell/HDP_pipeline.git

## Usage

### Input files and parameters

| Input      | Description |
| ----------- | ----------- |
| `mutation_matrix`      | Required input file, provided as a tab-delimited file (.tsv). The expected format is a matrix with one row per mutation type and one column per sample. Include the mutation types as the first column labelled as 'MutationType'. It is highly recommended to generate mutation matrices via [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator). For an example. please see [example_mutation_matrix.tsv](https://github.com/ahsr-cell/HDP_pipeline/blob/main/docs/example_input_data/example_mutation_matrix.tsv) for an example.        |
| `hierarchy`   | Required value, provided as a string. Options are `double`, `single` or `flat`         |
| `hierarchy_matrix`   | Optional input file, required if `hierarchy = double` or `hierarchy = single`. The expected format is a matrix with one row per sample ID (matching the input mutation matrix) and one column specifying hierarchy groupings. See [example_hierarchy_matrix.tsv](https://github.com/ahsr-cell/HDP_pipeline/blob/main/docs/example_input_data/example_hierarchy_matrix.tsv) for an example.         |
| `hierarchy_parameter1` and `hierarchy_parameter2`  | Optional values, required if `hierarchy = double` or `hierarchy = single`. This should be formatted exactly as the column name specifying hierarchy groupings in the input `hierarchy_matrix`. E.g., if a user provided [example_hierarchy_matrix.tsv](https://github.com/ahsr-cell/HDP_pipeline/blob/main/docs/example_input_data/example_hierarchy_matrix.tsv), `hierarchy_parameter1` would be `hierarchy_parameter1=Grouping` and `hierarchy_parameter2=Subgrouping`. Note that the algorithm does *not* sort the hierarchy key; thus, users are advised to sort their hierarchy key so that their sample ID and primary hierarchy parameter are columns in matching, alphanumeric order.              |
| `prior`   | Required value. Options are `true` or `false`          |
| `prior_matrix`   | Optional input file, required if `prior = true`. The expected format is a matrix with one column per mutational signature (e.g., SBS1, SBS5, etc.) and one column per mutation type (e.g., A[C>A]A). Include the mutation types as the first column labelled as 'MutationType'. Note that the prior matrix variant class should match the mutation matrix variant class; in other words, if the inputted mutation matrix is SBS96, the prior matrix should be SBS96. See [example_prior_matrix.tsv](https://github.com/ahsr-cell/HDP_pipeline/blob/main/docs/example_input_data/example_prior_matrix.tsv) for an example.         |
| `prior_pseudocounts`   | Optional value, provided as a string, required if `prior = true`, default is `1000`. Note, the pipeline will detect how many prior signatures are provided through the input prior matrix. If a user wants to assign the same amount of pseudocounts to each signature, they should provide a single value (i.e., if a user has SBS1, SBS5, and SBS88 and wants to assign 500 pseudocounts to each signature, they can simply provide `prior_pseudocounts=500`). If a user wants to assign different pseudocounts values to their signatures (i.e., SBS1 has 300 pseudocounts, SBS5 has 300 pseudocounts, and SBS88 has 500 pseudocounts), a user should provide the values separated by a comma, for example: `prior_pseudocounts=300,300,500`. See [below](https://github.com/ahsr-cell/HDP_pipeline?tab=readme-ov-file#pseudocounts) for further guidance on usage of pseudocounts in HDP.          |
| `analysis_type`   | Required value, provided as a string. Options are `analysis` or `testing`, default is `analysis`         |
| `mutation_context`   | Required value, provided as a string. Options are `SBS96`, `DBS78`, `ID83`, default is `SBS96`         |
| `plotting`   | Required value, provided as a string. Options are `true` or `false`, default is `true`         |
| `decompose`   | Required value, provided as a string. Options are `true` or `false`, default is `true`         |
| `outdir`   | Required path, specifying the location of output files generated by pipeline. See [`output.md`](docs/output.md) for the files/directories that will be contained within this directory.        |
| `user_defmemory`   | Optional value that can be provided if a user would like to specify the amount of memory to request for each HDP run. Note, by default, the pipeline will calculate the amount of memory required, based on the number of samples and highest mutation burden. By providing a value to `user_defmemory`, this process will be turned off.        |
| `burnin_iterations`[^*]   | Optional value, provided as a double. Used to change number of burn-in iterations conducted. Default is `30000`          |
| `posterior`[^*]   | Optional value, provided as a double. Used to change number of posterior samples collected. Default is `100`         |
| `posterior_iterations`[^*]   | Optional value, provided as a double. Used to change number of iterations conducted between each posterior sample. Default is `1000`         |
| `concentration_parameter`[^*]   | Optional value, provided as a double. Used to change number of iterations concentration parameter sampling. Default is `3`         |
| `chains`[^*]   | Optional value, provided as a double. Used to change number of posterior sampling chains and number of CPUs to use (for parallelisation). Default is `20`         |
| `threshold`[^*]   | Optional value, provided as a double. Used to change minimum number of mutations required in a sample. Default is `20`         |
[^*]: HDP run options that can be changed in analysis runs. Recommended to change only after conducting an initial analysis run with default settings and having validated reasonable cause for changing number of iterations.

### Pseudocounts
The HDP pipeline has functionality for prior pseudocounts. Pseudocounts are synthetic counts included at the start of each HDP structure, thereby providing an initial number of mutations directly attributable to the respective prior signature.

The number of pseudocounts provided to each prior signature is the ‘strength’. However, the number of pseudocounts to assign a prior signature is influenced by numerous factors, including but not limited to overall mutation burden, number of prior signatures, context and scientific question, etc.

Thus, it is recommended that users optimise their HDP run to test the run settings most appropriate for their data and research question. However, the following table has been generated as a guideline, with the understanding that the amount of pseudocounts assigned and ‘strength’ will be relative to dataset and context.

To illustrate this, suppose you are analysing 50 samples with moderate mutation burdens.

We know that SBS1 and SBS5 (clock-like signatures) are biologically expected in all tissues [(Boysen et al. 2025)](https://doi.org/10.1093/nar/gkaf1149). However, you could not detect them in samples with very low mutation counts.

Without pseudocounts: HDP might under-detect SBS1/5 in low-mutation samples because their weak signal is overshadowed by noisier signatures. With pseudocounts: By adding pseudocounts to SBS1/5 gives them a baseline presence, encouraging HDP to retain them.

If you have time, it is recommended you run HDP with different pseudocounts values as in the table below for the detailed example.

In general, fewer pseudocounts allows the model to freely discover mutation patterns, while in contrast, more pseudocounts lets the model lean on the provided prior signatures. On one hand, too few pseudocounts may ignore weak but biologically real priors. On the other hand, adding too many creates bias by forcing priors to present even when the data does not support them. Therefore, choosing the number of pseudocounts relies on the experience and the nature of your data. If you think of a signature should be there, then use high number of pseudocounts, and vice versa.

When you first use the pipeline, you can start with adding 100 pseudocounts. However, I would say you would not see much of the difference in the extraction nor attribution. One approach could be using the moderate prior strength, let’s say 400, then comparing the results to the run without prior. You will have a better idea of either increase or decrease pseudocounts based on biological interpretations and model stability.

| Relative Strength      | Pseudocounts | % of mutation data, *assuming a dataset of 20,000 real mutations* |
| ----------- | ----------- | ----------- |
| Weak      | 100 counts * 10 signatures = 1,000 pseudcounts        | 5% *(1,000/20,000)* |
| &#8593;      | 400 counts * 10 signatures = 4,000 pseudcounts        | 20% *(4,000/20,000)* |
| &#8595;      | 800 counts * 10 signatures = 8,000 pseudcounts        | 40% *(8,000/20,000)* |
| Strong      | 1200 counts * 10 signatures = 12,000 pseudcounts        | 60% *(12,000/20,000)* |


The pipeline can be run using (removing unused options):

```
     nextflow run /path/to/HDP_pipeline/main.nf \
     -profile <docker/singularity/.../institute> \
     -c /path/to/config_file \
     --mutational_matrix /path/to/mutation_matrix.tsv \
     --hierarchy <flat/single/double> \
     --hierarchy_matrix /path/to/hierarchy_key.tsv \
     --hierarchy_parameter1 <hierarchy_parameter1_column-name_in_hierarchy_key> \
     --hierarchy_parameter2 <hierarchy_parameter2_column-name_in_hierarchy_key> \
     --prior <true/false> \
     --prior_matrix /path/to/prior_matrix.tsv \
     --prior_pseudocounts <pseudocounts_values_to_assign> \
     --mutation_context <SBS96/DBS78/ID83> \
     --analysis_type <analysis/testing> \
     --outdir /path/to/outdir/ \
     --plotting <true/false> \
     --decompose <true/false> \
     --user_defmemory <desired_memory_inGB> \
     --numchains <number_of_desired_chains> \
     --threshold <minimum_mutation_threshold> \
     --burnin_iterations <number_of_desired_burnin_chains> \
     --posterior <number_of_desired_posterior_samples> \
     --posterior_iterations <number_of_desired_posterior_iterations> \
     --concentration_parameter <value_of_desired_concentration_parameter> \
     -resume
```

### Sanger Users
Sangers can run the pipeline using the wrapper script below. Note: users will need to remove/comment out options in the Nextflow run command that not required for their run:

```
#! /usr/bin/env bash

#BSUB -J HDP_pipeline_run
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -u user@sanger.ac.uk
#BSUB -q basement #set time requirements
#BSUB -n 1 #single node needed for Nextflow job submitter
#BSUB -M5000 #set memory requirements
#BSUB -R "select[mem>5000] rusage[mem=5000]" #set memory requirements

### Run inputs and options, change accordingly
mutational_matrix=/path/to/mutation_matrix.tsv
hierarchy=flat #Options include flat, single, double 

prior=false #if running prior, set to true
 
mutation_context=SBS96 #For SigProfilerPlotting, options include SBS96, DBS78, ID83
analysis_type=analysis #'testing' for test run, 'analysis' for full analysis run

hierarchy_matrix=/path/to/input/hierarchy_key.tsv #if running no hierarchy, delete/comment this line out 
hierarchy_parameter1=parameter1 #The column name of the primary hierarchy parameter
hierarchy_parameter2=parameter2 #The column name of secondary hierarchy parameter
prior=true #if running prior, set to true
prior_matrix=/path/to/input/prior_matrix.tsv #if running no prior, delete/comment this line out
prior_pseudocounts=1000 #If desired, the amount of pseudocounts to assign to prior signatures
outdir=/path/to/output #set to location of output files
plotting=true #set to false if you do not want to plot extracted signature spectra
decompose=true #set to false if you do not want to decompose extracted signatures 
user_defmemory=0 #Memory resources to request. Set if known, if not, leave as 0. Units are in Gigabytes.

numchains=20 #Number of posterior sampling chains to run
threshold=0 #Minimum number of mutations per sample 
burnin_iterations=30000 #Number of burn-in iterations
posterior=100 #Number of posterior samples
posterior_iterations=200 #Number of iterations conducted between posterior samples
concentration_parameter=3 #

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
     --prior ${prior} \
     --prior_matrix ${prior_matrix} \
     --prior_pseudocounts ${prior_pseudocounts} \
     --mutation_context ${mutation_context} \
     --analysis_type ${analysis_type} \
     --outdir ${outdir} \
     --plotting ${plotting} \
     --decompose ${decompose} \
     --user_defmemory ${user_defmemory} \
     --numchains ${numchains} \
     --threshold ${threshold} \
     --burnin_iterations ${burnin_iterations} \
     --posterior ${posterior} \
     --posterior_iterations ${posterior_iterations} \
     -resume
```
## Pipeline output

For more details about the output files and reports, please refer to the [`output.md`](docs/output.md) output documentation.

## Credits

The HDP pipeline was written by Andrew Ramos and Phuong Le. 

We thank the following people and teams for their assistance in the development of this pipeline:

Mimy Pham, Yichen Wang, Sarah Moody, and CASM IT

## Contributions and Support

Please feel free to contribute by either creating a pull request or create a new issue on this GitHub repo.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.