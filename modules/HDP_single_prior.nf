process HDP_single_prior {
    publishDir "${params.outdir}/HDP_chains", mode: "copy"

    input:
    path mutational_matrix
    path hierarchy_matrix
    path prior_matrix
    val analysis_type 
    val hierarchy_parameter1
    val burnin_iterations 
    val posterior 
    val posterior_space
    val threshold
    val chain_index
    tuple val(Sample_number), val(Mutation_burden), val(Memory_required)

    output:
    // path "HDP_chains", emit: HDP_chains
    path "hdp_chain_*.Rdata"
    
    script:
    """
    Rscript --vanilla ${projectDir}/bin/HDP_single_prior.R --hierarchy_matrix ${hierarchy_matrix} --hierarchy_parameter1 ${hierarchy_parameter1} --prior_matrix ${prior_matrix} --analysis_type ${analysis_type} --burnin_iterations ${burnin_iterations} --posterior ${posterior} --posterior_iterations ${posterior_space} --threshold ${threshold} --chain_index ${chain_index} ${mutational_matrix}
    """
}