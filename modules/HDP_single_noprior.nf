process HDP_single_noprior {
    publishDir "${params.outdir}/HDP_chains", mode: "copy"

    input:
    tuple val(chain_index), val(Sample_number), val(Mutation_burden), val(Memory_required)
    path mutational_matrix
    path hierarchy_matrix
    val analysis_type 
    val hierarchy_parameter1
    val burnin_iterations 
    val posterior 
    val posterior_space
    val threshold

    output:
    // path "HDP_chains", emit: HDP_chains
    path "hdp_chain_${chain_index}.Rdata"

    script:
    """
    Rscript --vanilla ${projectDir}/bin/HDP_single_noprior.R --hierarchy_matrix ${hierarchy_matrix} --hierarchy_parameter1 ${hierarchy_parameter1} --analysis_type ${analysis_type} --burnin_iterations ${burnin_iterations} --posterior ${posterior} --posterior_iterations ${posterior_space} --threshold ${threshold} --chain_index ${chain_index} ${mutational_matrix}
    """
}