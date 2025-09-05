process HDP_double_prior {
    publishDir "${params.outdir}/HDP_chains", mode: "copy"

    input:
    path mutational_matrix
    path hierarchy_matrix
    path prior_matrix
    val analysis_type
    val hierarchy_parameter1
    val hierarchy_parameter2
    val burnin_iterations 
    val posterior
    val posterior_space
    val threshold
    val chain_index

    output:
    // path "HDP_chains", emit: HDP_chains
    path "hdp_chain_*.Rdata"

    script:
    """
    Rscript --vanilla ${projectDir}/bin/HDP_double_prior.R -hierarchy ${hierarchy_matrix} -hp1 ${hierarchy_parameter1} -hp2 ${hierarchy_parameter1} --prior_matrix ${prior_matrix} -a ${analysis_type} -b ${burnin_iterations} -o ${posterior} -i ${posterior_space} -n ${chain_index} -t ${threshold} ${mutational_matrix}
    """
}