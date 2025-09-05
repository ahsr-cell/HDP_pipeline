process HDP_flat_prior {
    publishDir "${params.outdir}/HDP_chains", mode: "copy"

    input:
    path mutational_matrix
    path prior_matrix
    val analysis_type 
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
    Rscript --vanilla ${projectDir}/bin/HDP_flat_prior.R --prior_matrix ${prior_matrix} -a ${analysis_type} -b ${burnin_iterations} -o ${posterior} -i ${posterior_space} -t ${threshold} -n ${chain_index} ${mutational_matrix}
    """
}