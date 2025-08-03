process HDP_flat {
    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    path mutational_matrix
    path hierarchy_matrix
    path prior_matrix
    val analysis_type 
    val burnin_iterations 
    val posterior 
    val posterior_space
    val threshold
    val chain_index

    output:
    path "HDP_chains", emit: HDP_chains

    script:
    """
    Rscript --vanilla ${projectDir}/bin/HDP_flat.R -prior ${prior_matrix} -a ${analysis_type} -b ${burnin_iterations} -o ${posterior} -i ${posterior_space} -n ${chain_index} -t ${threshold} ${mutational_matrix}
    """
}