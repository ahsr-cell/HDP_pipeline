process HDP_flat_noprior {
    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    path mutational_matrix
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
    Rscript --vanilla ${projectDir}/bin/HDP_flat_noprior.R -a ${analysis_type} -b ${burnin_iterations} -o ${posterior} -i ${posterior_space} -t ${threshold} -n ${chain_index} ${mutational_matrix}
    """
}