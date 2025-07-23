process HDP_double {
    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    path mutational_matrix
    path hierarchy_matrix
    path prior_matrix
    val analysis_type 
    val burnin_iterations 
    val posterior 
    val posterior_space
    val n_iter
    val threshold

    output:
    path "HDP_chains", emit: HDP_chains

    script:
    """
    Rscript --vanilla ${projectDir}/bin/HDP_double.R -h ${hierarchy_matrix} -p ${prior_matrix} -a ${analysis_type} -b ${burnin_iterations} -o ${posterior} -i ${posterior_space} -n ${n_iter} -t ${threshold} ${mutational_matrix}
    """
}