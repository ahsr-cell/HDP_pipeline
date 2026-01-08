process HDP_double_noprior {
    publishDir "${params.outdir}/HDP_chains", mode: "copy"

    input:
    tuple val(chain_index), val(Sample_number), val(Mutation_burden), val(Memory_required)
    path mutational_matrix
    path hierarchy_matrix
    val analysis_type 
    val hierarchy_parameter1
    val hierarchy_parameter2
    val burnin_iterations 
    val posterior 
    val posterior_space
    val threshold
    
    output:
    // path "HDP_chains", emit: HDP_chains
    path "hdp_chain_${chain_index}.Rdata"

    script:
    """
    Rscript --vanilla ${projectDir}/bin/HDP_double_noprior.R --hierarchy_matrix ${hierarchy_matrix} -hp1 ${hierarchy_parameter1} -hp2 ${hierarchy_parameter1} -a ${analysis_type} -b ${burnin_iterations} -o ${posterior} -i ${posterior_space} -n ${chain_index} -t ${threshold} ${mutational_matrix}
    """
}