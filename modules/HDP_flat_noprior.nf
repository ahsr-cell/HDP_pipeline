process HDP_flat_noprior {
    publishDir "${params.outdir}/HDP_chains", mode: "copy"

    input:
    tuple val(chain_index), val(Sample_number), val(Mutation_burden), val(Memory_required)
    path mutational_matrix
    val analysis_type 
    val burnin_iterations 
    val posterior 
    val posterior_space
    val threshold
    val mutation_context

    output:
    // path "HDP_chains", emit: HDP_chains
    // path "**/hdp_chain_*.Rdata"
    path "hdp_chain_${chain_index}.Rdata"
    
    script:
    println(chain_index)
    """
    Rscript --vanilla ${projectDir}/bin/HDP_flat_noprior.R -a ${analysis_type} -b ${burnin_iterations} -o ${posterior} -i ${posterior_space} -t ${threshold} -n ${chain_index} -c ${mutation_context} ${mutational_matrix}
    """
}