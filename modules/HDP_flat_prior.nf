process HDP_flat_prior {
    publishDir "${params.outdir}/HDP_chains", mode: "copy"

    input:
    tuple val(chain_index), val(Sample_number), val(Mutation_burden), val(Memory_required)
    path mutational_matrix
    path prior_matrix
    val prior_pseudocount    
    val analysis_type 
    val burnin_iterations 
    val posterior 
    val posterior_space
    val threshold
    val chain_index

    output:
    // path "HDP_chains", emit: HDP_chains
    path "hdp_chain_${chain_index}.Rdata"

    script:
    """
    Rscript --vanilla ${projectDir}/bin/HDP_flat_prior.R --prior_matrix ${prior_matrix} --prior_pseudocount ${prior_pseudocount} -a ${analysis_type} -b ${burnin_iterations} -o ${posterior} -i ${posterior_space} -t ${threshold} -n ${chain_index} ${mutational_matrix}
    """
}