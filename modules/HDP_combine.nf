process HDP_combine {
    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    path mutational_matrix
    path HDP_all_chains
    val numchains
    val threshold
    val mutation_context
    
    output:
    path "HDP_ExtractedSigs", emit: deNovo_signaturesdir
    path "HDP_ExtractedSigs/HDP_deNovoSignatures.txt", emit: deNovo_extractedsigs
    path "HDP_ExtractedSigs/HDP_deNovoSigs_sigPADecomp.txt", emit: deNovo_extsigs_sigPA

    script:
    """
    Rscript --vanilla ${projectDir}/bin/HDP_combine.R ${mutational_matrix} --HDP_chains ${HDP_all_chains} --number_chains ${numchains} -t ${threshold} -c ${mutation_context}
    """
}