process HDP_combine {
    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    path HDP_chains
    path mutational_matrix
    
    output:
    path "deNovo_extractedsigs", emit: deNovo_signaturesdir
    path "deNovo_signatures/HDP_deNovoSignatures.txt", emit: deNovo_extractedsigs

    script:
    """
    Rscript --vanilla ${projectDir}/bin/HDP_combine.R ${mutational_matrix}
    """
}