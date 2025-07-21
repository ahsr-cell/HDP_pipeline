process HDP_flat {
    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    path mutational_matrix
    path hierarchy_matrix
    path prior_matrix
    val hierarchy
    val prior
    val analysis_type 
    val burning_iterations 
    val posterior 
    val posterior_space 

    output:
    path "deNovo_signatures", emit: deNovo_signaturesdir
    path "deNovo_signatures"

    script:
    """
    Rscript --vanilla ${projectDir}/bin/HDP_flat.R
    """
}