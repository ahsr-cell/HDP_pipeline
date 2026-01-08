process Prior_validation {

    input:
    path prior_matrix
    val prior_pseudocounts

    output:
    path "validation/prior_sig_pseudocounts.csv"
    
    script:
    """
    Rscript --vanilla ${projectDir}/bin/Prior_validation.R --prior_matrix ${prior_matrix} --prior_pseudocounts ${prior_pseudocounts}
    """
}