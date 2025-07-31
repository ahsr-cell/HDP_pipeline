process SigProfilerPlotting {

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input: 
    path deNovoSignatures_matrix
    val mutational_context

    output:
    path "Signature_Spectra"

    script:
    """
    SigProfilerPlotting.py --mutational_context ${mutational_context} --deNovoSignatures_matrix ${deNovoSignatures_matrix} --output_directory Signature_Spectra/
    """
}