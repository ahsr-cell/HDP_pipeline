process SigProfilerPlotting {

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input: 
    val mutation_context
    path deNovoSignatures_matrix

    output:
    path "Signature_Spectra"

    script:
    """
    SigProfilerPlotting.py --mutation_context ${mutation_context} --deNovoSignatures_matrix ${deNovoSignatures_matrix} --output_directory Signature_Spectra/
    """
}