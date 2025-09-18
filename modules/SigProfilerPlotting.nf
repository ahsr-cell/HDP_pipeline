process SigProfilerPlotting {

    publishDir "${params.outdir}", mode: "copy"

    input: 
    val mutation_context
    path deNovoSignatures_matrix

    output:
    path "Signature_Spectra"

    script:
    """
    export MPLCONFIGDIR=${workDir}
    rm -rf Signature_Spectra
    mkdir Signature_Spectra

    SigProfilerPlotting.py --mutation_context ${mutation_context} --deNovoSignatures_matrix ${deNovoSignatures_matrix} --output_directory /Signature_Spectra/
    """
}