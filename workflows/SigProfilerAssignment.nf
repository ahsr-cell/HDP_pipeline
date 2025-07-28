process SigProfilerAssignment {

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input: 
    path deNovoSignatures_matrix
    path mutational_matrix
    path output_directory

    output:
    path "SigProfilerDecomposition"

    script:
    """
    SigProfilerAssignment_wrapper.py --mutational_matrix ${mutational_matrix} --deNovoSignatures_matrix ${deNovoSignatures_matrix} --output_directory ${output_directory}
    """
}