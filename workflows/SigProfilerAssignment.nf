process SigProfilerAssignment {

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input: 
    path HDP_chains
    path mutational_matrix

    output:
    path "SigProfilerDecomposition"

    script:
    """
    SigProfilerAssignment.py --mutational_matrix ${mutational_matrix} --HDP_chains ${HDP_chains}
    """
}