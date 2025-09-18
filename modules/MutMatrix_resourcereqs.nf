process MutMatrix_resourcereqs {

    input:
    path mutational_matrix

    output:
    path "memory_requirements/memory_requirements.csv", emit: memory_reqs_matrix

    script:
    """
    Rscript --vanilla ${projectDir}/bin/MutMatrix_resourcereqs.R ${mutational_matrix}
    """
}