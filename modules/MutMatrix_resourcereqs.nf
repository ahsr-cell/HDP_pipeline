process MutMatrix_resourcereqs {

    input:
    path mutational_matrix
    val user_defmemory

    output:
    path "memory_requirements/memory_requirements.csv", emit: memory_reqs_matrix

    script:
    """
    Rscript --vanilla ${projectDir}/bin/MutMatrix_resourcereqs.R -umem ${user_defmemory} ${mutational_matrix}
    """
}