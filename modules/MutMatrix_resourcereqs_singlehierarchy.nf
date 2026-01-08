process MutMatrix_resourcereqs_singlehierarchy {

    input:
    path mutational_matrix
    path hierarchy_matrix
    val hierarchy_type
    val hierarchy_parameter1
    val user_defmemory

    output:
    path "memory_requirements/memory_requirements.csv", emit: memory_reqs_matrix

    script:
    """
    Rscript --vanilla ${projectDir}/bin/MutMatrix_resourcereqs_singlehierarchy.R --hierarchy_matrix ${hierarchy_matrix} --hierarchy_type ${hierarchy_type} --hierarchy_parameter1 ${hierarchy_parameter1} --user_defmemory ${user_defmemory} ${mutational_matrix}
    """
}