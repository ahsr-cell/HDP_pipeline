process MutMatrix_resourcereqs_doublehierarchy {

    input:
    path mutational_matrix
    path hierarchy_matrix
    val hierarchy_type
    val hierarchy_parameter1
    val hierarchy_parameter2
    val user_defmemory

    output:
    path "memory_requirements/memory_requirements.csv", emit: memory_reqs_matrix

    script:
    """
    Rscript --vanilla ${projectDir}/bin/MutMatrix_resourcereqs_doublehierarchy.R --hierarchy_matrix ${hierarchy_matrix} --hierarchy_type ${hierarchy_type} --hierarchy_parameter1 ${hierarchy_parameter1} --hierarchy_parameter2 ${hierarchy_parameter2} --user_defmemory ${user_defmemory} ${mutational_matrix}
    """
}