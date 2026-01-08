#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { MutMatrix_resourcereqs } from './modules/MutMatrix_resourcereqs.nf'
include { MutMatrix_resourcereqs_singlehierarchy } from './modules/MutMatrix_resourcereqs_singlehierarchy.nf'
include { MutMatrix_resourcereqs_doublehierarchy } from './modules/MutMatrix_resourcereqs_doublehierarchy.nf'
include { Prior_validation } from './modules/Prior_validation.nf'
include { HDP_flat_prior } from './modules/HDP_flat_prior.nf'
include { HDP_flat_noprior } from './modules/HDP_flat_noprior.nf'
include { HDP_single_prior } from './modules/HDP_single_prior.nf'
include { HDP_single_noprior  } from './modules/HDP_single_noprior.nf'
include { HDP_double_prior } from './modules/HDP_double_prior.nf'
include { HDP_double_noprior  } from './modules/HDP_double_noprior.nf'
include { HDP_combine } from './modules/HDP_combine.nf'
include { SigProfilerPlotting as SigPlt_Extracted } from './modules/SigProfilerPlotting.nf'
include { SigProfilerAssignment as SigPA_Extracted } from './modules/SigProfilerAssignment.nf'

//
// WORKFLOW: Run main analysis pipeline depending on user inputs
//

workflow {

    main:
    //
    // WORKFLOW: Full suite of analysis: HDP, SigProfilerPlotting, and SigProfilerAssignment
    //
    
    if (params.hierarchy == "double") {
        MutMatrix_resourcereqs_doublehierarchy(
            params.mutational_matrix,
            params.hierarchy_matrix,
            params.hierarchy_parameter1,
            params.hierarchy_parameter2,
            params.user_defmemory
        )

        memory_requirements_ch = MutMatrix_resourcereqs_doublehierarchy.out.memory_reqs_matrix
                        .splitCsv( header: true )
                        .map { row -> tuple( row.Sample_number, row.Mutation_burden, row.Memory_required )
                        }
    }
    if (params.hierarchy == "single") {
        MutMatrix_resourcereqs_singlehierarchy(
            params.mutational_matrix,
            params.hierarchy_matrix,
            params.hierarchy_parameter1,
            params.user_defmemory
        )
        memory_requirements_ch = MutMatrix_resourcereqs_singlehierarchy.out.memory_reqs_matrix
                        .splitCsv( header: true )
                        .map { row -> tuple( row.Sample_number, row.Mutation_burden, row.Memory_required )
                        }
    }
    if (params.hierarchy == "flat") {
        MutMatrix_resourcereqs(
        params.mutational_matrix,
        params.user_defmemory
        )
        
        memory_requirements_ch = MutMatrix_resourcereqs.out.memory_reqs_matrix
                        .splitCsv( header: true )
                        .map { row -> tuple( row.Sample_number, row.Mutation_burden, row.Memory_required )
                        }
    }

    if (params.prior == true) {
        Prior_validation(
        params.prior_matrix,
        params.prior_pseudocounts
        )
    }

    if (params.hierarchy == "double") {
        if (params.prior == true) {
            HDP_double_prior(
             Channel.of(1..params.numchains)
             .combine(memory_requirements_ch),
             params.mutational_matrix,
             params.hierarchy_matrix,
             params.prior_matrix,
             params.prior_pseudocounts,
             params.analysis_type, 
             params.hierarchy_parameter1,
             params.hierarchy_parameter2,
             params.burnin_iterations,
             params.posterior,
             params.posterior_space,
             params.threshold
             )
            // HDP_double_prior(
            //  params.mutational_matrix,
            //  params.hierarchy_matrix,
            //  params.prior_matrix,
            //  params.prior_pseudocount,
            //  params.analysis_type, 
            //  params.hierarchy_parameter1,
            //  params.hierarchy_parameter2,
            //  params.burnin_iterations,
            //  params.posterior,
            //  params.posterior_space,
            //  params.threshold,
            //  Channel.of(1..params.numchains),
            //  memory_requirements_ch
            //  )

             HDP_collected = HDP_double_prior.out.collect().map { file("${params.outdir}/HDP_chains") }
            
             HDP_combine(
                 params.mutational_matrix,
                 params.hierarchy_matrix,
                 HDP_collected,
                 params.numchains,
                 params.threshold,
                 params.mutation_context
             )
             
             if (params.plotting == true) {
                 if (params.decompose == true) {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
                 SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             } else {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
             }             
         } else {
             if (params.decompose == true) {
             SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             }
         }
        } else {
            HDP_double_noprior(
             Channel.of(1..params.numchains)
             .combine(memory_requirements_ch),
             params.mutational_matrix,
             params.hierarchy_matrix,
             params.analysis_type, 
             params.hierarchy_parameter1,
             params.hierarchy_parameter2,
             params.burnin_iterations,
             params.posterior,
             params.posterior_space,
             params.threshold
             )
            // HDP_double_noprior(
            //  params.mutational_matrix,
            //  params.hierarchy_matrix,
            //  params.analysis_type, 
            //  params.hierarchy_parameter1,
            //  params.hierarchy_parameter2,
            //  params.burnin_iterations,
            //  params.posterior,
            //  params.posterior_space,
            //  params.threshold,
            //  Channel.of(1..params.numchains),
            //  memory_requirements_ch
            //  )
    
             HDP_collected = HDP_double_noprior.out.collect().map { file("${params.outdir}/HDP_chains") }
            
             HDP_combine(
                 params.mutational_matrix,
                 params.hierarchy_matrix,
                 HDP_collected,
                 params.numchains,
                 params.threshold,
                 params.mutation_context
             )

             if (params.plotting == true) {
                 if (params.decompose == true) {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
                 SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             } else {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
             }             
         } else {
             if (params.decompose == true) {
             SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             }
         }
        }
     }
     if (params.hierarchy == "single") {
        if (params.prior == true) {
            HDP_single_prior(
             Channel.of(1..params.numchains)
             .combine(memory_requirements_ch),
             params.mutational_matrix,
             params.hierarchy_matrix,
             params.prior_matrix,
             params.prior_pseudocounts,             
             params.analysis_type,
             params.hierarchy_parameter1,
             params.burnin_iterations,
             params.posterior,
             params.posterior_space,
             params.threshold
             )
            // HDP_single_prior(
            //  params.mutational_matrix,
            //  params.hierarchy_matrix,
            //  params.prior_matrix,
            //  params.prior_pseudocount,             
            //  params.analysis_type,
            //  params.hierarchy_parameter1,
            //  params.burnin_iterations,
            //  params.posterior,
            //  params.posterior_space,
            //  params.threshold,
            //  Channel.of(1..params.numchains),
            //  memory_requirements_ch
            //  )
             HDP_collected = HDP_single_prior.out.collect().map { file("${params.outdir}/HDP_chains") }
            
             HDP_combine(
                 params.mutational_matrix,
                 params.hierarchy_matrix,
                 HDP_collected,
                 params.numchains,
                 params.threshold,
                 params.mutation_context
             )

             if (params.plotting == true) {
                 if (params.decompose == true) {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
                 SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             } else {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
             }             
         } else {
             if (params.decompose == true) {
             SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             }
         }
        } else {
            HDP_single_noprior(
             Channel.of(1..params.numchains)
             .combine(memory_requirements_ch),
             params.mutational_matrix,
             params.hierarchy_matrix,
             params.analysis_type,
             params.hierarchy_parameter1,
             params.burnin_iterations,
             params.posterior,
             params.posterior_space,
             params.threshold
            )
            
            // HDP_single_noprior(
            //  params.mutational_matrix,
            //  params.hierarchy_matrix,
            //  params.analysis_type,
            //  params.hierarchy_parameter1,
            //  params.burnin_iterations,
            //  params.posterior,
            //  params.posterior_space,
            //  params.threshold,
            //  Channel.of(1..params.numchains),
            //  memory_requirements_ch
            //  )

             HDP_collected = HDP_single_noprior.out.collect().map { file("${params.outdir}/HDP_chains") }
            
             HDP_combine(
                 params.mutational_matrix,
                 params.hierarchy_matrix,
                 HDP_collected,
                 params.numchains,
                 params.thresold,
                 params.mutation_context
             )

             if (params.plotting == true) {
                 if (params.decompose == true) {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
                 SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             } else {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
             }             
         } else {
             if (params.decompose == true) {
             SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             }
         }
        }
     }
     if (params.hierarchy == "flat") {
        if (params.prior == true) {
            HDP_flat_prior(
             Channel.of(1..params.numchains)
             .combine(memory_requirements_ch),
             params.mutational_matrix,
             params.prior_matrix,
             params.prior_pseudocounts,
             params.analysis_type, 
             params.burnin_iterations,
             params.posterior,
             params.posterior_space,
             params.threshold
             )

             
            // HDP_flat_prior(
            //  params.mutational_matrix,
            //  params.prior_matrix,
            //  params.prior_pseudocount,
            //  params.analysis_type, 
            //  params.burnin_iterations,
            //  params.posterior,
            //  params.posterior_space,
            //  params.threshold,
            //  Channel.of(1..params.numchains),
            //  memory_requirements_ch
            //  )
             HDP_collected = HDP_flat_prior.out.collect().map { file("${params.outdir}/HDP_chains") }
            
             HDP_combine(
                 params.mutational_matrix,
                 // params.hierarchy_matrix,
                 HDP_collected,
                 params.numchains,
                 params.threshold,
                 params.mutation_context
             )

             if (params.plotting == true) {
                 if (params.decompose == true) {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
                 SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             } else {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
             }             
         } else {
             if (params.decompose == true) {
             SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             }
         }
        } else {
             
            //  Channel.of(1..params.numchains)
            //  .combine(memory_requirements_ch)
            //  .view()

             HDP_flat_noprior(
             Channel.of(1..params.numchains)
             .combine(memory_requirements_ch),
             params.mutational_matrix,
             params.analysis_type, 
             params.burnin_iterations,
             params.posterior,
             params.posterior_space,
             params.threshold,
             params.mutation_context
             )

            //  HDP_flat_noprior(
            //  memory_requirements_ch,
            //  params.mutational_matrix,
            //  params.analysis_type, 
            //  params.burnin_iterations,
            //  params.posterior,
            //  params.posterior_space,
            //  params.threshold,
            //  Channel.of(1..params.numchains)
            //  )
             
            HDP_collected = HDP_flat_noprior.out.collect().map { file("${params.outdir}/HDP_chains") }
            
             HDP_combine(
                 params.mutational_matrix,
                 HDP_collected,
                 params.numchains,
                 params.threshold,
                 params.mutation_context
            )
            if (params.plotting == true) {
                 if (params.decompose == true) {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
                 SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             } else {
                 SigPlt_Extracted(
                    params.mutation_context,
                    HDP_combine.out.deNovo_extractedsigs
                 )
             }             
         } else {
             if (params.decompose == true) {
             SigPA_Extracted(
                 HDP_combine.out.deNovo_extsigs_sigPA,
                 params.mutational_matrix
             )
             }
         }
        }
     }
}
        //      if (params.plotting == true) {
        //          if (params.decompose == true) {
        //              SigPlt_Extracted(
        //                 params.mutation_context,
        //                 HDP_combine.out.deNovo_extractedsigs
        //              )
        //          SigPA_Extracted(
        //              HDP_combine.out.deNovo_extsigs_sigPA,
        //              params.mutational_matrix,
        //              params.outdir
        //          )
        //          } else {
        //          SigPlt_Extracted(
        //             params.mutation_context, 
        //             HDP_combine.out.deNovo_extractedsigs
        //          )
        //          }
        //  } else {
        //      if (params.decompose == true) {
        //      SigPA_Extracted(
        //          HDP_combine.out.deNovo_extsigs_sigPA,
        //          params.mutational_matrix
        //      )
        //      }
        //  }

 /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/