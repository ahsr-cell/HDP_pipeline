#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { HDP_flat } from './workflows/HDP_flat.nf'
include { HDP_single } from './workflows/HDP_single.nf'
include { HDP_double } from './workflows/HDP_double.nf'
include { HDP_combine } from './workflows/HDP_combine.nf'
include { SigProfilerPlotting as SigPlt_Extracted } from './workflows/SigProfilerPlotting.nf'
include { SigProfilerAssignment as SigPA_Extracted } from './workflows/SigProfilerAssignment.nf'

//
// WORKFLOW: Run main analysis pipeline depending on user inputs
//

workflow {

    main:
    //
    // WORKFLOW: Full suite of analysis: HDP, SigProfilerPlotting, and SigProfilerAssignment
    //
     if (params.hierarchy == "double") {
        prior_matrix_file = file(params.filter, checkIfExists:true)
         HDP_double(
             params.mutational_matrix,
             params.hierarchy_matrix,
             prior_matrix_file,
             params.analysis_type, 
             params.hierarchy_parameter1,
             params.hierarchy_parameter2,
             params.burnin_iterations,
             params.posterior,
             params.posterior_space,
             params.threshold,
             Channel.of(1..20)
             )
             HDP_combine(
              params.mutational_matrix,
              params.hierarchy_matrix,
              HDP_double.out.HDP_chains,
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
     if (params.hierarchy == "single") {
        prior_matrix_file = file(params.filter, checkIfExists:true)
         HDP_single(
             params.mutational_matrix,
             params.hierarchy_matrix,
             prior_matrix_file,
             //params.prior_matrix,             
             params.analysis_type,
             params.hierarchy_parameter1,
             params.burnin_iterations,
             params.posterior,
             params.posterior_space,
             params.threshold,
             Channel.of(1..20)
             )
             HDP_combine(
                 params.mutational_matrix,
                 params.hierarchy_matrix,
                 HDP_single.out.HDP_chains,
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
     if (params.hierarchy == "flat") {
        prior_matrix_file = file(params.filter, checkIfExists:true)
         HDP_flat(
             params.mutational_matrix,
             params.hierarchy_matrix,
             prior_matrix_file,
            // params.prior_matrix,
             params.analysis_type, 
             params.burnin_iterations,
             params.posterior,
             params.posterior_space,
             params.threshold,
             Channel.of(1..20)
             )
             HDP_combine(
                 params.mutational_matrix,
                 params.hierarchy_matrix,
                 HDP_flat.out.HDP_chains,
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
                     params.mutational_matrix,
                     params.outdir
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

 /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/