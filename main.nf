#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { HDP_flat } from './workflows/HDP_flat.nf'
include { HDP_single } from './workflows/HDP_single.nf'
include { HDP_double } from './workflows/HDP_double.nf'
include { HDP_combine } from './workflows/HDP_combine.nf'
include { SigProfilerPlotting as SigPlt_Extracted } from './workflows/SigProfilerPlotting.nf'
include { SigProfilerPlotting as SigPlt_LowConfidence } from './workflows/SigProfilerPlotting.nf'
include { SigProfilerAssignment as SigPA_Extracted } from './workflows/SigProfilerAssignment.nf'
include { SigProfilerAssignment as SigPA_LowConfidence } from './workflows/SigProfilerAssignment.nf'

//
// WORKFLOW: Run main analysis pipeline depending on user inputs
//

workflow {

    main:

    //
    // WORKFLOW: Full suite of analysis: HDP, SigProfilerPlotting, and SigProfilerAssignment
    //
    if (params.hierarchy == "double") {
        HDP_double(
            params.mutational_matrix,
            params.hierarchy_matrix,
            params.hierarchy_parameter1,
            params.hierarchy_parameter2,
            params.prior_matrix,
            params.analysis_type, 
            params.burnin_iterations,
            params.posterior,
            params.posterior_space,
            params.n_iter,
            params.threshold
            )
            HDP_combine(
             HDP_double.out.HDP_chains,
             params.mutational_matrix,
             params.hierarchy_matrix,
             params.threshold
            )
            if (params.plotting == true) {
                if (params.decompose == true) {
                SigPlt_Extracted(
                HDP_combine.out.deNovo_extractedsigs,
                params.mutational_context
            )
                SigPA_Extracted(
                HDP_combine.out.deNovo_extractedsigs,
                params.mutational_matrix
            )
//            if (HDP_combine.out.deNovo_lowconfsigs != null){
//                SigPlt_LowConfidence(
//                    HDP_double.out.deNovo_lowconfsigs,
//                    params.mutational_context
//                )
//                SigPA_LowConfidence(
//                    HDP_combine.out.deNovo_lowconfsigs,
//                    params.mutational_matrix
//                )
//                }
            } else {
                SigPlt_Extracted(
                    HDP_combine.out.deNovo_extractedsigs,
                    params.mutational_context
            )
//            if (HDP_combine.out.deNovo_lowconfsigs != null){
//                SigPlt_LowConfidence(
//                    HDP_combine.out.deNovo_lowconfsigs,
//                    params.mutational_context
//                )
//                }
            }             
        } else {
            if (params.decompose == true) {
            HDP_double(
            params.mutational_matrix,
            params.hierarchy_matrix,
            params.hierarchy_parameter1,
            params.hierarchy_parameter2,
            params.prior_matrix,
            params.analysis_type, 
            params.burnin_iterations,
            params.posterior,
            params.posterior_space,
            params.n_iter,
            params.threshold
            )
            HDP_combine(
             HDP_double.out.HDP_chains,
             params.mutational_matrix,
             params.hierarchy_matrix,
             params.threshold   
            )
            SigPA_Extracted(
                HDP_combine.out.deNovo_extractedsigs,
                params.mutational_matrix
            )
//            if (HDP_combine.out.deNovo_lowconfsigs != null) {
//                SigPA_LowConfidence(
//                    HDP_combine.out.deNovo_lowconfsigs,
//                    params.mutational_matrix
//                )
//            }
        }
        } 
    }
    if (params.hierarchy == "single") {
        HDP_single(
            params.mutational_matrix,
            params.hierarchy_matrix,
            params.hierarchy_parameter1,
            params.prior_matrix,
            params.analysis_type, 
            params.burnin_iterations,
            params.posterior,
            params.posterior_space,
            params.n_iter,
            params.threshold
            )
            HDP_combine(
                HDP_single.out.HDP_chains,
                params.mutational_matrix,
                params.hierarchy_matrix,
                params.threshold   
            )
            if (params.plotting == true) {
                if (params.decompose == true) {
                    SigPlt_Extracted(
                    HDP_combine.out.deNovo_extractedsigs,
                    params.mutational_context
                )
                SigPA_Extracted(
                    HDP_combine.out.deNovo_extractedsigs,
                    params.mutational_matrix
                )
//               if (HDP_combine.out.deNovo_lowconfsigs != null){
//                    SigPlt_LowConfidence(
//                        HDP_single.out.deNovo_lowconfsigs,
//                        params.mutational_context
//                    )
//                    SigPA_LowConfidence(
//                        HDP_combine.out.deNovo_lowconfsigs,
//                        params.mutational_matrix
//                    )
//                }
                } else {
                SigPlt_Extracted(
                    HDP_combine.out.deNovo_extractedsigs,
                    params.mutational_context
            )
//            if (HDP_combine.out.deNovo_lowconfsigs != null){
//                SigPlt_LowConfidence(
//                    HDP_combine.out.deNovo_lowconfsigs,
//                    params.mutational_context
//                )
//                }
            }             
        } else {
            if (params.decompose == true) {
            HDP_single(
            params.mutational_matrix,
            params.hierarchy_matrix,
            params.hierarchy_parameter1,
            params.prior_matrix,
            params.analysis_type, 
            params.burnin_iterations,
            params.posterior,
            params.posterior_space,
            params.n_iter,
            params.threshold
            )
            HDP_combine(
                HDP_single.out.HDP_chains,
                params.mutational_matrix,
                params.hierarchy_matrix,
                params.threshold   
            )
            SigPA_Extracted(
                HDP_combine.out.deNovo_extractedsigs,
                params.mutational_matrix
            )
//            if (HDP_combine.out.deNovo_lowconfsigs != null) {
//                SigPA_LowConfidence(
//                    HDP_combine.out.deNovo_lowconfsigs,
//                    params.mutational_matrix
//                )
//            }
        }
        }
    }
    if (params.hierarchy == "flat") {
        HDP_flat(
            params.mutational_matrix,
            params.hierarchy_matrix,
            params.prior_matrix,
            params.analysis_type, 
            params.burnin_iterations,
            params.posterior,
            params.posterior_space,
            params.n_iter,
            params.threshold
            )
            HDP_combine(
                HDP_flat.out.HDP_chains,
                params.mutational_matrix,
                params.hierarchy_matrix,
                params.threshold   
            )
            if (params.plotting == true) {
                if (params.decompose == true) {
                    SigPlt_Extracted(
                    HDP_combine.out.deNovo_extractedsigs,
                    params.mutational_context
                )
                SigPA_Extracted(
                    HDP_combine.out.deNovo_extractedsigs,
                    params.mutational_matrix
                )
//                if (HDP_combine.out.deNovo_lowconfsigs != null){
//                    SigPlt_LowConfidence(
//                        HDP_combine.out.deNovo_lowconfsigs,
//                        params.mutational_context
//                    )
//                    SigPA_LowConfidence(
//                        HDP_combine.out.deNovo_lowconfsigs,
//                        params.mutational_matrix
//                    )
//                }
                } else {
                SigPlt_Extracted(
                    HDP_combine.out.deNovo_extractedsigs,
                    params.mutational_context
            )
//            if (HDP_combine.out.deNovo_lowconfsigs != null){
//               SigPlt_LowConfidence(
//                    HDP_combine.out.deNovo_lowconfsigs,
//                    params.mutational_context
//                )
//                }
            }             
        } else {
            if (params.decompose == true) {
            HDP_flat(
            params.mutational_matrix,
            params.hierarchy_matrix,
            params.prior_matrix,
            params.analysis_type, 
            params.burnin_iterations,
            params.posterior,
            params.posterior_space,
            params.n_iter,
            params.threshold
            )
            HDP_combine(
                HDP_flat.out.HDP_chains,
                params.mutational_matrix,
                params.hierarchy_matrix,
                params.threshold
            )
            SigPA_Extracted(
                HDP_combine.out.deNovo_extractedsigs,
                params.mutational_matrix
            )
//            if (HDP_combine.out.deNovo_lowconfsigs != null) {
//                SigPA_LowConfidence(
//                    HDP_combine.out.deNovo_lowconfsigs,
//                    params.mutational_matrix
//                )
//            }
        }
        }
    }
}