#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { mSigHdp } from './workflows/mSigHdp.nf'
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
    // WORKFLOW: Full suite of analysis: mSigHdp, SigProfilerPlotting, and SigProfilerAssignment
    //
    if (params.plotting == true){
        if (params.decompose == true){
            mSigHdp(
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin_iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations,
                params.mutational_matrix
            )

            SigPlt_Extracted(
                mSigHdp.out.deNovo_extractedsigs,
                params.mutational_context
            )

            SigPA_Extracted(
                mSigHdp.out.deNovo_extractedsigs,
                params.mutational_matrix
            )
            if (mSigHdp.out.deNovo_lowconfsigs != null){
                SigPlt_LowConfidence(
                    mSigHdp.out.deNovo_lowconfsigs,
                    params.mutational_context
                )
                SigPA_LowConfidence(
                    mSigHdp.out.deNovo_lowconfsigs,
                    params.mutational_matrix
                )
            }
        }
        else {
            mSigHdp(
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin_iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations,
                params.mutational_matrix
            )
            SigPlt_Extracted(
                mSigHdp.out.deNovo_extractedsigs,
                params.mutational_context
            )
            if (mSigHdp.out.deNovo_lowconfsigs != null){
                SigPlt_LowConfidence(
                    mSigHdp.out.deNovo_lowconfsigs,
                    params.mutational_context
                )
            }
        }
    }
    else {
        if (params.decompose == true) {
            mSigHdp(
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin_iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations,
                params.mutational_matrix
            )
            SigPA_Extracted(
                mSigHdp.out.deNovo_extractedsigs,
                params.mutational_matrix
            )
            if (mSigHdp.out.deNovo_lowconfsigs != null){
                SigPA_LowConfidence(
                    mSigHdp.out.deNovo_lowconfsigs,
                    params.mutational_matrix
                )
            }
        }
        mSigHdp(
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin_iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations,
                params.mutational_matrix
            )
    } 
}