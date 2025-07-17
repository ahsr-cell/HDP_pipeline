#!/usr/bin/env python

import sigProfilerPlotting as sigPlt
import argparse
import sys

def get_arguments():
    parser = argparse.ArgumentParser(description = 'Get user defined SigProfilerPlotting options')
    parser.add_argument('--mutational_context', required = True, default = "SBS96",
                        help = 'Specify context of mutational matrix; options are SBS96 (default), SBS288, SBS1536, DBS78, or ID83.')
    parser.add_argument('--deNovoSignatures_matrix', required = True, 
                        help = 'Specify path to mSigHdp de novo signature matrix')
    return parser

def sigPlt_plotting(mutational_context, deNovoSignatures_matrix):
  print('Using SigProfilerPlotting to visualise mSigHdp results')
  project = "SigProfilerPlots"
  output_path = "Signature_Spectra"
  matrix_path = deNovoSignatures_matrix

  if mutational_context == 'SBS96' or mutational_context == 'SBS288' or mutational_context == 'SBS1536': 
    if mutational_context == 'SBS96':
      u_plot_type = '96'
    
    if mutational_context == 'SBS288':
      u_plot_type = '288'
    
    if mutational_context == 'SBS1536':
      u_plot_type = '1536'    
  sigPlt.plotSBS(matrix_path, output_path, project, u_plot_type, percentage=False)
    
  if mutational_context == 'DBS78':
      u_plot_type = '78'
  sigPlt.plotDBS(matrix_path, output_path, project, u_plot_type, percentage=False)
  
  if mutational_context == 'ID83':
      u_plot_type = '83'
  sigPlt.plotID(matrix_path, output_path, project, u_plot_type, percentage=False)
  
  print('SigProfilerPlotting completed. Output plots can be found in /Signature_Spectra')
  
def main(args):
   sigPlt_plotting(args.mutational_context, args.deNovoSignatures_matrix)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))