#!/usr/bin/env python3

import sigProfilerPlotting as sigPlt
import argparse
import sys
from pathlib import Path

def get_arguments():
    parser = argparse.ArgumentParser(description = 'Get user defined SigProfilerPlotting options')
    parser.add_argument('--mutation_context', required = True, default = "SBS96",
                        help = 'Specify context of mutational matrix; options are SBS96 (default), SBS288, SBS1536, DBS78, or ID83.')
    parser.add_argument('--deNovoSignatures_matrix', required = True, 
                        help = 'Specify path to mSigHdp de novo signature matrix')
    parser.add_argument('--output_directory', required = True,
                        help="Specify output directory")
    return parser

def sigPlt_plotting(mutation_context, deNovoSignatures_matrix, output_directory):
  print('Using SigProfilerPlotting to visualise mSigHdp results')
  project = "SigProfilerPlots"
  matrix_path = deNovoSignatures_matrix

  if mutation_context == 'SBS96' or mutation_context == 'SBS288' or mutation_context == 'SBS1536': 
    if mutation_context == 'SBS96':
      u_plot_type = '96'
    
    if mutation_context == 'SBS288':
      u_plot_type = '288'
    
    if mutation_context == 'SBS1536':
      u_plot_type = '1536'    
    
    sigPlt.plotSBS(matrix_path, output_directory, project, u_plot_type, percentage=True,
                   volume = Path(output_directory) / "pkl")
    
  if mutation_context == 'DBS78':
      u_plot_type = '78'
      sigPlt.plotDBS(matrix_path, output_directory, project, u_plot_type, percentage=True,
                   volume = Path(output_directory) / "pkl")
  
  if mutation_context == 'ID83':
      u_plot_type = '83'
      sigPlt.plotID(matrix_path, output_directory, project, u_plot_type, percentage=True,
                   volume = Path(output_directory) / "pkl")
  
  print('SigProfilerPlotting completed. Output plots can be found in /Signature_Spectra')
  
def main(args):
   sigPlt_plotting(args.mutation_context, args.deNovoSignatures_matrix, args.output_directory)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))