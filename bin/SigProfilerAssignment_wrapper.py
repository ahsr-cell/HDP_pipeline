#!/usr/bin/env python

from SigProfilerAssignment import Analyzer as Analyze
from pathlib import Path
import argparse
import sys

def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--mutational_matrix', required = True, 
                        help="")
    parser.add_argument('--deNovoSignatures_matrix', required = True,
                        help="Specify path to mSigHdp de novo signature matrix")
    parser.add_argument('--output_directory',required=True,
                        help='Specify path to mSigHdp to output directory')
    return parser

def sigPA_decompose(mutational_matrix, deNovoSignatures_matrix,output_directory):
    output = "SigProfilerDecomposition" 
    Analyze.decompose_fit(mutational_matrix, output, deNovoSignatures_matrix, genome_build="GRCh38",
                           verbose=False, make_plots = True,
                           volume = Path(output_directory) / "pkl"
                           )
    print('SigProfilerAssignment completed. Output plots can be found in /SigProfilerDecomposition')

def main(args):
   sigPA_decompose(args.mutational_matrix, args.deNovoSignatures_matrix, args.output_directory)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))    
