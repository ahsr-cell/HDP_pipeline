#!/usr/bin/env python

import SigProfilerAssignment
import argparse
import sys

def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--mutational_matrix', required = True, 
                        help="")
    parser.add_argument('--deNovoSignatures_matrix', required = True,
                        help="Specify path to mSigHdp de novo signature matrix")
    return parser

def sigPA_decompose(mutational_matrix, deNovoSignatures_matrix):
    output = "SigProfilerDecomposition" 
    SigProfilerAssignment.Analyzer.decompose_fit(mutational_matrix, output, deNovoSignatures_matrix, genome_build="GRCh38",
                           verbose=False, make_plots = True)
    print('SigProfilerAssignment completed. Output plots can be found in /SigProfilerDecomposition')

def main(args):
   sigPA_decompose(args.mutational_matrix, args.deNovoSignatures_matrix)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))    
