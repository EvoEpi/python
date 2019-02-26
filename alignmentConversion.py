import os
import sys
import argparse
from Bio import AlignIO

parser=argparse.ArgumentParser(description='Extract sequences from fasta file.')
parser.add_argument('--infile','-i',
                    type=argparse.FileType('r',encoding='UTF-8'),
                    required=True,
                    help='Infile alignment.')
parser.add_argument('--in_format','-if',
                    required=True,
                    help='Infile format e.g., "fasta", "phylip", "clustal".')
parser.add_argument('--outfile','-o',
                    type=argparse.FileType('w',encoding='UTF-8'),
                    required=True,
                    help='Outfile alignment.')
parser.add_argument('--out_format','-of',
                    required=True,
                    help='Outfile format e.g., "fasta", "phylip", "clustal". FYI "phylip-relaxed" is supported by RAxML.')
args=parser.parse_args()

alignments=AlignIO.parse(args.infile,args.in_format)
AlignIO.write(alignments,args.outfile,args.out_format)

args.outfile.close()
args.infile.close()
