import sys
import argparse
from Bio import SeqIO

parser=argparse.ArgumentParser(description='Extract sequences from fasta file.')
parser.add_argument('--ids','-ids',required=True,help='List of sequence ids to extract from fasta.')
parser.add_argument('--fasta','-fasta',required=True,help='Fasta file.')
parser.add_argument('--outfile','-out',type=argparse.FileType('w',encoding='UTF-8'),required=True,help='Outfile.')
args=parser.parse_args()

with open(args.ids,'r') as f:
	id=f.readlines()

id=[x.strip() for x in id]

for seq_record in SeqIO.parse(args.fasta,"fasta"):
	for i in id:
		if i==seq_record.id:
			args.outfile.write(">%s\n%s\n" % (seq_record.id, seq_record.seq))
