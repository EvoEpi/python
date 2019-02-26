import os
import sys
import argparse
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio import SeqIO

parser=argparse.ArgumentParser(description='Search for motif and reverse complement of motif in fasta file and record first and last bp position.')
parser.add_argument('--motif','-m',type=str,required=True,help='Motif fasta file')
parser.add_argument('--search','-s',required=True,help='Sequence fasta file to search.')
parser.add_argument('--outfile','-o',type=argparse.FileType('w',encoding='UTF-8'),required=True,help='Outfile (will be in bed format).')
args=parser.parse_args()

d=[]
for seq_motif in SeqIO.parse(args.motif,"fasta"):
	for seq in SeqIO.parse(args.search,"fasta"):
		results=SeqUtils.nt_search(str(seq.seq),seq_motif.seq)
		results_rc=SeqUtils.nt_search(str(seq.seq),seq_motif.seq.reverse_complement())
		for i in results[1:]:
			d.append({'search':seq.id,'motif':results[0],'first':i+1,'last':(i)+len(results[0])})
		for i in results_rc[1:]:
			d.append({'search':seq.id,'motif':results_rc[0],'first':i+1,'last':(i)+len(results_rc[0])})
df=pd.DataFrame(d)
df=df[['search','first','last','motif']]
df.to_csv(args.outfile,sep='\t',index=False)
