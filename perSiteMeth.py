import os
import sys
import argparse
import pandas as pd
pd.options.mode.chained_assignment=None

parser=argparse.ArgumentParser(description='From an allc file, gather per-site methylation levels of methylated C\'s and \
                                            output to separate files by context.')
parser.add_argument('--infile','-in',type=argparse.FileType('r',encoding='UTF-8'),
                                    required=True,
                                    help='Allc file (vcf-like).')
parser.add_argument('--coverage','-cov',type=int,
                                    required=True,
                                    help='Per-site sequencing coverage threshold.')
parser.add_argument('--taxa','-taxa',type=str,
                                    required=True,
                                    help='"plant" (CG, CHG, CHH) or "animal" (CG, CH)')
args=parser.parse_args()

df=pd.read_csv(args.infile,
                sep='\t',
                names=['chr','pos','str','con','mec','tot','bin'],
                dtype={'chr':str,'pos':int,'str':str,'con':str,'mec':int,'tot':int,'bin':int})

def per_site_mC(d,cov,context):
    for i in context:
        nd=d[(d['con']==i) & (d['tot']>=cov) & (d['bin']==1)]
        nd['per_site_wei']=nd['mec']/nd['tot']
        nd=nd[['chr','pos','str','per_site_wei']]
        nd.to_csv(i+'_per_site_methylation.tsv', sep='\t', index=False)

if args.taxa=="plant":
    df['con']=df['con'].replace(regex={r'CG[A|C|G|T|N]':'CG',
                                        r'C[A|C|T|N]G':'CHG',
                                        r'C[A|C|T|N][A|C|T|N]':'CHH'})
    per_site_mC(df,args.coverage,['CG','CHG','CHH'])
elif args.taxa=="animal":
    df['con']=df['con'].replace(regex={r'CG[A|C|G|T|N]':'CG',
                                        r'C[A|C|T|N][A|C|T|N]':'CH'})
    per_site_mC(df,args.coverage,['CG','CH'])
