import os
import sys
import argparse
import itertools
import pandas as pd
import pybedtools as pbt
import numpy as np
import statsmodels.api as sm
from functools import reduce
from scipy import stats
from collections import defaultdict
pd.options.mode.chained_assignment=None

parser=argparse.ArgumentParser(description='From an allc-gene intersect, gather methylation statistics.')
parser.add_argument('--allc','-a',type=argparse.FileType('r',encoding='UTF-8'),
                                    required=True,
                                    help='Allc (vcf-like).')
parser.add_argument('--bed','-b',type=argparse.FileType('r',encoding='UTF-8'),
                                    required=True,
                                    help='Four column bed file of attribute coordinates and ID.')
parser.add_argument('--taxa','-t',type=str,
                                    required=True,
                                    help='"animal" (CG, CH) or "plant" (CG, CHG, CHH)')
parser.add_argument('--sites','-s',type=int,
                                    required=True,
                                    help='Minumim number of cytosine sites.')
parser.add_argument('--coverage','-c',type=int,
                                    required=True,
                                    help='Per-site sequencing coverage threshold.')
parser.add_argument('--qvalue','-q',type=float,
                                    required=True,
                                    help='Q value threshold.')
parser.add_argument('--outfile','-o',type=argparse.FileType('w', encoding='UTF-8'),
                                    required=True,help='Outfile.')
args=parser.parse_args()

def allc_to_bed(d,bed=True):
    if args.taxa=='animal':
        d['con']=d['con'].replace(regex={r'CG[A|C|G|T|N]':'CG',r'C[A|C|T|N][A|C|G|T|N]':'CH'})
    elif args.taxa=='plant':
        d['con']=d['con'].replace(regex={r'CG[A|C|G|T|N]':'CG',r'C[A|C|T|N]G':'CHG',r'C[A|C|T|N][A|C|T|N]':'CHH'})
    new_col=d['pos']
    d.insert(loc=2,column='pos2',value=new_col)
    if bed is True:
        nd=pbt.BedTool.from_dataframe(d)
    return nd

def weighted_meth(d):
    return d['mec'].sum()/d['tot'].sum()

def methylation_stats(x):
    d={}
    d['con']=x['con'].unique()
    d['n']=x['start'].count()
    d['mes']=x['binom'].sum()
    d['mer']=x['mec'].sum()
    d['tor']=x['tot'].sum()
    d['wei']=x['mec'].sum()/x['tot'].sum()
    d['pval']=float(stats.binom.sf(d['mes']-1,d['n'],baseline[d['con']]))
    return pd.Series(d, index=['n','mes','mer','tor','wei','pval'])

def fdr(d,context):
    nd=d[(d['con']==context)]
    ar_gm=np.array(nd['att'][(nd['con']==context)])
    ar_pval=np.array(nd['pval'][(nd['con']==context)])
    mask=np.isfinite(ar_pval)
    pval_corrected=np.full(ar_pval.shape,np.nan)
    pval_corrected[mask]=sm.stats.multipletests(ar_pval[mask],
                                                alpha=0.05,method='fdr_bh',
                                                is_sorted=False)[1]
    equiv=dict(itertools.zip_longest(ar_gm,pval_corrected,fillvalue=None))
    nd["fdr"]=nd['att'].map(equiv)
    return nd

def bin(d,taxa,sites,coverage,qvalue):
    if taxa=='plant':
        conditions=[(d['fdr_cg']<=qvalue) & (d['fdr_chg']>qvalue) & (d['fdr_chh']>qvalue) & (d['n_cg']>=sites) & (d['tor_cg']/d['n_cg']>=coverage),
                    (d['fdr_chg']<=qvalue) & (d['fdr_chh']>qvalue) & (d['n_chg']>=sites) & (d['tor_chg']/d['n_chg']>=coverage),
                    (d['fdr_chh']<=qvalue) & (d['n_chh']>=sites) & (d['tor_chh']/d['n_chh']>=coverage),
                    (d['mes_cg']+d['mes_chg']+d['mes_chh'])==0 & (d['n_cg']+d['n_chg']+d['n_chh']>=sites)]
        choices=['gbM','CHG','CHH','UM']
        d['bin']=np.select(conditions, choices, default='NA')
        return d
    elif taxa=='animal':
        conditions=[(d['fdr_cg']<=qvalue) & (d['fdr_chg']>qvalue) & (d['fdr_chh']>qvalue) & (d['n_cg']>=sites) & (d['tor_cg']/d['n_cg']>=coverage),
                    (d['fdr_ch']<=qvalue) & (d['n_ch']>=sites) & (d['tor_chh']/d['n_chh']>=coverage),
                    (d['mes_cg']+d['mes_ch'])==0 & (d['n_cg']+d['n_ch']>=sites)]
        choices=['CG','CH','UM']
        d['bin']=np.select(conditions, choices, default='NA')
        return d

df_allc=pd.read_csv(args.allc,
                sep='\t',
                names=['chr','pos','str','con','mec','tot','binom'])

bed_allc=allc_to_bed(df_allc,bed=True)

df_bed=pd.read_csv(args.bed,
                sep='\t',
                names=['chr','start','end','cds'])

bed=pbt.BedTool.from_dataframe(df_bed)

df_intersect=pbt.BedTool.intersect(bed_allc,bed,wa=True,wb=True).to_dataframe().rename(index=str,columns={'chrom':'chr',
                                                                                                        'start':'start',
                                                                                                        'end':'end',
                                                                                                        'name':'strand',
                                                                                                        'score':'con',
                                                                                                        'strand':'mec',
                                                                                                        'thickStart':'tot',
                                                                                                        'thickEnd':'binom',
                                                                                                        'itemRgb':'chrAtt',
                                                                                                        'blockCount':'startAtt',
                                                                                                        'blockSizes':'endAtt',
                                                                                                        'blockStarts':'att'})

if args.taxa=="animal":
    baseline=df_intersect.groupby(['con']).apply(weighted_meth)
    df_stats=df_intersect.groupby(['att','con'],as_index=False).apply(methylation_stats).reset_index()
    df_cg=fdr(df_stats,'CG')
    df_ch=fdr(df_stats,'CH')
    df_cons=[df_cg,df_ch]
    df_final=reduce(lambda left,right: pd.merge(left,right,on='att'), df_cons)
    df_final.columns=['att','con_cg','n_cg','mes_cg','mer_cg','tor_cg','wei_cg','pval_cg','fdr_cg',
                        'con_ch','n_ch','mes_ch','mer_ch','tor_ch','wei_ch','pval_ch','fdr_ch']
    bin(df_final,args.taxa,args.sites,args.coverage,args.qvalue)
    df_final.to_csv(args.outfile, sep='\t', index=False)
elif args.taxa=="plant":
    baseline=df_intersect.groupby(['con']).apply(weighted_meth)
    df_stats=df_intersect.groupby(['att','con'],as_index=False).apply(methylation_stats).reset_index()
    df_cg=fdr(df_stats,'CG')
    df_chg=fdr(df_stats,'CHG')
    df_chh=fdr(df_stats,'CHH')
    df_cons=[df_cg,df_chg,df_chh]
    df_final=reduce(lambda left,right: pd.merge(left,right,on='att'), df_cons)
    df_final.columns=['att','con_cg','n_cg','mes_cg','mer_cg','tor_cg','wei_cg','pval_cg','fdr_cg',
                        'con_chg','n_chg','mes_chg','mer_chg','tor_chg','wei_chg','pval_chg','fdr_chg',
                        'con_chh','n_chh','mes_chh','mer_chh','tor_chh','wei_chh','pval_chh','fdr_chh']
    bin(df_final,args.taxa,args.sites,args.coverage,args.qvalue)
    df_final.to_csv(args.outfile, sep='\t', index=False)
