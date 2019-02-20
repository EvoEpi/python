import os
import sys
import argparse
import pandas as pd
import pybedtools as pbt
import numpy as np
pd.options.mode.chained_assignment=None

parser=argparse.ArgumentParser(description='Under construction.')
parser.add_argument('--allc','-i',type=argparse.FileType('r',encoding='UTF-8'),
                                    required=True,
                                    help='Allc file generated (preferably) by methylpy.')
parser.add_argument('--gff','-gf',type=argparse.FileType('r',encoding='UTF-8'),
                                    required=True,
                                    help='A gff of the genome WGBS reads were mapped to.')
parser.add_argument('--taxa','-t',type=str,
                                    required=True,
                                    help='"plant" (CG, CHG, CHH) or "animal" (CG, CH)')
parser.add_argument('--genome','-ge',type=str,
                                    required=True,
                                    help='Name of a tab separated text file containing the length of each contig/scaffold/chromosome.')
parser.add_argument('--feature','-fe',type=str,
                                    required=True,
                                    help='Name of feature in gff.')
parser.add_argument('--flanking','-fl',type=int,
                                    required=True,
                                    help='Number of base pairs of flanking.')
parser.add_argument('--windows','-w',type=int,
                                    required=True,
                                    help='Number of windows per upstream of feature, within feature, and downstream of feature.')
parser.add_argument('--outfile','-o',type=argparse.FileType('w', encoding='UTF-8'),
                                    required=True,
                                    help='Outfile.')
args=parser.parse_args()

def allc_to_bed(d,bed=True):
    if args.taxa=='animal':
        d['con']=d['con'].replace(regex={r'CG[A|C|G|T|N]':'CG',r'C[A|C|G|T|N][A|C|G|T|N]':'CH'})
    elif args.taxa=='plant':
        d['con']=d['con'].replace(regex={r'CG[A|C|G|T|N]':'CG',r'C[A|C|T|N]G':'CHG',r'C[A|C|T|N][A|C|T|N]':'CHH'})
    new_col=d['pos']
    d.insert(loc=2,column='pos2',value=new_col)
    if bed is True:
        nd=pbt.BedTool.from_dataframe(d)
    return nd

def gff_to_bed(d,feature,strand,bed=True):
    nd=d.loc[(df_gff['type']==feature) & (df_gff['strand']==strand)]
    nd=nd[['seqid','start','end','attributes']]
    if (bed is True and strand=='+'):
        nd=pbt.BedTool.from_dataframe(nd)
        feat_win_pd=pbt.BedTool.window_maker(nd,b=nd,n=args.windows,i='winnum').to_dataframe()
        us=pbt.BedTool.flank(nd,g=args.genome,l=args.flanking,r=0)
        ds=pbt.BedTool.flank(nd,g=args.genome,l=0,r=args.flanking)
        us_win=pbt.BedTool.window_maker(us,b=us,n=args.windows,i='winnum')
        ds_win=pbt.BedTool.window_maker(ds,b=ds,n=args.windows,i='winnum')
        us_win_pd=pbt.BedTool.to_dataframe(us_win)
        ds_win_pd=pbt.BedTool.to_dataframe(ds_win)
        feat_win_pd['name']=feat_win_pd['name'].apply(lambda x: x+args.windows)
        ds_win_pd['name']=ds_win_pd['name'].apply(lambda x: x+(args.windows*2))
        return {'us_win_pd_pos':us_win_pd,'feat_win_pd_pos':feat_win_pd,'ds_win_pd_pos':ds_win_pd}
    elif (bed is True and strand=='-'):
        nd=pbt.BedTool.from_dataframe(nd)
        feat_win_pd=pbt.BedTool.window_maker(nd,b=nd,n=args.windows,i='winnum',reverse=True).to_dataframe()
        ds=pbt.BedTool.flank(nd,g=args.genome,l=args.flanking,r=0)
        us=pbt.BedTool.flank(nd,g=args.genome,l=0,r=args.flanking)
        us_win=pbt.BedTool.window_maker(us,b=us,n=args.windows,i='winnum',reverse=True)
        ds_win=pbt.BedTool.window_maker(ds,b=ds,n=args.windows,i='winnum',reverse=True)
        us_win_pd=pbt.BedTool.to_dataframe(us_win)
        ds_win_pd=pbt.BedTool.to_dataframe(ds_win)
        feat_win_pd['name']=feat_win_pd['name'].apply(lambda x: x+args.windows)
        ds_win_pd['name']=ds_win_pd['name'].apply(lambda x: x+(args.windows*2))
        return {'us_win_pd_neg':us_win_pd,'feat_win_pd_neg':feat_win_pd,'ds_win_pd_neg':ds_win_pd}

def methylation_stats(d):
    nd={}
    nd['bin']=d['bin'].unique()
    nd['con']=d['con'].unique()
    nd['wei']=d['mec'].sum()/d['tot'].sum()
    return pd.Series(nd, index=['wei'])

df_allc=pd.read_csv(args.allc,
                sep='\t',
                names=['chr','pos','str','con','mec','tot','binom'])

df_gff=pd.read_csv(args.gff,
                sep='\t',
                names=['seqid','source','type','start','end','score','strand','phase','attributes'])

bed_pos=gff_to_bed(df_gff,feature=args.feature,strand='+',bed=True)

bed_neg=gff_to_bed(df_gff,feature=args.feature,strand='-',bed=True)

bed_allc=allc_to_bed(df_allc,bed=True)

concat=pd.concat([bed_pos['us_win_pd_pos'],bed_neg['us_win_pd_neg'],
                bed_pos['feat_win_pd_pos'],bed_neg['feat_win_pd_neg'],
                bed_pos['ds_win_pd_pos'],bed_neg['ds_win_pd_neg']],
                ignore_index=True)

concat_bed=pbt.BedTool.from_dataframe(concat).sort()

mapping=pbt.BedTool.intersect(bed_allc,concat_bed,wa=True,wb=True).to_dataframe().rename(index=str,columns={'chrom':'chr','start':'start',
                                                                                                            'end':'end','name':'strand',
                                                                                                            'score':'con','strand':'mec',
                                                                                                            'thickStart':'tot',
                                                                                                            'thickEnd':'binom',
                                                                                                            'itemRgb':'chrFeat',
                                                                                                            'blockCount':'startFeat',
                                                                                                            'blockSizes':'endFeat',
                                                                                                            'blockStarts':'bin'})

mapping_stats=mapping.groupby(['bin','con'],as_index=False).apply(methylation_stats).reset_index()

mapping_stats.to_csv(args.outfile,sep='\t',index=False)
