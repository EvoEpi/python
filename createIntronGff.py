import os
import sys
import argparse
import gffutils

parser=argparse.ArgumentParser(description='Creates a gff file of introns.')
parser.add_argument('--in_gff','-i',type=str,
                                    required=True,
                                    help='A gff[3] file.')
parser.add_argument('--out_gff','-o',type=str,
                                    required=True,
                                    help='Filename of intron gff outfile.')
args=parser.parse_args()

"""
function humbly taken from Dr. Chad E. Niederhuth's GitHub
https://github.com/niederhuth/python-scripts-and-functions/blob/master/functions_plant.py
"""
def intron_gff(features,output=(),remove_db=True):
    db=gffutils.create_db(features, 'temp_gff.db')
    gff_out=gffutils.gffwriter.GFFWriter(output, in_place=True)
    name="none"
    count=1
    for rec in db.create_introns():
        if rec['Parent'] == name:
            count=count+1
            rec.attributes['ID'] = [str(rec['Parent'])[2:-2] + '.intron' + str(count)]
            gff_out.write_rec(rec)
        else:
            name = rec['Parent']
            count = 1
            rec.attributes['ID'] = [str(rec['Parent'])[2:-2] + '.intron' + str(count)]
            gff_out.write_rec(rec)
    gff_out.close()
    if remove_db:
        os.remove('temp_gff.db')

intron_gff(args.in_gff,args.out_gff)
