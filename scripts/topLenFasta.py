#!/usr/bin/env python
import argparse
from Bio import SeqIO
import pandas as pd
import re, gzip

parser = argparse.ArgumentParser(description='get sequence length from fasta file')
parser.add_argument('-i','--infile', dest='infile', metavar='<fastq>',
                    type=str, default="", help="input file (default: stdin)")
parser.add_argument('-o','--outfile', dest='outfile', metavar='<fastq>',
                    type=str, default="", help="output file")
parser.add_argument('-n','--topNumber', dest='ntop', metavar='int',
                    type=int, default=5, help="Number of top reads")
parser.add_argument('--split', action='store_true', dest='split',
                    help="split fasta in to single file")
args = parser.parse_args()

infile = args.infile
outfile = args.outfile
ntop = args.ntop
top = pd.DataFrame()
out = {}

if args.infile == "stdin":
    infile = sys.stdin
elif re.search('.gz$',args.infile):
    infile = gzip.open(args.infile,'rt')
elif re.search('.bgz$',args.infile):
    infile = gzip.open(args.infile,'rt')
else:
    infile = open(args.infile,"r")

for seq_record in SeqIO.parse(infile, "fastq"):
    if len(top) < ntop:
        if len(top) == 0:
            top = pd.DataFrame([len(seq_record)], index=[seq_record.id])
            out[seq_record.id] = seq_record
        else:
            top = pd.concat([top, pd.DataFrame([len(seq_record)], index=[seq_record.id])])
            out[seq_record.id] = seq_record
    else:
        minL = top.loc[top[0] == top.min(0)[0],]
        if minL[0][0] < len(seq_record):
            top = top.drop(minL.index[0])
            del(out[minL.index[0]])
            top = pd.concat([top, pd.DataFrame([len(seq_record)], index=[seq_record.id])])
            out[seq_record.id] = seq_record

print("get top {} reads ".format(ntop))

if args.split:

    i = 1
    for record in out.values():
        SeqIO.write([record], "{}/read_{}_{}.fasta".format(outfile, i,len(record)), "fasta")
        i+=1
else:
    SeqIO.write(out.values(), "{}/top{}.fasta".format(outfile, ntop), "fasta")

print("Done")
