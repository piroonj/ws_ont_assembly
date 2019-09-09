#!/usr/bin/env python
from Bio import SeqIO
import sys

fname = sys.argv[1]
odir = sys.argv[2]
outFormat="fasta"
for record in SeqIO.parse(fname, 'fasta'):
    oname = "{:s}.{:d}.fasta".format(record.id,len(record.seq))
    SeqIO.write(record, "{}/{}".format(odir, oname), outFormat)


