#!/usr/bin/env python
import os
import sys
import re
import argparse
import random
import collections

parser = argparse.ArgumentParser(description="""           
Description
-----------
    screens samples for diagnostic SNPs from a sync file generated from samtools' pileup file with the script mpileup2sync.jar (can be found here: https://sourceforge.net/projects/popoolation2/)""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""

Author
-------------
Filip Wierzbicki

Prerequisites
-------------
    python version 3+
    for input samtools and mpileup2sync.jar from popoolation2
    diagnostic SNPs from diagnosticSNPs-finder.py

""")

parser.add_argument('--sync', type=argparse.FileType('r'), default=None,dest="file", required=True, help="A sync file")
parser.add_argument('--snp', type=argparse.FileType('r'), default=None,dest="diasnp", required=True, help="diagnostic SNPs")
parser.add_argument('--min-cov', type=float, required=False, dest="mincoverage", default=0.0, help="min coverage")
parser.add_argument('--min-freq', type=float, required=False, dest="minfreq", default=0.99, help="min freq")


args = parser.parse_args()
mcov=args.mincoverage
mfreq=args.minfreq

toit=collections.defaultdict((lambda:collections.defaultdict(lambda:[str])))

for line in args.file:
    ##Dmel_rhi	2	N	0:0:0:0:0:0	
    line=line.rstrip("\n")
    a=line.split("\t")
    chro=a[0]
    pos=int(a[1])
    if(a[3]=="0:0:0:0:0:0"):
        continue
    b=a[3].split(":")
    A,T,C,G,N,de=[int(f) for f in b]
    cov=A+T+C+G
    if(cov<mcov):
        continue
    elif (A/cov>=mfreq):
        toit[chro][pos][0]="A"
    elif (T/cov>=mfreq):
        toit[chro][pos][0]="T"
    elif (C/cov>=mfreq):
        toit[chro][pos][0]="C"
    elif (G/cov>=mfreq):
        toit[chro][pos][0]="G"
    else:
        toit[chro][pos][0]="NA"
    

for line in args.diasnp:
    #Dsim_rhi 35 A T
    line=line.rstrip("\n")
    a=line.split(" ")
    chro=a[0]
    pos=int(a[1])
    ref=a[2]
    alt=a[3]
    for c,tmp in toit.items():
        if c==chro:
            for p,ch in tmp.items():
                if p==pos:
                    if ref==ch[0]:
                        print(chro,pos,ref,alt,ch[0],"clean")
                    elif alt==ch[0]:
                        print(chro,pos,ref,alt,ch[0],"contamination")
                    else:
                        print(chro,pos,ref,alt,ch[0],"error/polymorph")

