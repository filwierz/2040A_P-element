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
    finds diagnostic SNPs from a sync file generated from samtools' pileup file with the script mpileup2sync.jar (can be found here: https://sourceforge.net/projects/popoolation2/)""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""

Author
-------------
Filip Wierzbicki

Prerequisites
-------------
    python version 3+
    for input samtools and mpileup2sync.jar from popoolation2

""")

parser.add_argument('--sync', type=argparse.FileType('r'), default=None,dest="file", required=True, help="A sync file")
parser.add_argument('--min-cov', type=float, required=False, dest="mincoverage", default=0.0, help="min coverage")
parser.add_argument('--selfIDs', type=str, required=True, dest="sID", default=None, help="comma-separated, ordered integer IDs of samples that are the reference species")
parser.add_argument('--nonselfIDs', type=str, required=True, dest="nsID", default=None, help="comma-separated, ordered integer IDs of samples that are the non-reference species")

args = parser.parse_args()
sIDs=args.sID
nsIDs=args.nsID
mcov=args.mincoverage

s=sIDs.split(",")
n=nsIDs.split(",")

s1,s2,s3=[int(i) for i in s]
n1,n2,n3=[int(i) for i in n]

##note:
#../Canton-S.sort.bam ../Iso1.sort.bam ../Oregon-R.sort.bam Mod6.sort.bam.....
#Dmel_rhi	2	N	0:0:0:0:0:0	12:0:0:0:0:0	28:0:0:0:0:0	0:0:0:0:0:0	0:0:0:0:0:0	0:0:0:0:0:0
#column1 gene
#column2 position
#reference character if provided to samtools
#column>3: A-count:T-count:C-count:G-count:N-count:deletion-count

ptopr=collections.defaultdict((lambda:collections.defaultdict(lambda:[0,0])))

for line in args.file:
    line=line.rstrip("\n")
    a=line.split("\t")
    chro=a[0]
    pos=int(a[1])

    rSNP="XYZ"
    nRSNP="XYZ"

    rs=[]
    ns=[]

    for i in s1,s2,s3:
        j=i+2
        if(a[j]=="0:0:0:0:0:0"):
            break
        b=a[j].split(":")
        A,T,C,G,N,de=[int(f) for f in b]
        if (N>0) or (de>0):
            break
        cov=A+T+C+G
        if(cov<mcov):
            break
        if (A/cov==1):
            rSNP="A"
        if (T/cov==1):
            rSNP="T"
        if (C/cov==1):
            rSNP="C"
        if (G/cov==1):
            rSNP="G"
        rs.append(rSNP)
    

    for i in n1,n2,n3:
        j=i+2
        if(a[j]=="0:0:0:0:0:0"):
            break
        b=a[j].split(":")
        A,T,C,G,N,de=[int(f) for f in b]
        if (N>0) or (de>0):
            break
        cov=A+T+C+G
        if(cov<mcov):
            break
        if (A/cov==1):
            nrSNP="A"
        if (T/cov==1):
            nrSNP="T"
        if (C/cov==1):
            nrSNP="C"
        if (G/cov==1):
            nrSNP="G"
        ns.append(nrSNP)

    if rSNP=="XYZ" or nrSNP=="XYZ":
        continue
    if len(set(rs)) != 1:
        continue
    if len(set(ns)) != 1:
        continue

    if (rSNP!=nrSNP):
        print(chro,pos,rSNP,nrSNP)



 



