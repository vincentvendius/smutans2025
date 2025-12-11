#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:03:03 2023

@author: vincentven

Takes in a sorted contig list with number of reads and damage estimates
prints out the contigs list with only unique contigs
and a metric for damage esmination
The calculation being done for each set of contig values:
#sum(n){number-of-reads__n*damage__n}/sum(n){number-of-reads__n}

"""
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py damagefile.tsv damagefilecorrected.tsv")
    sys.exit(1)


damagevector=open(sys.argv[1],"r")
outfile=open(sys.argv[2],"w")
#outvector=open("/home/vincentven/Downloads/testoutdmg.tsv","w")
contigdmg=dict()
contigreads=dict()
damsum=0
readsum=0
for line in damagevector:

    linelist=line[:-1].split("\t")
    contig=linelist[0]
    nreads=linelist[1]
    damage=linelist[2]
    if contig in contigdmg:        
        contigdmg[contig]+=float(damage)*int(nreads)
        contigreads[contig]+=int(nreads)
    else:
        contigdmg[contig]=float(damage)*int(nreads)
        contigreads[contig]=int(nreads)

newfile="contig\tdamage\n"
for key in contigdmg:
    finaldmg=contigdmg[key]/contigreads[key]
    newfile+=key+"\t"+str(finaldmg)+"\n"


outfile.write(newfile)

damagevector.close()
outfile.close()
