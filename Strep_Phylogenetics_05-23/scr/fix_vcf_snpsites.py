#!/usr/bin/env python
#
# Author: Martin Sikora, 2022-01-27

# fix_vcf_snpsites.py
#
# fix missing genotype coding in snp-sites output
#----------------------------------------------------------------------


from __future__ import print_function
import sys,re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-r", action="store", dest="ref_contig", default="1", type = str, help="reference contig name")
args = parser.parse_args()

def main():

    ## process vcf file
    for line in sys.stdin:

        ## print vcf header lines as is 
        if line.strip().startswith("#"):
            if line.strip().startswith("##contig="):
                print(re.sub(r"ID=.*,","ID="+args.ref_contig+",",line.strip()))
            print(line.strip())
            continue
        
        ## extract fields and replace ref contig name
        fields = line.strip().split()
        annot = fields[:9]
        annot[0] = args.ref_contig

        ## replace "*" allele genotypes with "."
        alleles = fields[4].split(",")  
        
        if "*" in alleles:
            idx = alleles.index("*")
            gt = ["." if f == str(idx + 1) else f for f in fields[9:]]
        else:
            gt = fields[9:]
            
        ## write output, 1 allele per
        print('\t'.join(annot), '\t'.join(gt), sep = '\t')

if __name__ == '__main__':
    main()
