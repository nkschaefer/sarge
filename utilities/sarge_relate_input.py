#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import gzip
import random
"""
sarge_relate_input.py
Convert a set of SARGE input files to Relate input files.
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inbase", "-i", help="SARGE base file name", required=True)
    parser.add_argument("--outbase", "-o", help="Relate base (output) file name", required=True)
    return parser.parse_args()

def main(args):
    """Main method"""
    options = parse_args()
    
    if options.inbase == options.outbase:
        print("ERROR: must use different input and output base name", file=sys.stderr)
        exit(1)
    
    allelefile_given = False
    allelefile = None
    try:    
        allelefile = open('{}.alleles'.format(options.inbase), 'r')
        allelefile_given = True
    except:
        allelefile_given = False
    
    sitesfile = open("{}.sites".format(options.inbase), 'r')
    genofile = gzip.open("{}.geno.gz".format(options.inbase), 'r')
    
    snp_id = 1
    
    haps_out = open('{}.haps'.format(options.outbase), 'w')
    samps_out = open('{}.samples'.format(options.outbase), 'w')
    map_out = open('{}.map'.format(options.outbase), 'w')
    
    print("pos COMBINED_rate Genetic_Map", file=map_out)
    print("ID_1 ID_2 missing", file=samps_out)
    print("0 0 0.0", file=samps_out)
    hfile = open('{}.haps'.format(options.inbase), 'r')
    prevhap = None
    
    for ind, hline in enumerate(hfile):
        hline = hline.rstrip()
        if ind % 2 == 0:
            prevhap = hline
        else:
            print("{} {} 0.0".format(prevhap, hline), file=samps_out)
    samps_out.close()
    
    prevpos = None
    for gline in genofile:
        sline = sitesfile.next()
        gline = gline.rstrip()
        sline = sline.rstrip()
        
        anc = None
        der = None
        if allelefile_given:
            aline = allelefile.next()
            aline = aline.rstrip()
            anc, der = aline.split('\t')
        else:
            anc, der = random.sample(['A', 'C', 'G', 'T'], 2)
        
        chrom, pos = sline.split('\t')
        if pos != prevpos:
            print("{} SNP{} {} {} {} ".format(chrom, snp_id, pos, anc, der) + " ".join(list(gline)), file=haps_out)
            print("{} 1.0 {}".format(pos, int(pos)/1e6), file=map_out)
            snp_id += 1
        prevpos = pos
    
    sitesfile.close()
    
    if allelefile_given:
        allelefile.close()
    
    haps_out.close()
    map_out.close()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

