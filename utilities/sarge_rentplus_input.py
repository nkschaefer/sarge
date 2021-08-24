#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import gzip
"""
sarge_rentplus_input.py
Converts SARGE format input data to Rent+ format.
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    return parser.parse_args()

def main(args):
    """Main method"""
    if len(args) < 3:
        print("USAGE: sarge_rentplus_input.py basename out", file=sys.stderr)
        exit(1)
    
    basename = args[1]
    out = args[2]
    
    haps = []
    first = True
    
    sites = []
    s = open("{}.sites".format(basename), 'r')
    for line in s:
        line = line.rstrip()
        dat = line.split('\t')
        sites.append(dat[1])
    s.close()
    
    outf = open(out, "w")
    print(" ".join(sites), file=outf)
    
    g = gzip.open("{}.geno.gz".format(basename), 'r')
    for line in g:
        line = line.rstrip()
        if first:
            for ind in xrange(0, len(line)):
                haps.append([])
            first = False
        for ind, elt in enumerate(line):
            haps[ind].append(elt)
    
    g.close()
    
    for hap in haps:
        print("".join(hap), file=outf)
    
    outf.close()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

