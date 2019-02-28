#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import gzip
"""
Given SARGE input files and a propagation distance, finds all gaps (i.e. centromeres)
where the files can be safely split into multiple files to run in parallel.
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prop_dist", "-p", type=int, help="Propagation distance in bp", required=True)
    parser.add_argument("--geno", "-g", help="SARGE input genotype file", required=True)
    parser.add_argument("--sites", "-s", help="SARGE input sites file", required=True)
    parser.add_argument("--output_prefix", "-o", help="base name for new output files to create", required=True)
    return parser.parse_args()

def main(args):
    """Main method"""
    options = parse_args()
    f = open(options.sites, 'r')
    fg = gzip.open(options.geno, 'r')
    
    out_index = 1
    geno_out = gzip.open("{}.{}.geno.gz".format(options.output_prefix, out_index), "w")
    sites_out = open("{}.{}.sites".format(options.output_prefix, out_index), "w")
    
    prevpos = None
    for line in f:
        line = line.rstrip()
        chrom, pos = line.split('\t')
        pos = int(pos)
        geno_line = fg.next()
        geno_line = geno_line.rstrip()
        if prevpos is not None and pos - prevpos > 6*options.prop_dist:
            # Time for a new output file
            out_index += 1
            geno_out.close()
            sites_out.close()
            geno_out = gzip.open("{}.{}.geno.gz".format(options.output_prefix, out_index), "w")
            sites_out = open("{}.{}.sites".format(options.output_prefix, out_index), "w")
        
        # Print stuff to current file.
        print(geno_line, file=geno_out)
        print(line, file=sites_out)
        prevpos = pos
        
    f.close()
    fg.close()
    
    geno_out.close()
    sites_out.close()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

