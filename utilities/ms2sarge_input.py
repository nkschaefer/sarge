#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import subprocess
import numpy
import argparse
from collections import defaultdict, OrderedDict
import re
"""
Given output from an MS simulation, converts to SARGE input format.
"""

def parse_args():
    parser = argparse.ArgumentParser(description="__doc__")
    parser.add_argument("--ms", "-m", help="The output from MS, with -T specified", required=True)
    parser.add_argument("--output_base", "-o", help="The base name for output files", required=True)
    parser.add_argument("--num_bases", "-n", help="Number of bases in MS simulation", required=True, type=int)
    return parser.parse_args()

def main(args):
    """Main method"""
    options = parse_args()
    basename = options.output_base
    
    indvfile = open('{}.haps'.format(basename), 'w')
   
    # Regular expression to match tree lines
    treeline = re.compile(r'^\[([0-9]+)\]([\(\)\,:0-9;\.]+)$')
    
    # Regular expression to match "segsites" line
    segsitesline = re.compile(r'^segsites: ([0-9]+)$')
    
    # Regular expression to match variant lines
    varline = re.compile(r'^[01]+$')
    
    # Regular expression to match "positions" line
    posline = re.compile(r'^positions:')
    
    inf = open(options.ms, 'r')
    
    tree_pos = 0
    
    positions = []
    positions_inverted = {}
    
    seqname = "seq"
    
    locsfile = open('{}.sites'.format(basename), 'w')
    
    variants = []
    
    # Store trees, keyed by positions at which they are valid.
    trees_pos = OrderedDict()
    
    samp_idx = 1
    hap_idx = 1

    for line in inf:
        line = line.rstrip()
        # Determine what type of line this is.
        if treeline.match(line):
            pass        
        elif segsitesline.match(line):
            match = segsitesline.match(line)
        elif posline.match(line):
            #match = posline.match(line)
            #positions_text = match.group(1).split()
            positions_text = line.split()[1:]
            for pos in positions_text:
                pos = float(pos)
                pos_bases = int(pos * options.num_bases)
                positions.append(pos_bases)
                positions_inverted[pos_bases] = len(positions)-1
                print("{}\t{}".format(seqname, pos_bases), file=locsfile)

        elif varline.match(line):
            # Rows: haplotypes
            # Columns: sites
            match = varline.match(line)
            variants.append(list(match.group(0)))
            
            print("Human-{}-{}".format(samp_idx, hap_idx), file=indvfile)
            hap_idx += 1
            if hap_idx > 2:
                hap_idx = 1
                samp_idx += 1
    
    indvfile.close()

    # We now have all variants. Transpose and write to file.
    variants = numpy.array(variants).T
    
    argfile = open("{}.geno".format(basename), 'w')
    for siteIndex in xrange(0, variants.shape[0]):
        print("".join(variants[siteIndex,:]), file=argfile)
    
    argfile.close()
    inf.close()
    locsfile.close()
    
    subprocess.call(['gzip', '{}.geno'.format(basename)])
    
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

