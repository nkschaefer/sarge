#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import re
from collections import defaultdict
"""
sarge_relate_output.py
After running Relate on converted SARGE input files (use sarge_relate_input.py), 
print tree information to stdout in a format that can be piped to 
test/treedat2sarge to create SARGE-format output files.
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basename", "-b", help="Base name", required=True)
    return parser.parse_args()

def get_leaves_below(cur_node_ind, children, nhaps, leaves):
    if cur_node_ind < nhaps:
        leaves.append(str(cur_node_ind))
            
def print_tree(nodes, children, cur_node_ind, nhaps):
    tree = ""

    leaves_below = []
    get_leaves_below(cur_node_ind, children, nhaps, leaves_below)
    tree += ",".join(leaves_below) + "|"

    # Branch length
    tree += str(nodes[cur_node_ind]) + "|"
    
    if cur_node_ind in children:
        nchildren = 0
        for child_ind in children[cur_node_ind]:
            nchildren += 1
        
        tree += str(nchildren) + "|"
        for child_ind in children[cur_node_ind]:
            tree += print_tree(nodes, children, child_ind, nhaps)
    else:
        tree += "0|"
    
    return tree
        
def main(args):
    """Main method"""
    options = parse_args()
    
    # Get SNP IDs    
    snps = []
    chrom = None
    f = open('{}.haps'.format(options.basename), 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split()
        if chrom is None:
            chrom = dat[0]
        snps.append(int(dat[2]))
    f.close()
    
    prevpos = None
    treematch = re.compile(r'([\-0-9]+):\(([0-9\.]+) ([0-9\.]+) [0-9]+ [0-9]+\)')
    nhaps = None
    
    nodes = None
    children = None
    
    f = open('{}.anc'.format(options.basename), 'r')
    for ind, line in enumerate(f):
        line = line.rstrip()
        if ind == 0:
            nhaps = int(line.split()[1])
        if ind >= 2:
            
            pos = int(line.split(' ')[0][:-1])
            
            if prevpos is not None:
                prevtree = print_tree(nodes, children, len(nodes)-1, nhaps)
                for site_ind in xrange(prevpos, pos):
                    print("{}\t{}\t{}".format(chrom, snps[site_ind], prevtree))
            
            nodes = []
            children = defaultdict(list)
            for match in treematch.finditer(line):
                node_ind = len(nodes)
                nodes.append(float(match.group(2)))
                children[int(match.group(1))].append(node_ind)
                num_muts = float(match.group(3))
            
            prevpos = pos
            
    f.close()
    
    if prevpos is not None:
        prevtree = print_tree(nodes, children, len(nodes)-1, nhaps)
        for site_ind in xrange(prevpos, len(snps)-1):
            print("{}\t{}\t{}".format(chrom, snps[site_ind], prevtree))
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

