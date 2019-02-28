#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
from Bio import Phylo
from cStringIO import StringIO
"""
view_trees.py
Takes Newick format trees (or the output of trees2newick -- chrom<tab>pos<tab>newick)
    from stdin and displays them. Requires BioPython and a working X11
    connection.
"""

def main(args):
    """Main method"""
    for line in sys.stdin:
        line = line.rstrip()
        data = line.split('\t')
        newick = ""
        if len(data) == 1:
            newick = data[0]
        else:
            newick = data[2]
        
        if len(data) > 1:
            print("{}\t{}".format(data[0], data[1]), file=sys.stderr)
        
        tree = Phylo.read(StringIO(newick), "newick")
        Phylo.draw(tree)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
