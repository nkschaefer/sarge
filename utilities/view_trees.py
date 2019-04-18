#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
from Bio import Phylo
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle
from cStringIO import StringIO
import argparse
"""
view_trees.py
Takes Newick format trees (or the output of trees2newick -- chrom<tab>pos<tab>newick)
    from stdin and displays them. Requires BioPython and a working X11
    connection.
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pops", "-p", help=\
"A file mapping haplotype IDs to populations", required=False, default=None)
    parser.add_argument("--colors", "-c", help=\
"A file mapping population IDs to hex color codes", required=False, default=None)
    return parser.parse_args()

def layout(node):
    if node.is_leaf():
        f = AttrFace("name", fsize=20)
        faces.add_face_to_node(f, node, 0, position="aligned")
        
def main(args):
    """Main method"""
    
    # Check if the user has decided to color nodes by populations.
    colors_given = False
    pop2col = {}
    hap2col = {}
    nodestyles = {}
    options = parse_args()
    if (options.pops is not None and options.colors is None) or \
        (options.colors is not None and options.pops is None):
        print("ERROR: to color leaves by population, you must provide a file \
mapping haplotype IDs to population IDs, and another file mapping population \
IDs to hexadecimal color codes.", file=sys.stderr)
        exit(1)
    elif options.pops is not None and options.colors is not None:
        colors_given = True
        f = open(options.colors, 'r')
        for line in f:
            line = line.rstrip()
            if line != "":
                pop, col = line.split('\t')
                pop2col[pop] = col
        f.close()
        f = open(options.pops, 'r')
        for line in f:
            line = line.rstrip()
            if line != "":
                hap, pop = line.split('\t')
                if pop in pop2col:
                    hap2col[hap + "-1"] = pop2col[pop]
                    hap2col[hap + "-2"] = pop2col[pop]
                    ns1 = NodeStyle()
                    ns1['bgcolor'] = pop2col[pop]
                    ns2 = NodeStyle()
                    ns2['bgcolor'] = pop2col[pop]
                    nodestyles[hap + '-1'] = ns1
                    nodestyles[hap + '-2'] = ns2
        f.close()
    
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
        
        tree = Tree(newick)
        
        if colors_given:
            for leaf in tree:

                if leaf.name in nodestyles:
                    leaf.set_style(nodestyles[leaf.name])
        
        ts = TreeStyle()
        ts.mode = 'c'      
        ts.layout_fn = layout  
        ts.show_leaf_name = False
        tree.show(tree_style=ts)
        #tree = Phylo.read(StringIO(newick), "newick")
        #Phylo.draw(tree, str, True, True, None, None, {"S_Kusunda-2-2": "#FF0000"})
        #Phylo.draw(tree, label_colors=({"Altai-1": "#FF0000"}))

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
