#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
from Bio import Phylo
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace
from cStringIO import StringIO
import matplotlib
from matplotlib.colors import rgb2hex, LinearSegmentedColormap
import argparse
"""
view_trees.py
Takes Newick format trees (or the output of trees2newick -- chrom<tab>pos<tab>newick)
    from stdin and displays them. Requires BioPython and a working X11
    connection.
    
    
    NOTE: if ETE3 gives you an error importing the TextFace module,
    make sure you've installed Qt5 (sudo apt-get install python-qtpy).
    If Qt4 is still installed, you have to run sudo apt-get remove python-qt4
    to avoid conflicts.
    see here https://github.com/etetoolkit/ete/issues/195 
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pops", "-p", help=\
"A file mapping haplotype IDs to populations", required=False, default=None)
    parser.add_argument("--colors", "-c", help=\
"A file mapping population IDs to hex color codes", required=False, default=None)
    return parser.parse_args()

def layout_default(node):
    if node.is_leaf():
        f = AttrFace("name", fsize=20)
        faces.add_face_to_node(f, node, 0, position="branch-right")
        
def layout_popcols(node):
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
    
    grad = LinearSegmentedColormap.from_list('gradient2n', ['#003ac8', '#c80700'])
    
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
        
        if len(data) > 3:
            # Parse metadata in 4th field.
            meta = data[3].split(';')
            if meta[0] != "":
                max_pers = -1
                clades_pers = []
                for elt in meta:
                    leaves, persistence = elt.split(':')
                    persistence = int(persistence)
                    leaves_arr = leaves.split(',')
                    
                    if (len(leaves_arr) > 1):
                        if persistence > max_pers:
                            max_pers = persistence
                        clades_pers.append((leaves_arr, persistence))
                        
                for tup in clades_pers:
                    node = tree.get_common_ancestor(tup[0])
                    # Add persistence value to node.
                    ns = NodeStyle()
                    ns['size'] = int(round((tup[1]/max_pers)*40))
                    col_rgb = grad(tup[1]/max_pers)
                    ns['fgcolor'] = rgb2hex(col_rgb)
                    #ns['fgcolor'] = '#A8CFEA'
                    node.set_style(ns)
                    node.add_face(TextFace(str(tup[1]), fsize=15), column=0, position="branch-right")
            
            if len(data) > 4:
                der_allele = data[4]
                leaves_arr = der_allele.split(',')
                if (len(leaves_arr) > 1):
                    node = tree.get_common_ancestor(leaves_arr)
                    style = NodeStyle()
                    style['fgcolor'] = "#0f0f0f"
                    style['size'] = 0
                    style['vt_line_color'] = '#ff0000'
                    style['hz_line_color'] = '#ff0000'
                    style['vt_line_width'] = 8
                    style['hz_line_width'] = 8
                    style['vt_line_type'] = 0
                    style['hz_line_type'] = 0
                    node.set_style(style)
                
        ts = TreeStyle()
        ts.mode = 'c'      
        # Hide normal leaf names if we're showing population colors - it'll be
        # more useful to show stuff aligned at the tips in that case
        
        #if colors_given:
        if True:
            ts.layout_fn = layout_popcols
            ts.show_leaf_name = False
        else:
            # Otherwise, just increase the default font size.
            ts.layout_fn = layout_default
            ts.show_leaf_name = False
            
        tree.show(tree_style=ts)
        
        #tree = Phylo.read(StringIO(newick), "newick")
        #Phylo.draw(tree, str, True, True, None, None, {"S_Kusunda-2-2": "#FF0000"})
        #Phylo.draw(tree, label_colors=({"Altai-1": "#FF0000"}))

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
