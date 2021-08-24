#! /usr/bin/env python3
import sys
import argparse
import tsinfer
import gzip
import subprocess
import random
import time
import tsdate
"""
sarge_tsinfer.py

Takes SARGE-format input files, runs tsinfer, and prints information to stdout.
This can then be piped to test/treedat2sarge to create SARGE-format output files
for analysis.
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    return parser.parse_args()

def print_tree(tree, node):
    clade = []
    
    for leaf in tree.leaves(node):
        if not tree.is_internal(leaf):
            clade.append(str(leaf))
    
    print(",".join(clade), end="|")
    # Add branch length
    print(tree.branch_length(node), end="|")
    numchildren = 0
    for child in tree.children(node):
        numchildren += 1
    print(numchildren, end="|")
    for child in tree.children(node):
        print_tree(tree, child)
       
def print_tree_wrapper(tree):
    # Get root haplotype name
    root = None
    nodes = tree.nodes(order='preorder')
    for n in nodes:
        root = n
        break
    print_tree(tree, root)
    print("")
    
def parse_args():
    parser = argparse.ArgumentParser(description="Run tsinfer on SARGE input files, and print \
tree information to stdout that can be piped to test/treedat2sarge to create SARGE-format output files.")
    parser.add_argument("--basename", "-b", help="The base file name of all SARGE input files", \
        required=True)
    parser.add_argument("--N", "-n", help="Effective population size (for tsdate)", required=True, type=int)
    parser.add_argument("--mu", "-m", help="Mutation rate (for tsdate)", required=True, type=float)
    return parser.parse_args()

def main(args):
    """Main method"""
    
    options = parse_args()
    basename = options.basename
    
    sites_f = open("{}.sites".format(basename), "r")
    geno_f = gzip.open("{}.geno.gz".format(basename), "r")
    alleles_f = None
    has_alleles = False
    try:
        alleles_f = open("{}.alleles".format(basename), "r")
        has_alleles = True
    except:
        has_alleles = False
        
    p_out = subprocess.check_output(['tail', '-1', '{}.sites'.format(basename)])

    seqlen = int(p_out.decode('ascii').rstrip().split('\t')[1])
    print("seq length: {}".format(seqlen), file=sys.stderr)
    seqlen += 1
    sampdat = tsinfer.SampleData(sequence_length=seqlen)

    sites_list = []

    print("Loading data...", file=sys.stderr)

    prevpos = None

    for geno_line in geno_f:
        sites_line = sites_f.readline()
        
        alleles_line = None
        if has_alleles:
            alleles_line = alleles_f.readline()
            alleles_line = alleles_line.rstrip()
            
        geno_line = geno_line.decode('ascii').rstrip()
        sites_line = sites_line.rstrip()
        
    
        chrom, pos = sites_line.split('\t')
    
        pos = int(pos)
        
        alleles = ['A', 'C', 'G', 'T']
        
        anc, der = random.sample(alleles, 2)
        
        if has_alleles:
            anc, der = alleles_line.split('\t')
        
        genotypes = list(geno_line)
    
        for ind, elt in enumerate(genotypes):
            genotypes[ind] = int(elt)
    
        if pos != prevpos:
            sampdat.add_site(pos, genotypes, [anc, der])
            sites_list.append([chrom, pos])
    
        prevpos = pos
    
    sampdat.finalise()
    out_time = open("{}_tsinfer.time".format(basename), 'w')
    print("Inferring trees...", file=sys.stderr)
    t1 = time.time()
    treeseq = tsinfer.infer(sampdat)
    treeseq = treeseq.simplify(keep_unary=False)
    treeseq = tsdate.date(treeseq, Ne=options.N, mutation_rate=options.mu)
    t2 = time.time()
    print("Done", file=sys.stderr)
    print(t2-t1, file=out_time)
    out_time.close()
    
    for ind, tree in enumerate(treeseq.trees()):
        start, end = tree.interval
        while len(sites_list) > 0 and sites_list[0][1] >= start and sites_list[0][1] < end:
            print("{}\t{}\t".format(sites_list[0][0], sites_list[0][1]), end="")
            print_tree_wrapper(tree)
            del sites_list[0]
        
    sites_f.close()
    geno_f.close()
    if has_alleles:
        alleles_f.close()
    

if __name__ == "__main__" :
    sys.exit(main(sys.argv))

