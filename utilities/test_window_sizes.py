#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import gzip
import subprocess
import os
from os.path import dirname, abspath
import math
"""
Tests propagation distances to determine an appropriately large (but minimum) 
one.
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geno", "-g", help="The .geno.gz file for SARGE input.", required=True)
    parser.add_argument("--sites", "-s", help="The .sites file for SARGE input.", required=True)
    parser.add_argument("--min", "-m", help="The minimum window size to test", required=True, type=int)
    parser.add_argument("--max", "-M", help="The maximum window size to test (will \
break once an appropriate size is found)", required=True, type=int)
    parser.add_argument("--step", "-S", help="The number of bases to add for each trial", type=int, required=True)
    return parser.parse_args()

def main(args):
    """Main method"""
    # Parse arguments
    options = parse_args()
    
    # Get root directory
    sarge_root = dirname(dirname(abspath(__file__)))
    
    print("winsize\tmean\tsd")
    
    devnull = open(os.devnull, 'w')
    
    for winsize in range(options.min, options.max, options.step):
        p1 = subprocess.Popen(['{}/bin/sarge'.format(sarge_root), '-i', options.geno, \
            '-s', options.sites, '-p', str(winsize), '-o', '-', '-r', 'test.recomb'], \
            stdout=subprocess.PIPE, stderr=devnull)
        p2 = subprocess.Popen(['{}/bin/haplens'.format(sarge_root), '-o', '-', \
            '-p', str(winsize)], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=devnull)
        p3 = subprocess.Popen(['{}/bin/dump_haplens'.format(sarge_root)], stdin=p2.stdout,\
            stdout=subprocess.PIPE, stderr=devnull)
        p4 = subprocess.Popen(['head', '-2000'], stdin=p3.stdout, stdout=subprocess.PIPE)
        p5 = subprocess.Popen(['tail', '-1000'], stdin=p4.stdout, stdout=subprocess.PIPE)
        meansum = 0
        meantot = 0
        out, err = p5.communicate()
        nums = []
        for line in out.split('\n'):
            line = line.strip()
            if line != "":
                dat = line.split('\t')
                if dat[4] != "persistence":
                    meansum += int(dat[4])
                    meantot += 1
                    nums.append(int(dat[4]))
        if meantot == 0:
            meantot = 1
        avg = meansum/meantot
        varsum = 0 
        for num in nums:
            varsum += math.pow(num-avg, 2)
        print("{}\t{}\t{}".format(winsize, avg, math.sqrt(varsum/meantot)))
        
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

