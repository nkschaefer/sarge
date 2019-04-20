#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
from scipy.special import gegenbauer, hyp2f1, binom
import scipy
import numpy
import math
from multiprocessing import Pool
import argparse
"""
calc_ages_freqs.py

Given a population size (in haplotypes, i.e. 2N if diploid) and a number of bins, computes a matrix
of probabilities for frequencies of alleles at various frequencies (rows)
after various numbers of generations (as percentages of 2*TMRCA with outgroup).
Writes to a file that can be read in and used by the selection scan component.
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pop_size", "-N", help=\
        "The effective population size of the population. Report as total number of chromosomes, \
not individuals. In other words, if Ne in humans is 12,000, report 24000 here.",
        type=int,
        required=True)
    parser.add_argument("--sample_size", "-n", help=\
        "The number of chromosomes in your ARG data set. In other words, if you have 300 \
humans, then report 600 here.",
        type=int,
        required=True)
    #parser.add_argument("--freq_bins", "-fb", help=\
    #    "The number of bins into which to split the frequencies",\
    #    type=int,
    #    required=False,
    #    default=50)
    parser.add_argument("--gen_bins", "-gb", help=\
        "The number of bins into which to split the TMRCAs (generations in the past)",\
        type=int,
        required=False,
        default=50)
    parser.add_argument("--threads", "-t", help=\
        "The number of threads to use for calculations.",
        type=int,
        required=False,
        default=4)
    parser.add_argument("--split_time", "-s", help=\
        "The number of generations since the split with the outgroup (i.e. the \
entire height of the tree)", \
        type=int,
        required=True)
    return parser.parse_args()

def phi_diffusion(x, t, N, p):
    """
    Uses the diffusion approximation to the genetic drift equation to 
    calculate the probability of all possible allele frequencies in a population
    of size N.
    
    Reference for full solution including loss and fixation:
    http://www.sciencedirect.com/science/article/pii/S0022519307001981
    Reference for original solution (excluding loss & fixation) by Kimura, 1955:
    http://www.pnas.org/content/41/3/144.full.pdf
    
    Parameters:
        x: the allele frequency for which to calculate probability (eg 5/2N)
        t: the number of generations for which drift has been allowed to
            take place (eg 10)
        N: the population size
        p: the starting allele frequency (should be 1/2N)
    
    Returns:
        a tuple of loss_coeff, prob, fix_coeff
            loss_coeff is the weight of the Dirac delta function at allele
            frequency = 0 for this time point
            prob is the probability density at the given value of x
            fix_coeff is the weight of the Dirac delta function at allele
            frequency = 1 for this time point
    """
    tau = t/(N)
    
    # Stopping threshold. When a number being added to a sum is less than
    # this value, the infinite sum will be considered "converged."
    delta = 1e-20
    
    loss_tot = 0
    stop = False
    n = 0
    
    # http://docs.scipy.org/doc/scipy/reference/special.html
    
    # Coefficient of dirac delta function for allele loss
    while not stop:
        lambdacoeff = ((n+1)*(n+2))/2
        newsum = p * (2*n+3) * hyp2f1(-n, n+3, 2, p) *\
            math.exp(-lambdacoeff * tau)
        if newsum < delta:
            stop = True
        else:
            loss_tot += newsum
        n += 1
    loss_coeff = 1 - loss_tot
    loss_coeff *= (1-p)
    
    fix_tot = 0
    stop = False
    n = 0
    
    # Coefficient of dirac delta function for allele fixation
    while not stop:
        lambdacoeff = ((n+1)*(n+2))/2
        newsum = (1-p) * (2*n+3) * math.pow(-1, n+1) *\
            hyp2f1(-n, n+3, 2, p) * math.exp(-lambdacoeff * tau)
        if newsum < delta:
            stop = True
        else:
            fix_tot += newsum
        n += 1
    fix_coeff = 1 - fix_tot
    fix_coeff *= p
    
    # Main part of function (equivalent to Kimura's 1955 solution)
    
    tot = 0
    stop = False
    n = 0
    
    while not stop:
        lambdacoeff = ((n+1)*(n+2))/2
        newsum = (2*n+3) * (n+1) * (n+2) *\
            hyp2f1(-n, n+3, 2, p) * hyp2f1(-n, n+3, 2, x) *\
            math.exp(-lambdacoeff * tau)
        if newsum < delta:
            stop = True
        else:
            tot += newsum
        n += 1
    tot *= p*(1-p)
    
    return (loss_coeff, fix_coeff, tot)

def phi_diffusion_t(t, N, n, nbins=None):
    """
    Calculates probabilities of all possible allele frequencies at a given
    time point. Uses the above function (phi_diffusion) based on the diffusion
    approximation to the Wright-Fisher genetic drift problem. Since numbers
    can be weird in some cases, scales all results to ensure probabilities
    sum to 1.
    
    Parameters:
        t: the number of generations that genetic drift is allowed to take place
        N: the population size.
    
    Returns:
        a one-dimensional Numpy array of probabilities of various allele 
            frequencies, ranging from 0 (loss) to 2*N (fixation).
    """
    
    # We want one bin for every possible frequency (clade size), without accounting for
    # loss and fixation.
    nbins = n-1
    
    # Determine width of bins in terms of number of chromosomes in the population,
    # so we can accurately integrate the PDF.
    binwidth = 1/nbins
    
    # Define data structure to store results
    freq_vector = numpy.array([0] * nbins, dtype=numpy.float64)
    
    # These variables will store the weights of the Dirac delta functions
    # corresponding to allele loss and fixation.
    loss_coeff = 0
    fix_coeff = 0
    
    # Do calculation for every possible (discrete) allele frequency not 
    # corresponding to loss or fixation
    for bin_index in xrange(0, nbins):
        # The above function (Diffusion approximation) is a probability 
        # density function, so its results can only be interpreted across
        # ranges of values. Therefore, we will evaluate the PDF 0.5/(2*N) units
        # to the left and to the right of the desired frequency, and then
        # simulate integration by calculating the area under the curve via
        # the trapezoid rule. 
        
        bin_lower = (bin_index+1)/n - (0.5/n)
        bin_upper = (bin_index+1)/n + (0.5/n)
        
        #bin_lower = bin_index * binwidth - 0.5 * binwidth
        #bin_upper = bin_index * binwidth + 0.5 * binwidth
        
        loss_coeff, fix_coeff, prob_dens_lower = \
            phi_diffusion(bin_lower, t, N, 1/N)
        loss_coeff, fix_coeff, prob_dens_upper = \
            phi_diffusion(bin_upper, t, N, 1/N)
            
        if prob_dens_lower < prob_dens_upper:
            smaller = prob_dens_lower
            larger = prob_dens_upper
        else:
            smaller = prob_dens_upper
            larger = prob_dens_lower
        
        # Use trapezoid rule to get area under the curve and hence probability.
        p = (smaller * binwidth) + 0.5 * binwidth * (larger-smaller)
        
        # Store result.
        freq_vector[bin_index] = p

    # Build a "legend" where each index corresponds to an index in freq_vector
    # and denotes a population allele frequency.
    freq_legend = [0] * nbins
    for bin_index in xrange(0, nbins):
        freq_legend[bin_index] = "{}".format(bin_index+1)
        
    return (freq_vector, freq_legend)
    
def main(args):
    """Main method"""
    
    options = parse_args()
    
    pool = Pool(options.threads)
    threads = []
    
    #for start in xrange(0, 1, stepsize):
    #breaks = numpy.logspace(-options.gen_bins,1,options.gen_bins, base=10, endpoint=False)
    breaks = numpy.logspace(numpy.log10(2), numpy.log10(options.split_time*2), options.gen_bins, base=10, endpoint=False)
    
    for ind in xrange(0, len(breaks)):
        start = (breaks[ind])
        if ind < len(breaks)-1:
            binmid = ((breaks[ind+1])-start)/2 + start
        else:
            binmid = (1-start)/2 + start
            
        # This value is a % of the divergence to outgroup represented by this
        # bin.
        #gens = int(round(binmid*options.split_time*2))
        gens = numpy.round(binmid)
        thread = pool.apply_async(phi_diffusion_t, [gens, options.pop_size, options.sample_size])
        threads.append(thread)
        #vec, legend = phi_diffusion_t(gens, options.pop_size, nbins=options.bins)
    
    pool.close()
    pool.join()
    
    header_printed = False
    
    for thread_index, thread in enumerate(threads):
        vec, legend = thread.get()
        vec = vec/numpy.sum(vec)
        if not header_printed:
            for ind, elt in enumerate(legend):
                legend[ind] = str(elt)
            print("tmrca_bin\t{}".format("\t".join(legend)))
            header_printed = True
        vec2 = []
        for ind, elt in enumerate(vec):
            vec2.append(str(elt))
        endpt = 1
        if thread_index < len(breaks)-1:
            endpt = breaks[thread_index+1]
        print("{} {}\t{}".format(breaks[thread_index]/options.split_time, \
            endpt/options.split_time, "\t".join(vec2)))

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
