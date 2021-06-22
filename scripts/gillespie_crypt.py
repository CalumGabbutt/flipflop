#!/usr/bin/env python3

import numpy as np
import pandas as pd
import random
from ticktockclock import ticktock
import os
import sys
import copy

import argparse

def choose_cell_ring(S):
    # Function to choose 2 adjacent cells at random, assuming a ring geometry 
    # (i.e. the first and last cells in the array are neighbours)

    # Choose one cell at random
    cell = np.random.randint(S)

    # We must account for the periodic boundary conditions
    if cell == 0:
        if random.random() < 0.5:
            cell2 = S - 1
        else:
            cell2 = 1
    elif cell == S - 1:
        if random.random() < 0.5:
            cell2 = S - 2
        else:
            cell2 = 0
    else:
        if random.random() < 0.5:
            cell2 = cell - 1
        else:
            cell2 = cell + 1

    return cell, cell2

def choose_cell_mixed(S):
    # Function to choose to 2 cells at random, assuming a well-mixed geometry 
    # (i.e. sampling without replacement)
    return np.random.choice(np.arange(S), size=2, replace=False)

def methylationSim(allele, lam, mu, gamma, S, age, mode='mixed'):

    if mode not in ['mixed', 'ring']:
        raise ValueError("mode must be either 'ring' or 'mixed'")

    # Create a 3D array representing the stemcell niche,
    # where the first index determines the cell, the second index 
    # determines the DNA strand and the third index determines the 
    # CpG site
    dna = [copy.deepcopy(allele) for i in range(2)]
    dna = np.stack(dna)
    stemcells = np.stack([copy.deepcopy(dna) for i in range(S)])

    # Use a Gillepsie algorithm to model individual CpG changing methylation 
    # status and the stem cells replacing each other
    n = np.shape(dna)[1]

    t = 0
    counter = 0
    while t < age:
        # Draw the time until the next replacement event, assuming the 
        # replacements are Poisson distributed
        twait = 1 / (lam * S) * np.log(1 / random.random())

        if t + twait < age:

            t += twait 

            # Count the number of methylated/non-methylated sites
            num_on = np.count_nonzero(stemcells)
            num_off = 2*n*S - num_on

            # Create a copy of the stemcells array to act as a mask
            boolean_mask = copy.deepcopy(stemcells)

            # Calculate the probability that a site is methylated at time (t+twait),
            # given that the site was methylated at time t
            p_given_on = mu / (gamma+mu) + (1 - mu/(gamma+mu))*np.exp(-(mu+gamma)*twait)
            # For each site methylated at time t, draw from a random distribution 
            # with the above probability to see whether it is still methylated
            stemcells[boolean_mask] = np.random.choice([True, False], size=num_on, p=[p_given_on, 1-p_given_on])

            # Calculate the probability that a site is methylated at time (t+twait),
            # given that the site was not methylated at time t
            p_given_off = mu / (gamma+mu) * (1 - np.exp(-(mu+gamma)*twait))
            # For each site unmethylated at time t, draw from a random distribution 
            # with the above probability to see whether it is methylated
            stemcells[~boolean_mask] = np.random.choice([True, False], size=num_off, p=[p_given_off, 1-p_given_off])

            # Now that we have taken care of any mutations that have occured, 
            # pick 2 cells at random, one to replace the other
            if mode == 'mixed':
                expanding_cell, retracting_cell = choose_cell_mixed(S)
            elif mode == 'ring':
                expanding_cell, retracting_cell = choose_cell_ring(S)
            
            stemcells[retracting_cell] = copy.deepcopy(stemcells[expanding_cell])

            counter += 1

        else:
            twait = age - t
            t += twait

            # Count the number of methylated/non-methylated sites
            num_on = np.count_nonzero(stemcells)
            num_off = 2*n*S - num_on

            # Create a copy of the stemcells array to act as a mask
            boolean_mask = copy.deepcopy(stemcells)

            # Calculate the probability that a site is methylated at time (t+twait),
            # given that the site was methylated at time t
            p_given_on = mu / (gamma+mu) + (1 - mu/(gamma+mu))*np.exp(-(mu+gamma)*twait)
            # For each site methylated at time t, draw from a random distribution 
            # with the above probability to see whether it is still methylated
            stemcells[boolean_mask] = np.random.choice([True, False], size=num_on, p=[p_given_on, 1-p_given_on])

            # Calculate the probability that a site is methylated at time (t+twait),
            # given that the site was not methylated at time t
            p_given_off = mu / (gamma+mu) * (1 - np.exp(-(mu+gamma)*twait))
            # For each site unmethylated at time t, draw from a random distribution 
            # with the above probability to see whether it is methylated
            stemcells[~boolean_mask] = np.random.choice([True, False], size=num_off, p=[p_given_off, 1-p_given_off])

            break

    return stemcells

def noise(beta_true, delta, eta, kappa):

    # Linear transform on beta to account for array saturating (e.g. shifts mean of
    # lower peak from 0 to delta and upper peak from 1 to eta)
    beta_rescale = ticktock.rescale_beta(beta_true, delta, eta)

    # Apply noise model (for each true beta value, draw the measured beta
    # value from a beta distribution with a mean equal to the true beta value
    # and a sample size equal to kappa)
    beta_sample = ticktock.beta_rvs(beta_rescale, kappa)

    return beta_sample

def generateAllele(n):
    # Generate a random methylation pattern           
    allele = np.random.choice([True, False], size=n)

    allele = np.full(n, False)
    idx = np.arange(n)
    np.random.shuffle(idx)
    allele[idx[:n//2]] = True
    return allele

def runSimsCounts(lam, mu, gamma, S, n, N, age, mode='mixed'):

    unique_values = np.linspace(0, 1, 2*S + 1)
    counts = np.zeros((N, 2*S+1), dtype=int)
    for i in range(N):
        allele = generateAllele(n)
        stemcells = methylationSim(allele, lam, mu, gamma, S, age, mode=mode)

        mean_beta = stemcells.mean(axis=(0, 1))
        counts[i, :] = [(abs(mean_beta - uniq) < 0.001).sum() for uniq in unique_values]

    return counts

def main():
    parser = argparse.ArgumentParser(description='Simulate synthetic crypt')
    parser.add_argument('outputfile', type=str, default='~', 
                        help='path to csv file in which to store output')
    parser.add_argument('mode', type=str, default='mixed', 
                        help='stem cell geometry (must be "ring" or "mixed" - default "mixed")')
    parser.add_argument('-samplename', action='store', default='testsample', type=str,
                        help='samplename of the generated beta values (column header)')
    parser.add_argument('-N', default=100, type=int, help='number of crypts to average over')                  
    parser.add_argument('-n', default=2000, type=int, help='number of CpG loci')
    parser.add_argument('-S', default=5, type=int, help='number of stem cells')
    parser.add_argument('-lam', default=1.0, type=float, 
                        help='replacement rate per stem cell')
    parser.add_argument('-mu', default=0.05, type=float, 
                        help='methylation rate per CpG site per year')      
    parser.add_argument('-gamma', default=0.05, type=float, 
                        help='demethylation rate per CpG site per year')  
    parser.add_argument('-delta', default=0.04, type=float, 
                        help='offset from zero') 
    parser.add_argument('-eta', default=0.92, type=float, 
                        help='offset from one')   
    parser.add_argument('-kappa', default=100, type=float, 
                        help='''
                        sample size of beta distribution (see 
                        https://en.wikipedia.org/wiki/Beta_distribution#Mean_and_sample_size)
                        ''')   
    parser.add_argument('-age', default=30, type=float, 
                        help='age of patient')  
    # Execute the parse_args() method
    args = parser.parse_args()

    outputfile = args.outputfile
    samplename = args.samplename
    mode = args.mode
    N = int(args.N)
    n = int(args.n)
    S = int(args.S)
    lam = args.lam
    mu = args.mu
    gamma = args.gamma
    delta = args.delta
    eta = args.eta
    kappa = args.kappa
    age = args.age

    if mode not in ['mixed', 'ring']:
        raise ValueError("mode must be either 'ring' or 'mixed'")
    if (delta >= 1) | (delta <= 0):
        raise Exception('delta must be between 0 and 1')
    if (eta >= 1) | (eta <= 0):
        raise Exception('eta must be between 0 and 1')
    if (kappa <= 0):
        raise Exception('kappa must be greater than 0')

    counts = runSimsCounts(lam, mu, gamma, S, n, N, age, mode=mode)
    mean = np.round(counts.mean(axis=0)).astype(int)
    beta_arr = np.repeat(np.linspace(0, 1, 2*S+1), mean)
    beta_noisy = noise(beta_arr, delta, eta, kappa)

    beta = pd.DataFrame({samplename:beta_noisy})

    beta.to_csv(outputfile)

if __name__ == "__main__":
    main()