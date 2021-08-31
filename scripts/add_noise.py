#!/usr/bin/env python3

import pandas as pd
import numpy as np
from flipflop import flipflop
import argparse

def noise(beta_true, delta, eta, kappa):

    # Linear transform on beta to account for array saturating (e.g. shifts mean of
    # lower peak from 0 to delta and upper peak from 1 to eta)
    beta_rescale = flipflop.rescale_beta(beta_true, delta, eta)

    # Apply noise model (for each true beta value, draw the measured beta
    # value from a beta distribution with a mean equal to the true beta value
    # and a sample size equal to kappa)
    beta_sample = flipflop.beta_rvs(beta_rescale, kappa)

    return beta_sample

def main():

    # Initialise argparser
    parser = argparse.ArgumentParser(description='Apply noise model to a csv of "true" beta values')
    parser.add_argument('inputfile', type=str,
                        help='path to csv file containing "true" beta values')
    parser.add_argument('outputfile', type=str, 
                        help='path to csv file in which to store output')
    parser.add_argument('-delta', default=0.04, type=float, 
                        help='offset from zero (default 0.04)') 
    parser.add_argument('-eta', default=0.92, type=float, 
                        help='offset from one (default 0.92)')   
    parser.add_argument('-kappa', default=100, type=float, 
                        help='sample size of beta distribution (default 100)\n(see https://en.wikipedia.org/wiki/Beta_distribution#Mean_and_sample_size)')   
    parser.add_argument('--index', action='store_true', default=False, dest='index',
                        help='indicate whether the first column of the inputfile is the index of the csv (default False)')

    # Execute the parse_args() method
    args = parser.parse_args()

    inputfile = args.inputfile
    outputfile = args.outputfile
    delta = float(args.delta)
    eta = float(args.eta)
    kappa = float(args.kappa)
    index = args.index

    if (delta >= 1) | (delta <= 0):
        raise Exception('delta must be between 0 and 1')
    if (eta >= 1) | (eta <= 0):
        raise Exception('eta must be between 0 and 1')
    if (kappa <= 0):
        raise Exception('kappa must be greater than 0')

    # Load csv file as dataframe
    if index:
        beta_true = pd.read_csv(inputfile, index_col=0)
    else:
        beta_true = pd.read_csv(inputfile)

    # Add noise to the true beta values
    beta_sample = noise(beta_true, delta, eta, kappa)

    # Save the noisy beta values
    if index:
        beta = pd.DataFrame(beta_sample, columns=beta_true.columns, index=beta_true.index)
        beta.to_csv(outputfile)
    else:
        beta = pd.DataFrame(beta_sample, columns=beta_true.columns)
        beta.to_csv(outputfile, index=False)

if __name__ == "__main__":
    main()