#!/usr/bin/env python3

import pandas as pd
import numpy as np
from flipflop import flipflop
import argparse

def draw_beta(S, lam, mu, gamma, delta, eta, kappa, age, n):

    ProbDist = flipflop.runModel(S, lam, mu, gamma, age)

    # Randomly draw n beta values from the ProbDist distribution (discrete dist.)
    beta_sample = np.random.choice(np.linspace(0, 1, 2*S+1), size=n, p=ProbDist)

    # Linear transform on beta to account for array saturating (e.g. shifts mean of
    # lower peak from 0 to delta and upper peak from 1 to eta)
    beta_rescale = flipflop.rescale_beta(beta_sample, delta, eta)

    # Apply noise model (for each true beta value, draw the measured beta
    # value from a beta distribution with a mean equal to the true beta value
    # and a sample size equal to kappatrue)
    return flipflop.beta_rvs(beta_rescale, kappa)


def main():
    parser = argparse.ArgumentParser(description='Simulate synthetic crypt')
    parser.add_argument('outputfile', type=str, default='~', 
                        help='path to csv file in which to store output')
    parser.add_argument('-samplename', action='store', default='testsample', type=str,
                        help='samplename of the generated beta values (column header)')
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
    n = int(args.n)
    S = int(args.S)
    lam = args.lam
    mu = args.mu
    gamma = args.gamma
    delta = args.delta
    eta = args.eta
    kappa = args.kappa
    age = args.age

    beta_noisy = draw_beta(S, lam, mu, gamma, delta, eta, kappa, age, n)
    beta = pd.DataFrame({samplename:beta_noisy})

    beta.to_csv(outputfile, index=False)

if __name__ == "__main__":
    main()
