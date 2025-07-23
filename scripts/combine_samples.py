#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.special import logsumexp, softmax, gammaln, logit, expit
from scipy import linalg   
import joblib
from time import time
import arviz as az
import dynesty
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from flipflop import flipflop
import glob

import argparse

def main():
    parser = argparse.ArgumentParser(description='Combine the individual posteriors for each S value.')
    parser.add_argument('datafile', type=str,
                        help='path to csv containing beta values')
    parser.add_argument('patientinfofile', type=str,
                        help='path to csv containing patientinfo')
    parser.add_argument('outputdir', type=str, default='~', 
                        help='path to folder in which to store output')
    parser.add_argument('sample', type=str,
                        help='samplename of beta array (must be a col in datafile index in patientinfo)')

    # Execute the parse_args() method
    args = parser.parse_args()

    datafile = args.datafile
    patientinfofile = args.patientinfofile
    outputdir = args.outputdir
    sample = args.sample

    outsamplesdir = os.path.join(outputdir, sample, 'posterior')

    outfinaldir = os.path.join(outputdir, sample, 'outfinal')
    os.makedirs(outfinaldir, exist_ok=True)

    beta_values = pd.read_csv(datafile, index_col = 0)
    patientinfo = pd.read_csv(patientinfofile, keep_default_na=False, index_col = 0) 

    beta = beta_values[sample].dropna().values
    age = patientinfo.loc[sample, 'age']

    n = len(beta)

    outsamples = os.path.join(outsamplesdir, f'{sample}_sample.pkl')

    with open(outsamples, 'rb') as f:
        res = joblib.load(f)

    posterior =  dynesty.utils.resample_equal(res.samples, softmax(res.logwt))

    df = pd.DataFrame(posterior, columns = ['lam', 'mu', 'gamma', 'nu', 'zeta',
                                            'delta', 'eta', 'kappa', 'S'])
    
    df.to_csv(os.path.join(outfinaldir, "posterior.csv"), index=False)

    df_S = df['S'].astype(int).value_counts()
    df_S.to_csv(os.path.join(outfinaldir, "prob_of_S.csv"), index=False)

    sns.set_style('white') 
    sns.set_context("paper", font_scale=1.6)

    fig, ax = plt.subplots()
    ax.bar(df_S.index, df_S.values, color=sns.xkcd_rgb["denim blue"])
    sns.despine()
    plt.xlabel("Stem Cell Number (S)")
    plt.ylabel("Probability")
    plt.tight_layout()
    plt.savefig(os.path.join(outfinaldir, "probability_S.png"), dpi = 300)
    plt.close()


    Ndraws = 3000

    if posterior.shape[0] > Ndraws:
        random_rows = np.random.choice(posterior.shape[0], 
                        replace = False, size = Ndraws)
        final_posterior = posterior[random_rows, :]
    else: 
        final_posterior = posterior
        Ndraws = posterior.shape[0]

    beta_hat = np.empty((1, Ndraws, n))
    LL = np.empty((1, Ndraws, n))

    progress_ints = (np.arange(0.1, 1.1, 0.1) * Ndraws - 1).astype(int)  
    counter = 10
    for i in range(final_posterior.shape[0]):
        (lamsample, musample, gammasample, nusample, zetasample, 
            deltasample, etasample, kappasample, s) = final_posterior[i, :]

        LL[0, i, :] = flipflop.loglikelihood_perpoint(final_posterior[i, :], beta, age)

        ProbDist = flipflop.runModel(s, lamsample, musample, gammasample,
                                        nusample, zetasample, age)

        k_sample = np.random.choice(np.arange(0, 2*s+1), size=n, p=ProbDist)
        beta_sample = k_sample / (2*s)

        beta_sample = flipflop.rescale_beta(beta_sample, deltasample, etasample)

        beta_hat[0, i, :] = flipflop.beta_rvs(beta_sample, kappasample)
        
    fig, ax = plt.subplots()      
    plt.hist(beta, np.linspace(0, 1, 100), density=True, alpha=0.4, linewidth=0) 
    plt.hist(np.ravel(beta_hat), np.linspace(0, 1, 100), density=True, alpha=0.4, linewidth=0) 
    plt.legend(("Data", "Posterior predictive"))
    plt.xlabel("Fraction Methylated (Beta)")
    plt.ylabel("Probability density")
    sns.despine()
    plt.tight_layout()
    plt.savefig("{}/posterior_predictive.png".format(outfinaldir), dpi = 300)
    plt.close()

    inference = az.from_dict(posterior={col:final_posterior[:, j] 
                                        for j, col in enumerate(df.columns)},
                            observed_data={'beta':beta},
                            posterior_predictive={'beta_hat':beta_hat}, 
                            sample_stats={"log_likelihood":LL}
                                    ) 

    az.to_netcdf(inference, "{}/inference.nc".format(outfinaldir))

    pairs = az.plot_pair(inference, var_names=('lam', 'mu', 'gamma', 
                                            'nu', 'zeta', 'delta', 'eta',
                                            'kappa')) 
    plt.savefig('{}/plot_pairs.png'.format(outfinaldir), dpi=300)
    plt.close()

    az.plot_loo_pit(inference, y='beta', y_hat='beta_hat', ecdf=True)
    plt.savefig('{}/plot_loo_pit_ecdf.png'.format(outfinaldir), dpi=300)
    plt.close()

    sns.set_context("paper", font_scale=1.0)
    az.plot_loo_pit(inference, y='beta', y_hat='beta_hat')
    plt.ylabel('Leave One Out - Probability Integral Transform')
    plt.xlabel('Cumulative Density Function')
    plt.savefig('{}/plot_loo_pit.png'.format(outfinaldir), dpi=300)
    plt.close()

if __name__ == "__main__":
    main()