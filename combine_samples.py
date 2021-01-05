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
from ticktockclock import ticktock
import glob

import argparse

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

outsampleslist = glob.glob(os.path.join(outsamplesdir, 'sample_*.pkl'))

S = list()
results = dict()
for outsamples in outsampleslist:
    s = int(outsamples.split('.pkl')[0].split('_')[-1])
    try:
        with open(outsamples, 'rb') as f:
            res = joblib.load(f)

        results[s] = res

        S.append(s)
        print(s)
    except EOFError:
        print('sample_{}.pkl is not a correctly formatted pickle file'.format(s))


S.sort()

n = len(beta)

logZs = np.array([results[s].logz[-1] for s in S])
Nsamples = np.array([results[s].niter for s in S])

prob_s = softmax(logZs)

print('\nS:P(S)')
for i, s in enumerate(S):
    print('{}:{:.3e}'.format(s, prob_s[i]))

df = pd.DataFrame({'S':S, 'prob':prob_s})
df.S.astype(int)
df.to_csv(os.path.join(outfinaldir, "prob_of_S.csv"), index=False)

sns.set_style('white') 
sns.set_context("paper", font_scale=2.0)

fig, ax = plt.subplots()
sns.barplot(x=S, y=prob_s, color=sns.xkcd_rgb["denim blue"])
sns.despine()
plt.xlabel("Stem Cell Number")
plt.ylabel("Probability")
for label in ax.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(outfinaldir, "probability_S.png"), dpi = 300)
plt.close()


Ndraws = 3000
Ssamples = np.random.choice(S, size=Ndraws, p=prob_s)


final_posterior = np.empty((Ndraws, 8))
final_posterior[:, -1] = Ssamples
beta_hat = np.empty((1, Ndraws, n))
LL = np.empty((1, Ndraws, n))

progress_ints = (np.arange(0.1, 1.1, 0.1) * Ndraws - 1).astype(int)  
counter = 10
for i in range(Ndraws):

    if i in progress_ints:
        print('{}% complete'.format(counter))
        counter += 10

    s = Ssamples[i]
    posterior =  dynesty.utils.resample_equal(results[s].samples, softmax(results[s].logwt))
    random_row = np.random.randint(posterior.shape[0])
    final_posterior[i, :7] = posterior[random_row, :7]
    lamsample, musample, gammasample, deltasample, etasample = final_posterior[i, :5]

    kappasample = posterior[random_row, 7:]

    LL[0, i, :] = ticktock.loglikelihood_perpoint(posterior[random_row, :], beta, s, age)

    ProbDist = ticktock.runModel(s, lamsample, musample, gammasample, age)

    k_sample = np.random.choice(np.arange(0, 2*s+1), size=n, p=ProbDist)
    beta_sample = k_sample / (2*s)

    beta_sample = ticktock.rescale_beta(beta_sample, deltasample, etasample)

    beta_hat[0, i, :] = ticktock.beta_rvs(beta_sample, kappasample[k_sample])
    

with open(os.path.join(outfinaldir, "finalposterior.pkl"), 'wb') as f:
    joblib.dump(final_posterior, f)

df = pd.DataFrame({'lam':final_posterior[:,0],  
                    'mu':final_posterior[:,1], 
                    'gamma':final_posterior[:,2], 
                    'delta':final_posterior[:,3],
                    'eta':final_posterior[:,4], 
                    'kappamean':final_posterior[:,5],
                    'kappadisp':final_posterior[:,6],
                    'S':Ssamples})
df.to_csv(os.path.join(outfinaldir, "finalposterior.csv"), index=False)

fig, ax = plt.subplots()      
plt.hist(beta, np.linspace(0, 1, 100), density=True, alpha=0.4, linewidth=0) 
plt.hist(np.ravel(beta_hat), np.linspace(0, 1, 100), density=True, alpha=0.4, linewidth=0) 
plt.legend(("Data", "Posterior predictive"))
plt.xlabel("Fraction methylated")
plt.ylabel("Probability density")
sns.despine()
plt.tight_layout()
plt.savefig("{}/posterior_predictive.png".format(outfinaldir), dpi = 300)
plt.close()

inference = az.from_dict(posterior={'lam':final_posterior[:,0],  
                                    'mu':final_posterior[:,1], 
                                    'gamma':final_posterior[:,2], 
                                    'delta':final_posterior[:,3],
                                    'eta':final_posterior[:,4], 
                                    'kappamean':final_posterior[:,5],
                                    'kappadisp':final_posterior[:,6],
                                    'S':Ssamples}, 
                        observed_data={'beta':beta},
                        posterior_predictive={'beta_hat':beta_hat}, 
                        sample_stats={"log_likelihood":LL}
                                ) 

az.to_netcdf(inference, "{}/inference.nc".format(outfinaldir))

pairs = az.plot_pair(inference, var_names=('lam', 'mu', 'gamma', 'delta', 'eta', 'kappamean', 'kappadisp')) 
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