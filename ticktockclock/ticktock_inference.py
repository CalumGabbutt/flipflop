import os
import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.special import logsumexp, softmax, gammaln, logit, expit
from scipy import linalg   
import joblib
from time import time
from ticktock_model import runModel, beta_lpdf
import dynesty
from dynesty import NestedSampler

def rescale_beta(beta, delta, eta):
    # Linear transform of beta values from between 
    # 0 and 1 to between delta and eta
    return (eta - delta) * beta + delta

def beta_convert_params(mu, kappa):
    """
    Convert mean/dispersion parameterization of a beta distribution to the ones scipy supports

    """

    if np.any(kappa <= 0):
        raise Exception("kappa must be greater than 0")
    elif np.any(mu <= 0) or np.any(mu >= 1):
        raise Exception("mu must be between 0 and 1")
    
    alpha = kappa * mu 
    beta = kappa * (1- mu)

    return alpha, beta


def beta_ppf(y, mean, kappa):

    alpha, beta = beta_convert_params(mean, kappa)

    return stats.beta.ppf(y, alpha, beta)

def beta_rvs(mean, kappa, **kwargs):

    alpha, beta = beta_convert_params(mean, kappa)

    return stats.beta.rvs(alpha, beta, **kwargs)


def truncnormal_convert_params(mean, std, lb=0.0, ub=1.0):
    """
    Convert mean/dispersion parameterization of a beta distribution to the ones scipy supports

    """

    if np.isfinite(lb):
        a = (lb - mean) / std
    else:
        a = lb
    
    if np.isfinite(ub):
        b = (ub - mean) / std
    else: 
        b = ub

    return a, b


def truncnormal_ppf(y, mean, std, lb=0.0, ub=1.0):

    a, b = truncnormal_convert_params(mean, std, lb, ub)

    return stats.truncnorm.ppf(y, a, b, loc=mean, scale=std)

    
def prior_transform(flat_prior, scales):
    # priors for parameters [lam_S, mu, gamma, delta, eta, kappamean, kappadisp, kappa[]]
    prior = np.empty(np.shape(flat_prior))
    prior[:3] = stats.halfnorm.ppf(flat_prior[:3], scale=scales)

    # priors on delta and eta
    prior[3] = stats.beta.ppf(flat_prior[3], 5, 95)
    prior[4] = stats.beta.ppf(flat_prior[4], 95, 5)

    # mean and dispersion hyperpriors of kappa 
    prior[5] = stats.halfnorm.ppf(flat_prior[5], scale = 500)
    prior[6] = stats.halfnorm.ppf(flat_prior[6], scale = 50)

    # hierarchical priors on kappa[] 
    prior[7:] = truncnormal_ppf(flat_prior[7:], prior[5], prior[6], lb=0.0, ub=np.inf)
    return prior

def loglikelihood_perpoint(params, y, S, age):

    # unpack parameters
    lam, mu, gamma, delta, eta, kappamean, kappadisp = params[:7]
    kappa = params[7:]

    # calucalate p(z|lam, mu, gamma)
    ProbDist = runModel(S, lam, mu, gamma, age)
    lProbDist = np.log(ProbDist)

    # repeat the above values to create a (2S+1, n) array
    lpk = np.repeat(lProbDist[:, np.newaxis], np.shape(y)[0], axis=1)

    # linear transform of "true" beta peaks
    betatrue = rescale_beta(np.linspace(0, 1, 2*S+1), delta, eta)

    # transform the mean and shape parameters to the ones scipy supports
    alpha, beta = beta_convert_params(betatrue, kappa)

    # calculate log(p(y|z))
    lpk += beta_lpdf(y, alpha, beta)

    # p(y|lam, mu, gamma) = sum(p(y|z)*p(z|lam, mu, gamma))
    # on the log scale this requires logsumexp 
    LL = logsumexp(lpk, axis = 0)

    return LL


def loglikelihood(params, beta, S, age):
    return np.sum(loglikelihood_perpoint(params, beta, S, age))


def run_inference(beta, age, S, nlive=1500, lamscale=1.0, muscale=0.05, gammascale=0.05, print_progress=True):
    # set the std of the halfnormal priors on lam, mu, gamma
    scales = np.array([lamscale, muscale, gammascale])
    ndims = 7 + 2*S +1
    
    # create dummy functions which takes only the params as an argument to pass
    # to dynesty
    prior_function = lambda flat_prior: prior_transform(flat_prior, scales)
    loglikelihood_function = lambda params: loglikelihood(params, beta, S, age)

    t0 = time()
    mode = 'rwalk'
    print('Performing {} sampling'.format(mode))
    sampler = NestedSampler(loglikelihood_function, prior_function, ndims,
                            bound='multi', sample=mode, nlive=nlive)
    if print_progress=True:
        sampler.run_nested(print_progress=True)
    else:
        sampler.run_nested(print_progress=False)
    res = sampler.results
    t1 = time()

    timenestle = int(t1-t0)
    print("\nTime taken to run 'Dynesty' is {} seconds".format(timenestle))

    print(res.summary())

    return res