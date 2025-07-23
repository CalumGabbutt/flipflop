import os
import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.special import logsumexp, softmax, gammaln, logit, expit
from scipy import linalg   
import joblib
from time import time
import flipflop_model
import dynesty
from dynesty import NestedSampler

def runModel(S, lam, mu, gamma, nu, zeta, age):
    return flipflop_model.runModel(S, lam, mu, gamma, nu, zeta, age)

def beta_lpdf(y, alpha, beta):
    return flipflop_model.beta_lpdf(y, alpha, beta)

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

def lognormal_convert_params(mu, sigma):
    """
    Convert mean/std parameterization of a beta distribution to the ones 
    scipy supports
    """

    if np.any(mu <= 0):
        raise Exception("mu must be greater than 0")
    elif np.any(sigma <= 0) :
        raise Exception("sigma must be greater than 0")
    
    log_mean = np.log(mu) - 0.5 * np.log1p(sigma**2 / mu**2)
    log_std = np.sqrt(np.log1p(sigma**2 / mu ** 2))

    return log_mean, log_std

def lognormal_ppf(x, mu, sigma, **kwargs):
    """
    Calculate the ppf of a lognormal distribution with mean/std specified
    """
    log_mean, log_std = lognormal_convert_params(mu, sigma)
    return stats.lognorm.ppf(x, log_std, scale=np.exp(log_mean), **kwargs)


def prior_transform(flat_prior, scales, Smin, Smax):
    # priors for parameters [lam, mu, gamma, nu, zeta, delta, eta, kappa, S]
    prior = np.empty(np.shape(flat_prior))
    prior[:5] = stats.halfnorm.ppf(flat_prior[:5], scale=scales)

    # priors on delta and eta
    prior[5] = stats.beta.ppf(flat_prior[5], 5, 95)
    prior[6] = stats.beta.ppf(flat_prior[6], 95, 5)

    # mean and dispersion hyperpriors of kappa 
    prior[7] = lognormal_ppf(flat_prior[7], 100, 30)

    prior[8] = stats.randint.ppf(flat_prior[8], Smin, Smax + 1)

    return prior

def loglikelihood_perpoint(params, y, age):

    # unpack parameters
    lam, mu, gamma, nu, zeta, delta, eta, kappa, S = params
    S = int(S)

    # calucalate p(z|lam, mu, gamma)
    ProbDist = runModel(S, lam, mu, gamma, nu, zeta, age)
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


def loglikelihood(params, beta, age):
    return np.sum(loglikelihood_perpoint(params, beta, age))


def run_inference(
    beta, 
    age, 
    nlive=1500, 
    lamscale=1.0, 
    muscale=0.05, 
    gammascale=0.05, 
    nuscale=0.05, 
    zetascale=0.05, 
    Smin=3,
    Smax=20,
    verbose=False
):

    # set the std of the halfnormal priors on lam, mu, gamma
    scales = np.array([lamscale, muscale, gammascale, nuscale, zetascale])
    ndims = 9
    
    # create dummy functions which takes only the params as an argument to pass
    # to dynesty
    prior_function = lambda flat_prior: prior_transform(flat_prior, scales, Smin, Smax)
    loglikelihood_function = lambda params: loglikelihood(params, beta, age)

    t0 = time()

    print('Performing Dynesty sampling')
    sampler = NestedSampler(loglikelihood_function, prior_function, ndims,
                            bound='multi', sample='rwalk', nlive=nlive)
    sampler.run_nested(print_progress=verbose)
    res = sampler.results

    t1 = time()
    timesampler = int(t1-t0)

    print("\nTime taken to run Dynesty is {} seconds".format(timesampler))

    print(res.summary())

    return res
