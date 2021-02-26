import pandas as pd
from ticktockclock import ticktock
import os
import sys
import joblib

import argparse

parser = argparse.ArgumentParser(description='Run ticktock Bayesian inference.')
parser.add_argument('datafile', type=str,
                    help='path to csv containing beta values')
parser.add_argument('patientinfofile', type=str,
                    help='path to csv containing patientinfo')
parser.add_argument('outputdir', type=str, default='~', 
                    help='path to folder in which to store output')
parser.add_argument('sample', type=str,
                    help='samplename of beta array (must be a col in datafile index in patientinfo)')
parser.add_argument('S', type=int,
                    help='stem cell number to evaluate')
parser.add_argument('-nlive', default=1500, dest='nlive', type=int,
                    help='number of live points in dynesty sampler (default:1500)')
parser.add_argument('--verbose', action='store_true', default=False, dest='verbose')
parser.add_argument('-mode', default='dynesty', type=str,
                    help='which nested sampling tool to use ["dynesty", "ultranest"] (default "dynesty")')

# Execute the parse_args() method
args = parser.parse_args()

datafile = args.datafile
patientinfofile = args.patientinfofile
outputdir = args.outputdir
sample = args.sample
S = args.S
nlive = args.nlive
verbose = args.verbose
mode = args.mode

if mode not in ["dynesty", "ultranest"]:
    raise Exception('mode argument must be in ["dynesty", "ultranest"]')

outsamplesdir = os.path.join(outputdir, sample, 'posterior')
outsamples = os.path.join(outsamplesdir, 'sample_{}.pkl'.format(S))

os.makedirs(outsamplesdir, exist_ok=True)

beta_values = pd.read_csv(datafile, index_col = 0)
patientinfo = pd.read_csv(patientinfofile, keep_default_na=False, index_col = 0) 

beta = beta_values[sample].dropna().values
age = patientinfo.loc[sample, 'age']

res = ticktock.run_inference(beta, age, S, nlive=nlive, 
                            verbose=verbose, mode=mode)

with open(outsamples, 'wb') as f:
    joblib.dump(res, f)
