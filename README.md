# ticktockclock
Python code to accompany our paper "Reconstructing Contemporary Human Stem Cell Dynamics with Oscillatory Molecular Clocks"

This package is designed for the inference of stem cell dynamics from methylation array data (e.g. Illumina EPIC arrays) of glandular tissue. 

# Installation 
`ticktockclock` is compatible with Python 3.6. It requires `numpy` (basic maths functions), `scipy` (special maths functions), `pandas` (dataframes manipulation), `dynesty` (nested sampling - Bayesian inference), `joblib` (pickling large data files, preserving dynesty structure),`matplotlib` (plotting), `seaborn` (plotting), `arviz` (Bayesian plotting), `cython` (accelerated computation)

The package can be installed directly from a local copy of the Github repo by running:

`python3 setup.py install`

# Usage
The functions used to perform the Bayesian inference can be imported in python using the command
`from ticktockclock import ticktock`

A number of key functions, in particular calculating key components of the log-likelihood function which is used as an input to dynesty, are computationally expensive and so have been written in Cython.

An example script `run_inference.sh` has been included to run the Bayesian inference framework described in the manuscript on one of the methylation arrays included in the paper, reproducing some of the plots included in the supplementary information of the manuscript. However, due to the computational expense of the inference (the analysis in the paper used 1500 live-points in the nested sampling step and consisted of 4 independent chains, initialised randomly), the analysis was in fact performed in parallel on a HPC. 

To run the example script, run the following code snippet from the command line, whilst located in the main directory copied from the Gitgub repo:
`./run_inference.sh`

(if a permission denied error is returned, run the command `chmod +x run_inference.sh` to grant access to the script)

To adapt the script to be applied to different datasets, simply change the variable definitions in the `./run_inference.sh` script:
-The datafile variable must point to the path of a csv file containing the beta values of CpG loci that have been identified as oscillatory (a histogram of the ticktock CpGs beta values should have a w-shape, as discussed in the manuscript). 
-The patientinfofile variable must point to the path of a csv, indexed by the samplenames used as the column headers in the datafile and containing an `age` column. S
-The samplename variable must be a string corresponding to one of the column headers in the datafile csv and one of the indices of the patientinfofile csv.
-Smin and Smax determine the range of stem cell numbers to be considered by the inference pipeline

The pipeline consists of a sequential for loop which uses the dynesty nested sampling algorithm to calculate the posterior and Bayesian evidence for each stem cell number in a given range. The posterior probability for S can then be found by applying Bayes rule (see manuscript for details). 

Following the calculation of the composite posterior, a series of plots are generated to evaluate the quality of the fitting.
