# ticktockclock
Python code to accompany our paper "Reconstructing Contemporary Human Stem Cell Dynamics with Oscillatory Molecular Clocks"

This package is designed for the inference of stem cell dynamics from methylation array data (e.g. Illumina EPIC arrays) of glandular tissue. 

# Installation 
`ticktockclock` is compatible with Python 3.6. It requires `numpy` (basic maths functions), `scipy` (special maths functions), `pandas` (analysis of dataframes), `dynesty` (nested sampling - Bayesian inference), `joblib` (pickling large data files, preserving dynesty structure),`matplotlib` (plotting), `seaborn` (plotting), `arviz` (Bayesian plotting)

The package can be installed directly from a local copy of the Github repo by running:

`python3 setup.py install`

