# flipflop
Python code to accompany our paper "Cell lineage tracing with molecular clocks based on fluctuating DNA methylation"

This package is designed for the inference of stem cell dynamics from methylation array data (e.g. Illumina EPIC arrays) of glandular tissue. 

# Installation 
`flipflop` is compatible with Python 3.6 or 3.7 It requires `numpy` (basic maths functions), `scipy` (special maths functions), `pandas` (dataframes manipulation), `dynesty` (nested sampling - Bayesian inference), `joblib` (pickling large data files, preserving dynesty structure),`matplotlib` (plotting), `seaborn` (plotting), `arviz` (Bayesian plotting), `cython` (accelerated computation)

The package can be installed directly from a local copy of the Github repo. We reccommend installing flipflop in a virtual environment, using venv, and pip to install the dependencies (conda can also be used):

```
git clone https://github.com/CalumGabbutt/flipflop.git
cd flipflop
python3 -m venv flipflopenv
source flipflopenv/bin/activate
pip install -r requirements.txt
python3 setup.py install
```

# Usage
The functions used to perform the Bayesian inference can be imported in python using the command
`from flipflop import flipflop`

A number of key functions, in particular calculating key components of the log-likelihood function which is used as an input to dynesty, are computationally expensive and so have been written in Cython.

An example script `run_inference.sh` has been included to run the Bayesian inference framework described in the manuscript on one of the methylation arrays included in the paper, reproducing some of the plots included in the supplementary information of the manuscript. However, due to the computational expense of the inference (the analysis in the paper used 1500 live-points in the nested sampling step and consisted of 4 independent chains, initialised randomly), the analysis was in fact performed in parallel on a HPC. 

To run the example script, run the following code snippet from the command line, whilst located in the main directory copied from the Gitgub repo:
`./run_inference.sh`

(if a permission denied error is returned, run the command `chmod +x run_inference.sh` to grant access to the script)

To adapt the script to be applied to different datasets, simply change the variable definitions in the `./run_inference.sh` script:
* The datafile variable must point to the path of a csv file containing the beta values of CpG loci that have been identified as oscillatory (a histogram of the flipflop CpGs beta values should have a w-shape, as discussed in the manuscript). 
* The patientinfofile variable must point to the path of a csv, indexed by the samplenames used as the column headers in the datafile and containing an `age` column.
* The samplename variable must be a string corresponding to one of the column headers in the datafile csv and one of the indices of the patientinfofile csv.
* Smin and Smax determine the range of stem cell numbers to be considered by the inference pipeline
* nlive (optional) is the number of live points used in the nested sampling step, more live points means less sampling related uncertainty, but takes longer

The pipeline consists of a sequential for loop which uses the dynesty nested sampling algorithm to calculate the posterior and Bayesian evidence for each stem cell number in a given range. The posterior probability for S can then be found by applying Bayes rule (see manuscript for details). 

Following the calculation of the composite posterior, a series of plots are generated to evaluate the quality of the fitting.

# Fluctuating CpG Selection
It is likely that fluctuating CpG (fCpG) loci are tissue-specific, hence the script `select_cpgs.py` has been included to identify a set of novel Tick-tock CpG loci from a new dataset (containing multiple single-gland samples for each patient). The code follows the same logic as in the manuscript, we first filter out Type I probes, CpG loci located on CpG Islands or Shores, or on non-autosomal chromosomes. Next we filter for probes with a low signal intensity (hence poor precision). For each patient we calculate the standard deviation at each CpG locus across all the samples for that patient, and select the set of CpG loci with the top `p`% (default 5%) average heterogeneity as putative fCpG loci. The average methylation fraction (beta value) across the entire patient is expected to follow a distribution centred approximately at 50%, however to enrich for fCpG loci where the methylation and demethylation rates are similar, we select only CpG loci with an average beta value between 0.4 and 0.6. 

To run this script, simply specify the name of the output csv file (`outputfile`), the path to the beta values, methylated and unmethylated intensities (`betafile`, `methylatedfile`, `unmethylatedfile`), and a patientinfo file (`patientinfo`) containing at a minimum each of the samples as the first column and a column labelled `patient`. Additionally, if only a subset of the samples are to be included, two additional optional arguments may be used - if only a specific tissue type is to be included, then the keyword `--tissue` can be set, this should refer to the labels in the `patientinfo` column `location`. Similarly, if only samples with a specific disease state are to be included, the keyword `--disease` can be set, which should refer to the labels in the `patientinfo` column `disease`. For further information, type `select_cpgs.py --help`.

Example usage:
`select_cpgs.py beta_flipflop_colon.csv beta_values-noob_all.csv Methylated_intensities-noob_all.csv Unmethylated_intensities-noob_all.csv patientinfo.csv -t Colon -d -`