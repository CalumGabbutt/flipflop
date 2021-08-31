#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import os
import flipflop
import inspect

def filter_noise(beta, M, U, r_min=500):
    fil_index = M+U >= r_min

    if len(np.shape(beta)) > 1:
        fil_index = fil_index.all(axis=1)

    M_fil = M[fil_index]
    U_fil = U[fil_index]
    beta_fil = beta[fil_index]

    return beta_fil, M_fil, U_fil

def subset_sites(beta, M, U, samples, patients, patientinfo, p=5):

    patientinfo = patientinfo.loc[samples]

    beta_fil, M_fil, U_fil = filter_noise(beta, M, U, r_min=1200) 

    mean_beta = np.mean(beta_fil, axis = 1)

    patient_std = pd.DataFrame({})
    for patient in patients:
        patient_samples = list(patientinfo.loc[patientinfo['patient']==patient, :].index)  
        patient_std[patient] = np.std(beta_fil[patient_samples], axis = 1)

    mean_std = np.mean(patient_std, axis = 1)
    
    thresh = np.percentile(mean_std, 100-p)
    index = (mean_std > thresh) & (0.4 < mean_beta) & (mean_beta < 0.6)

    return index[index == True].index 

def main():
    parser = argparse.ArgumentParser(description='Select tick-tock CpG Loci')
    parser.add_argument('outputfile', type=str,
                        help='path to csv file in which to store output')
    parser.add_argument('betafile', type=str, 
                        help='path to csv file containing noob beta values')
    parser.add_argument('methylatedfile', type=str,
                        help='path to csv file containing noob M values')
    parser.add_argument('unmethylatedfile', type=str,
                        help='path to csv file containing noob U values')
    parser.add_argument('patientinfofile', type=str,
                        help='''
                        path to csv file containing patient details (must 
                        contain sample IDs that are the columns of betafile
                        and a column titled "patient")
                        ''')    
    parser.add_argument('-e', '--manifestfile', default=None, 
                        help='path to Illumina EPIC manifest csv file (default MethylationEPIC_v-1-0_B4.csv)')
    parser.add_argument('-c', '--crossreactivefile', default=None, 
                        help='path to csv file containing CpGs to exclude (default 13059_2016_1066_MOESM1_ESM.csv)')
    parser.add_argument('-t', '--tissue', default=None, 
                        help='''
                        str determining which samples to include - if not None then 
                        patientinfo must contain a column titled "location" 
                        containing at least one entry corresponding to this option
                        (default None)
                        ''')
    parser.add_argument('-d', '--disease', default=None, 
                        help='''
                        str determining which samples to include - if not None then 
                        patientinfo must contain a column titled "disease" 
                        containing at least one entry corresponding to this option
                        (default None)
                        ''')
    parser.add_argument('-p', '--percent', default=5, 
                        help='''
                        float specifying the top p%% of most heterogeneous CpG loci
                        to include in the 
                        ''')
    # Execute the parse_args() method
    args = parser.parse_args()

    outputfile = args.outputfile
    betafile = args.betafile
    methylatedfile = args.methylatedfile
    unmethylatedfile = args.unmethylatedfile
    patientinfofile = args.patientinfofile
    manifestfile = args.manifestfile
    crossreactivefile = args.crossreactivefile
    tissue = args.tissue
    disease = args.disease
    p = float(args.percent)

    packagedir = os.path.dirname(inspect.getfile(flipflop))
    data_directory = os.path.join(packagedir, "files")

    if manifestfile is None:
        manifestfile = os.path.join(data_directory, 'MethylationEPIC_v-1-0_B4.csv')
    if crossreactivefile is None:
        crossreactivefile = os.path.join(data_directory, '13059_2016_1066_MOESM1_ESM.csv')

    if (p <= 0) | (p >= 100):
        raise Exception("p must be number between 0 and 100")

    beta_all = pd.read_csv(betafile, index_col = 0)
    M_all = pd.read_csv(methylatedfile, index_col = 0)
    U_all = pd.read_csv(unmethylatedfile, index_col = 0)

    patientinfo = pd.read_csv(patientinfofile, index_col = 0, keep_default_na = False) 
    crossreactive = pd.read_csv(crossreactivefile, index_col = 0)  
    manifest = pd.read_csv(manifestfile, index_col = 0, low_memory = False) 

    difference = beta_all.index.difference(crossreactive.index)

    beta_all = beta_all.loc[difference]
    M_all = M_all.loc[difference]
    U_all = U_all.loc[difference]

    manifest = manifest.loc[manifest.index.intersection(beta_all.index), :]

    good_sites = manifest.loc[(manifest['Relation_to_UCSC_CpG_Island'].isna() 
                & ~(manifest['CHR'].str.contains('X|Y')) 
                & (manifest['Infinium_Design_Type'] == 'II')), :].index

    if tissue is None and disease is None:
        samples = patientinfo.index
    elif disease is None:
        samples = patientinfo.loc[patientinfo['location']==tissue, :].index  
    elif tissue is None:
        samples = patientinfo.loc[patientinfo['disease']==disease, :].index  
    else:
        samples = patientinfo.loc[(patientinfo['disease']==disease) &
                                        (patientinfo['location']==tissue), :].index  

    patients = pd.unique(patientinfo.loc[samples, 'patient'])  

    beta_loc = beta_all.loc[good_sites, samples]
    M_loc = M_all.loc[good_sites, samples]
    U_loc = U_all.loc[good_sites, samples]

    ticktock_index = subset_sites(beta_loc, M_loc, U_loc, samples, patients, patientinfo, p)

    beta = beta_all.loc[ticktock_index, samples]

    beta.to_csv(outputfile)

if __name__ == "__main__":
    main()
