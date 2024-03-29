#!/bin/sh

datafile=data/beta_flipflop.csv
patientinfofile=data/patientinfo.csv
outputdir=examples
samplename=PA4
Smin=3
Smax=7
nlive=300

mkdir -p ${outputdir}

for ((i=$Smin;i<=$Smax;i++)); 
do
    inference.py $datafile $patientinfofile $outputdir $samplename $i -nlive $nlive --verbose
done 
combine_samples.py $datafile $patientinfofile $outputdir $samplename
