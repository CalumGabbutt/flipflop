#!/bin/sh

datafile=data/beta_ticktock.csv
patientinfofile=data/patientinfo.csv
outputdir=examples
samplename=PA5
Smin=3
Smax=5
nlive=300
verbose=False

mkdir -p ${outputdir}

for ((S=$Smin;i<=$Smax;i++)); 
do
    python3 inference.py $datafile $patientinfofile $outputdir $samplename $S --nlive $nlive --verbose $verbose
done 
