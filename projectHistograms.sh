#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage of the script:"
  echo "$0 inputFile outputFile"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  exit
fi

INPUT=$1    # Name of the input file
OUTPUT=$2   # Name of the output file

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the projection code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/projectJetBackgroundHistograms.C

# Project all histograms
root -l -b -q 'plotting/projectJetBackgroundHistograms.C("'${INPUT}'","'${OUTPUT}'",15)'

# Put the placeholder string back to the histogram projection file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/projectJetBackgroundHistograms.C
