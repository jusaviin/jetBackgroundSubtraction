#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage of the script:"
  echo "$0 inputList"
  echo "inputList = Name of file containing a list of files from which histograms are projected"
  exit
fi

INPUTLIST=$1    # Name of the input file list

# For outputfile, remove the input prefix and add output prefix
while read -r FILENAME; do
  OUTPUTFILENAME="eventPlaneCorrelation/jetBackgroundHistograms_${FILENAME:27}"
  ./projectHistograms.sh $FILENAME $OUTPUTFILENAME
done < $INPUTLIST
