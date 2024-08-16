#!/bin/bash

if [ "$#" -lt 1 ]; then
  echo "Usage of the script:"
  echo "$0 inputList [-h histograms]"
  echo "inputList = Name of file containing a list of files from which histograms are projected"
  echo "-h histograms = Binary code for histograms that are loaded from the file. Default: 31"
  echo "  Bit 0 = Load event information histograms (to set: 1)"
  echo "  Bit 1 = Load jet histograms (to set: 2)"
  echo "  Bit 2 = Load jet pT closure histograms (to set: 4)"
  echo "  Bit 3 = Load jet pT response matrices (to set: 8)"
  echo "  Bit 4 = Load jet-event plane corerlation histograms (to set: 16)"
  exit
  exit
fi

INPUTLIST=$1    # Name of the input file list
shift           # Shift the argument to read the optional arguments

# Read the optional arguments
while getopts ":h:" opt; do
case $opt in
h) HISTOGRAMS="$OPTARG"
;;
\?) echo "Invalid option -$OPTARG" >&2
exit 1
;;
esac
done

# Set default values to optional arguments if they are not given
HISTOGRAMS=${HISTOGRAMS:-31}

# For outputfile, remove the input prefix and add output prefix
while read -r FILENAME; do
  OUTPUTFILENAME="eventPlaneCorrelation/jetBackgroundHistograms_${FILENAME:27}"
  ./projectHistograms.sh $FILENAME $OUTPUTFILENAME -h $HISTOGRAMS  
done < $INPUTLIST
