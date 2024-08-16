#/bin/bash

if [ "$#" -lt 2 ]; then
  echo "Usage of the script:"
  echo "$0 inputFile outputFile [-h histograms]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  echo "-h histograms = Binary code for histograms that are loaded from the file. Default: 31"
  echo "  Bit 0 = Load event information histograms (to set: 1)"
  echo "  Bit 1 = Load jet histograms (to set: 2)"
  echo "  Bit 2 = Load jet pT closure histograms (to set: 4)"
  echo "  Bit 3 = Load jet pT response matrices (to set: 8)"
  echo "  Bit 4 = Load jet-event plane corerlation histograms (to set: 16)"
  exit
fi

INPUT=$1    # Name of the input file
OUTPUT=$2   # Name of the output file
shift 2     # Shift the arguments by 2 to read the optional arguments

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

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the projection code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/projectJetBackgroundHistograms.C

# Project all histograms
root -l -b -q 'plotting/projectJetBackgroundHistograms.C("'${INPUT}'","'${OUTPUT}'",'${HISTOGRAMS}')'

# Put the placeholder string back to the histogram projection file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/projectJetBackgroundHistograms.C
