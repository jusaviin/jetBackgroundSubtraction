# Jet background subtraction

Code to study phi-dependent background subtraction algorithms for jets

## Producing HiForests

To produce HiForests you can subsequently analyze with this code, check the instructions from https://github.com/jusaviin/cmssw/tree/flowSubtractionUpdate

## The idea behind the analysis

We are trying to improve the phi-dependent background subtraction for jets. The underlying event contains long range correlations, knows as the anisotropic flow. This created long range correlation structures in phi. Thus if we just average the undeerlying event energy density in eta and subtract it from the jet energy, we will get higher pT jets reconstructed in the upwards fluctuations of the underlying event and lower pT jets in the downwards fluctuations of the underlying event. The performance of the underlying event subtraction algorithm can be studied with the help of Pythia+Hydjet simulation. When the hard Pythia event is embedded into the soft underlying event background created by Hydjet, there is no preferred direction for the jets. This means that there is no correlation between the event plane angle calculated from the Hydjet particles, and the direction of Pythia jets. However, reconstructing the jets from all particle flow candidates will create correlations between jet direction and the event plane angle. We want to find a configuration for the background subtraction algorithm that minimizes these correlations. This is what is done in the code in this repository.

## Structure of this repository

This repository contains files needed to run the analysis described above, and to make plots from the analysis. The `src` folder contains all the files for the analysis, while the `plotting` folder has the code needed to postprocess the analysis file and produce the plots. The `crab` folder has necessary configuration to run the analysis code on CRAB. The jet energy corrections for reconstructed jets are done on analysis level, with the correction files located in `jetEnergyCorrections` folder. The files on the top directory are

`Makefile`: File for easy compilation of the analysis code

`cardJetBackground.input`: Configuration given to the analysis code

`jetBackgroundAnalysis.cxx`: Main file of the analysis code

`makeAnalysisTar.sh`: Script for making a tar ball of all analysis file for CRAB running

`projectHistograms.sh`: Script to project one dimensional histogram from the THnSparses the analysis code provides

## Running the analysis

### Local analysis

For running the analysis locally, I assume that you have first produced some HiForest files following the instruction in the other repository linked above, and downloaded some test files to your local computer. Then you will need to make a text file, that lists the paths to all the downloaded test files. Let's say the name of this file list is `testFileList.txt`.

1. Compile the code
   ```
   make
   ```
2. Run the code
   ```
   ./jetBackgroundAnalysis testFileList.txt cardJetBackground.input veryCoolData.root 0 true
   ```
   This will produce a file named `veryCoolData.root` that contains the jet-event plane correlation histograms as THnSparses. You can learn what the different arguments mean by running `./jetBackgroundAnalysis` without arguments.
3. Compile the plotting code
   ```
   cd plotting
   make
   cd ..
   ```
4. Project the one dimensional histograms from the produced files
   ```
   ./projectHistograms.sh veryCoolData.root jetEventPlaneCorrelationTest.root
   ```
   Again, you can learn the meaning of different arguments of the script by running it without arguments: `./projectHistograms.sh`.
5. Add the created `jetEventPlaneCorrelationTest.root` file as an input file in the beginning of the plotting macro `plotting/fitJetEventPlaneVn.C`. There is also some other configuration in this file that you can check to plot exactly the figures you want.
6. Run the plotting macro to see the jet-event plane correlation
   ```
   root -l plotting/fitJetEventPlaneVn.C
   ```

### CRAB analysis

For the CRAB analysis, you will need a CMSSW area on `lxplus`. Any version of CMSSW will do, since to CMSSW functionality other than CRAB is used. For example, for `CMSSW_13_3_3`, you can create this are with command `cmsrel CMSSW_13_3_3` in your work area in `lxplus`. Once this area is created, follow these instructions

1. Copy all the files from `crab` folder to your newly created CMSSW area:
   ```
   scp crab/*.* username@lxplus.cern.ch:path/to/work/area/CMSSW_13_3_3/src/flowSubtractionImprovement/
   ```
2. Package the analysis code to a tar ball and copy it to the work area on `lxplus`
   ```
   ./makeAnalysisTar.sh
   scp jetBackgroundAnalysis.tar.gz username@lxplus.cern.ch:path/to/work/area/CMSSW_13_3_3/src/flowSubtractionImprovement/
   ```
3. In your work area in `lxplus`, configure `cardJetBackground.input` and `crabFlowSubtractionStudy.py` according to your needs. 
4. Send the jobs to CRAB
   ```
   crab submit -c crabFlowSubtractionStudy.py
   ```
5. Wait until the jobs finish. Then merge and download the output files. I assume that the merged output file is named as `myOutputFile.root`
6. Project the one dimensional histograms from the produced files
   ```
   ./projectHistograms.sh myOutputFile.root jetEventPlaneCorrelation.root
   ```
   You can learn the meaning of different arguments of the script by running it without arguments: `./projectHistograms.sh`.
7. Add the created `jetEventPlaneCorrelation.root` file as an input file in the beginning of the plotting macro `plotting/fitJetEventPlaneVn.C`. There is also some other configuration in this file that you can check to plot exactly the figures you want.
8. Run the plotting macro to see the jet-event plane correlation
   ```
   root -l plotting/fitJetEventPlaneVn.C
   ```
