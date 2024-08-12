#include "JetBackgroundCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JetBackgroundHistogramManager.h"

#include <bitset>

/*
 * Macro for projecting the histograms needed in the jet background analysis from the THnSparses
 *
 *  Arguments:
 *   TString inputFileName = File from which the histograms are read
 *   const char* outputFileName = If we are producing output file, name of the output file
 *   int histogramSelection = If > 0, select a preset group of histograms. Intended to be used for easier production of output files.
 */
void projectJetBackgroundHistograms(TString inputFileName, const char* outputFileName, int histogramSelection){

  // Print the file name to console
  cout << "Projecting histograms from " << inputFileName.Data() << endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";
  
  // Choose which figure sets to draw
  bool loadEventInformation = false;
  bool loadJets = false;
  bool loadJetPtClosure = false;
  bool loadJetEventPlaneHistograms = false;
  
  /*
   * Loading only selected histograms. Done with bitwise check of an integer
   *
   *  Bit 0 = Load event information histograms (to set: 1)
   *  Bit 1 = Load jet histograms (to set: 2)
   *  Bit 2 = Load jet pT closure histograms (to set: 4)
   *  Bit 3 = Load jet-event plane corerlation histograms (to set: 8)
   */
  if(histogramSelection > 0){
    std::bitset<4> bitChecker(histogramSelection);
    loadEventInformation = bitChecker.test(0);
    loadJets = bitChecker.test(1);
    loadJetPtClosure = bitChecker.test(2);
    loadJetEventPlaneHistograms = bitChecker.test(3);
  }
  
  // ====================================================
  //  Binning configuration for the projected histograms
  // ====================================================
  
  // Option to read all the binning information from JetBackgroundCard used to create the file
  const bool readCentralityBinsFromFile = true;
  const bool readJetPtBinsFromFile = true;
  
  // If not reading the bins from the file, manually define new bin borders
  const int nCentralityBins = 4;
  const int nJetPtBins = 7;
  double centralityBinBorders[nCentralityBins+1] = {4,14,34,54,94};   // Bin borders for centrality
  double jetPtBinBorders[nJetPtBins+1] = {120,140,160,180,200,300,500,5020}; // Bin borders for jet pT in jet-event plane correlation histograms
  
  // Projected bin range
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnJetPtBin = 0;
  int lastDrawnJetPtBin = nJetPtBins-1;
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Open the input file
  TFile* inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  JetBackgroundCard* card = new JetBackgroundCard(inputFile);

  // If we manually define bin borders, check that they match with the ones on the file or give a warning
  bool binFound;
  if(!readCentralityBinsFromFile){
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      binFound = false;
      for(int iCardBin = 0; iCardBin < card->GetNCentralityBins(); iCardBin++){
        if(TMath::Abs(centralityBinBorders[iCentrality] - card->GetLowBinBorderCentrality(iCardBin)) < 0.2501){
          binFound = true;
          break;
        }
        if(TMath::Abs(centralityBinBorders[iCentrality] - card->GetHighBinBorderCentrality(iCardBin)) < 0.2501){
          binFound = true;
          break;
        }
      } // Card centrality bin loop 
      if(!binFound){
        cout << "WARNING: " << centralityBinBorders[iCentrality] << " is not a centrality bin border in the original file." << endl;
        cout << "Original centrality bin borders are: ";
        for(int iCardBin = 0; iCardBin < card->GetNCentralityBins(); iCardBin++){
          cout << card->GetLowBinBorderCentrality(iCardBin) << ", ";
        }
        cout << card->GetHighBinBorderCentrality(card->GetNCentralityBins()-1) << endl;
        cout << "If you are sure you know what you are doing, you can proceed. Otherwise make sure the bin borders match." << endl;
      } // Message to be printed if bin is not found
    } // Centrality loop
  }

  // If we change the binning, save the new binning to the card
  if(!readCentralityBinsFromFile) card->AddVector(JetBackgroundCard::kCentralityBinEdges,nCentralityBins+1, centralityBinBorders);
  if(!readJetPtBinsFromFile) card->AddVector(JetBackgroundCard::kJetPtBinEdges,nJetPtBins+1, jetPtBinBorders);
  
  // Add information about the used input files to the card
  card->AddFileName(inputFileName);
  
  // The git hash here will be replaced by the latest commit hash by projectHistograms.sh script
  const char* gitHash = "GITHASHHERE";
  card->AddProjectionGitHash(gitHash);
  
  // ==================================== //
  //     JetBackgroundHistogramManager    //
  // ==================================== //
    
  // Create and setup a new histogram manager to project and handle the histograms
  JetBackgroundHistogramManager* histograms = new JetBackgroundHistogramManager(inputFile, card);
  
  // Set which histograms to project from the input file
  histograms->SetLoadEventInformation(loadEventInformation);
  histograms->SetLoadJetHistograms(loadJets);
  histograms->SetLoad2DHistograms(true);
  histograms->SetLoadJetPtClosureHistograms(loadJetPtClosure);
  histograms->SetLoadJetEventPlaneHistograms(loadJetEventPlaneHistograms);

  // Set the binning information
  histograms->SetCentralityBins(readCentralityBinsFromFile,nCentralityBins,centralityBinBorders,true);
  if(!readCentralityBinsFromFile) histograms->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  histograms->SetJetPtBins(readJetPtBinsFromFile,nJetPtBins,jetPtBinBorders,true);
  if(!readJetPtBinsFromFile) histograms->SetJetPtBinRange(firstDrawnJetPtBin,lastDrawnJetPtBin);
  
  // Project the one dimensional histograms from the THnSparses
  histograms->LoadHistograms();
  
  // Save the histograms to an output file
  histograms->Write(outputFileName,fileWriteMode);
  
}
