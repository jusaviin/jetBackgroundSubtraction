#include "AlgorithmLibrary.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "JetBackgroundCard.h"
#include "JetBackgroundHistogramManager.h"

// Capitalize all characters after a space
void capitalizeAfterSpace(TString &str){
  bool capitalizeNext = false;

  // Loop over all characters in the string
  for (int i = 0; i < str.Length(); i++){

    // If the flag to capitalize a character is set, check that the character we try to capitalize is a letter
    if(capitalizeNext && isalpha(str[i])){
      // If it is a letter, do capitalization
      str[i] = toupper(str[i]);
      capitalizeNext = false;
    } else if (str[i] == ' ') {
      // If we found a space, set up a flag that the next letter should be capitalized
      capitalizeNext = true;
    } // If for capitalization
  } // Character loop in the string
}

/*
 * Macro for drawing jet pT response matrices from different configurations
 *
 *  TString inputFileList = If defined, read the input files and legend strings from this file. If not, use the manually defined file names and legend strings 
 */
void drawJetPtResponseMatrix(TString inputFileList = ""){

  // Define vectors for input files and legend string corresponding to said files
  std::vector<TFile*> inputFile;
  std::vector<TString> jetLegendString;
  std::vector<TString> saveNameString;
  TString saveComment;
  AlgorithmLibrary *jackOfAllTrades = new AlgorithmLibrary();

  // If a text file is provided as input, read the input files and legend string from there. Otherwise use manually defined ones
  if(inputFileList.EndsWith(".txt")){

    std::tie(inputFile, jetLegendString, saveNameString, saveComment) = jackOfAllTrades->ReadFileList(inputFileList);
    cout << "We did read the stuff" << endl;

    // If there was an error in loading the files, exit the program
    if(inputFile.at(0) == NULL){
      cout << "File loading failed! Cannot execute the code." << endl;
      return;
    }

  } else {

    inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetBackgroundHistograms_defaultFlow_2024-08-10.root"));
    jetLegendString.push_back("Default flow");
    saveNameString.push_back("DefaultFlow");
    inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetBackgroundHistograms_noProbabilityCuts_2024-08-09.root"));
    jetLegendString.push_back("No probability cuts");
    saveNameString.push_back("NoProbabilityCuts");
  }

  // Create a vector of cards from all input files
  std::vector<JetBackgroundCard*> cardVector;
  for(auto thisFile : inputFile){
    cardVector.push_back(new JetBackgroundCard(thisFile));
  }

  // Find the number of files
  const int nFiles = inputFile.size();
  cout << "Found " << nFiles << " files" << endl;
  
  // Select the analyzed centrality bin range
  // Default centrality bins: 4, 14, 34, 54, 94
  std::vector<std::pair<int,int>> analyzedCentralityBin;
  analyzedCentralityBin.push_back(std::make_pair(4,14));
  analyzedCentralityBin.push_back(std::make_pair(14,34));
  analyzedCentralityBin.push_back(std::make_pair(34,54));
  analyzedCentralityBin.push_back(std::make_pair(54,94));

  bool normalizeRows = true;  // Normalize reco distribution for each gen bin to 1
  bool logZ = true; // Logarithmic z-axis
  
  bool saveFigures = true;
  
  // Fing the number of bins in the files
  const int nCentralityBins = cardVector.at(0)->GetNCentralityBins();

  // Initialize the response matrix histograms to NULL
  TH2D* hResponseMatrix[nFiles][nCentralityBins];
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      hResponseMatrix[iFile][iCentrality] = NULL;
    } // Centrality loop
  } // File loop

  // Initialize histogram managers from each input file
  std::vector<JetBackgroundHistogramManager*> histograms;
  JetBackgroundHistogramManager* manager;

  for(auto thisFile : inputFile){
    manager = new JetBackgroundHistogramManager(thisFile);

    // Load the jet pT response matrices
    manager->SetLoadJetPtResponseMatrix(true);
    manager->LoadProcessedHistograms();

    // Add the histogram manager with properly loaded histograms to the manager of histogram managers
    histograms.push_back(manager);
  }

  // Read the histograms from the histogram managers
  int iCentrality, iCentralityMatched;
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(auto centralityBin : analyzedCentralityBin){
      iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
      iCentralityMatched = cardVector.at(iFile)->FindBinIndexCentrality(centralityBin);

      hResponseMatrix[iFile][iCentrality] = histograms.at(iFile)->GetHistogramJetPtResponseMatrix(iCentralityMatched);

      // If selected, normalize the matrix
      if(normalizeRows){
        cout << "Normalizinf rows" << endl;
        jackOfAllTrades->NormalizeRows(hResponseMatrix[iFile][iCentrality]);
        cout << "Survived" << endl;
      }
    } // Centrality loop
  } // File loop  

  // ======================== //
  //     Draw the figures     //
  // ======================== //
    
  JDrawer* drawer = new JDrawer();
  drawer->SetRightMargin(0.12);
  drawer->SetLeftMargin(0.15);
  drawer->SetLogZ(logZ);
  
  TString centralityString;
  TString compactCentralityString;
  
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(auto centralityBin : analyzedCentralityBin){
      iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
    
      centralityString = Form("Cent: %d-%d%%", centralityBin.first, centralityBin.second);
      compactCentralityString = Form("_C=%d-%d", centralityBin.first, centralityBin.second);

      hResponseMatrix[iFile][iCentrality]->GetZaxis()->SetRangeUser(1e-5,0.4);

      // Draw the histogram to a canvas
      drawer->DrawHistogram(hResponseMatrix[iFile][iCentrality], "Reco jet p_{T} (GeV)", "Gen jet p_{T} (GeV)", Form("%s, %s",centralityString.Data(), jetLegendString.at(iFile).Data()), "colz");
     
      // Option to save the figures 
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetPtResponseMatrix%s%s.pdf", saveNameString.at(iFile).Data(), compactCentralityString.Data()));
      }
    } // Centrality loop
  } // File loop
}
