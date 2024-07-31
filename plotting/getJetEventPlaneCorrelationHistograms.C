#include "JetBackgroundCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)

/*
 *  Macro for projecting the jet-event plane correlation histograms from THnSparses
 *
 *  const char* inputFileName = Name of the file from which the histograms are projected
 *  const char* outputFileName = Name of the file to which the extracted histograms go
 */
void getJetEventPlaneCorrelationHistograms(const char* inputFileName, const char* outputFileName){

  // Open the data file
  TFile *inputFile = TFile::Open(inputFileName);
  
  // Configuration
  const int nEventPlaneOrder = 3;

  // Add information about the used input files to the card
  JetBackgroundCard* card = new JetBackgroundCard(inputFile);
  card->AddFileName(inputFileName);

  // The git hash here will be replaced by the latest commit hash by projectHistograms.sh script
  const char* gitHash = "GITHASHHERE";
  card->AddProjectionGitHash(gitHash);
  
  // Read the histogram with the given name from the file
  THnSparseD *inclusiveJetEventPlaneArray[nEventPlaneOrder];
  THnSparseD *leadingJetEventPlaneArray[nEventPlaneOrder];
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    inclusiveJetEventPlaneArray[iOrder] = (THnSparseD*) inputFile->Get(Form("inclusiveJetEventPlaneOrder%d", iOrder+2));
    leadingJetEventPlaneArray[iOrder] = (THnSparseD*) inputFile->Get(Form("leadingJetEventPlaneOrder%d", iOrder+2));
  }
  
  // If cannot find histogram, inform that it could not be found and return null
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    if(inclusiveJetEventPlaneArray[iOrder] == nullptr || leadingJetEventPlaneArray[iOrder] == nullptr){
      cout << "Could not find histograms of order " << iOrder << ". Will not compute." << endl;
      return;
    }
  }
  
  // Apply the restrictions in the set of axes
  const int nCentralityBins = 4;
  TH1D *hInclusiveJetEventPlane[nEventPlaneOrder][nCentralityBins];
  TH1D *hLeadingJetEventPlane[nEventPlaneOrder][nCentralityBins];
  
  TString newName;
  
  // Regular centrality based binning
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      inclusiveJetEventPlaneArray[iOrder]->GetAxis(2)->SetRange(iCentrality+2,iCentrality+2);
      newName = Form("inclusiveJetEventPlaneOrder%d_C%d", iOrder+2, iCentrality);
      
      hInclusiveJetEventPlane[iOrder][iCentrality] = (TH1D*) inclusiveJetEventPlaneArray[iOrder]->Projection(0);
      hInclusiveJetEventPlane[iOrder][iCentrality]->SetName(newName);
      
      leadingJetEventPlaneArray[iOrder]->GetAxis(2)->SetRange(iCentrality+2,iCentrality+2);
      newName = Form("leadingJetEventPlaneOrder%d_C%d", iOrder+2, iCentrality);
      
      hLeadingJetEventPlane[iOrder][iCentrality] = (TH1D*) leadingJetEventPlaneArray[iOrder]->Projection(0);
      hLeadingJetEventPlane[iOrder][iCentrality]->SetName(newName);
      
    }
  }
  
  // Reset the range in histogram arrays
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    inclusiveJetEventPlaneArray[iOrder]->GetAxis(2)->SetRange(0,0);
    leadingJetEventPlaneArray[iOrder]->GetAxis(2)->SetRange(0,0);
  }
  
  // Save the histograms to the output file
  TFile *outputFile = new TFile(outputFileName,"UPDATE");
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      hInclusiveJetEventPlane[iOrder][iCentrality]->Write("",TObject::kOverwrite);
      hLeadingJetEventPlane[iOrder][iCentrality]->Write("",TObject::kOverwrite);
    }
  }

  // Save also the card information to file
  if(!gDirectory->GetDirectory("JCard")) card->Write(outputFile);
  
  outputFile->Close();
}
