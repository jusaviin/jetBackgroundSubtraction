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
  
  // Find the number of centrality and jet pT bins, and create one dimensional histograms
  const int nCentralityBins = inclusiveJetEventPlaneArray[0]->GetAxis(2)->GetNbins();
  const int nJetPtBins = inclusiveJetEventPlaneArray[0]->GetAxis(1)->GetNbins();
  TH1D *hInclusiveJetEventPlane[nEventPlaneOrder][nCentralityBins][nJetPtBins+1];
  TH1D *hLeadingJetEventPlane[nEventPlaneOrder][nCentralityBins][nJetPtBins+1];
  
  TString newName;
  
  // Project the histograms in the defined centrality and jet pT bins
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      inclusiveJetEventPlaneArray[iOrder]->GetAxis(2)->SetRange(iCentrality+1,iCentrality+1);
      newName = Form("inclusiveJetEventPlaneOrder%d_C%d", iOrder+2, iCentrality);
      
      hInclusiveJetEventPlane[iOrder][iCentrality][nJetPtBins] = (TH1D*) inclusiveJetEventPlaneArray[iOrder]->Projection(0);
      hInclusiveJetEventPlane[iOrder][iCentrality][nJetPtBins]->SetName(newName);
      
      leadingJetEventPlaneArray[iOrder]->GetAxis(2)->SetRange(iCentrality+1,iCentrality+1);
      newName = Form("leadingJetEventPlaneOrder%d_C%d", iOrder+2, iCentrality);
      
      hLeadingJetEventPlane[iOrder][iCentrality][nJetPtBins] = (TH1D*) leadingJetEventPlaneArray[iOrder]->Projection(0);
      hLeadingJetEventPlane[iOrder][iCentrality][nJetPtBins]->SetName(newName);

      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){

        inclusiveJetEventPlaneArray[iOrder]->GetAxis(1)->SetRange(iJetPt+1,iJetPt+1);
        newName = Form("inclusiveJetEventPlaneOrder%d_C%dJ%d", iOrder+2, iCentrality, iJetPt);
      
        hInclusiveJetEventPlane[iOrder][iCentrality][iJetPt] = (TH1D*) inclusiveJetEventPlaneArray[iOrder]->Projection(0);
        hInclusiveJetEventPlane[iOrder][iCentrality][iJetPt]->SetName(newName);

        leadingJetEventPlaneArray[iOrder]->GetAxis(1)->SetRange(iJetPt+1,iJetPt+1);
        newName = Form("leadingJetEventPlaneOrder%d_C%dJ%d", iOrder+2, iCentrality, iJetPt);
      
        hLeadingJetEventPlane[iOrder][iCentrality][iJetPt] = (TH1D*) leadingJetEventPlaneArray[iOrder]->Projection(0);
        hLeadingJetEventPlane[iOrder][iCentrality][iJetPt]->SetName(newName);
        
      }

      // Reset the jet pT range in the histogram arrays
      inclusiveJetEventPlaneArray[iOrder]->GetAxis(1)->SetRange(0,0);
      leadingJetEventPlaneArray[iOrder]->GetAxis(1)->SetRange(0,0);
    }
  }
  
  // Reset the centrality range in histogram arrays
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    inclusiveJetEventPlaneArray[iOrder]->GetAxis(2)->SetRange(0,0);
    leadingJetEventPlaneArray[iOrder]->GetAxis(2)->SetRange(0,0);
  }
  
  // Save the histograms to the output file
  TFile *outputFile = new TFile(outputFileName,"UPDATE");
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      hInclusiveJetEventPlane[iOrder][iCentrality][nJetPtBins]->Write("",TObject::kOverwrite);
      hLeadingJetEventPlane[iOrder][iCentrality][nJetPtBins]->Write("",TObject::kOverwrite);
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        hInclusiveJetEventPlane[iOrder][iCentrality][iJetPt]->Write("",TObject::kOverwrite);
        hLeadingJetEventPlane[iOrder][iCentrality][iJetPt]->Write("",TObject::kOverwrite);
      }
    }
  }

  // Save also the card information to file
  if(!gDirectory->GetDirectory("JCard")) card->Write(outputFile);
  
  outputFile->Close();
}
