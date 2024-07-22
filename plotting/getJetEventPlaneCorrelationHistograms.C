void getJetEventPlaneCorrelationHistograms(){

  // Open the data file
  TFile *inputFile = TFile::Open("veryCoolData.root");
  
  // Configuration
  const int nEventPlaneOrder = 3;
  
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
  
  // Save the histogram to a file
  TFile *outputFile = new TFile("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_smallTest_genJets_2024-06-26.root","UPDATE");
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      hInclusiveJetEventPlane[iOrder][iCentrality]->Write("",TObject::kOverwrite);
      hLeadingJetEventPlane[iOrder][iCentrality]->Write("",TObject::kOverwrite);
    }
  }
  
  outputFile->Close();
}
