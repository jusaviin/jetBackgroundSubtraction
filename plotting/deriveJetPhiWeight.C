#include "JetBackgroundHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JetBackgroundCard.h"
#include "JDrawer.h"

/*
 * Macro for deriving a weight that will make jet phi distribution flat
 */
void deriveJetPhiWeight(){
  
  // Open the files and check that they exist
  TString fileName = "eventPlaneCorrelation/jetBackgroundHistograms_defaultFlow_2024-08-10.root";
  TFile* inputFile = TFile::Open(fileName);
    
  if(inputFile == NULL){
    cout << "Error! The file " << fileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the eventPlaneCorrelation/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  JetBackgroundCard* card = new JetBackgroundCard(inputFile);
  
  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Read the number of centrality bins from the card
  const int nCentralityBins = card->GetNCentralityBins();

  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  JetBackgroundHistogramManager* histograms = new JetBackgroundHistogramManager(inputFile);
  
  // Histograms for jet phi
  TH1D* hJetPhi[nCentralityBins];         // Jet phi distribution
  double averageLevel[nCentralityBins];   // Average level in histograms
  TH1D* hPhiCorrection[nCentralityBins];  // Correction depending on phi
  TH1D* hCorrectedPhi[nCentralityBins];   // Corrected phi distribution
  
  // Initialize the jet histograms to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    hJetPhi[iCentrality] = nullptr;
    averageLevel[iCentrality] = 0;
    hPhiCorrection[iCentrality] = nullptr;
    hCorrectedPhi[iCentrality] = nullptr;
  }
  
  // Prepare the corrections
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){

    // Read the histograms from the file
    hJetPhi[iCentrality] = histograms->GetHistogramInclusiveJetPhi(iCentrality);

    // For each histogram, fit a zeroth order polynomial to determine the overall level
    hJetPhi[iCentrality]->Fit("pol0","Q0");
    averageLevel[iCentrality] = hJetPhi[iCentrality]->GetFunction("pol0")->GetParameter(0);

    // Prepare the correction such that in each bin, multiplying the content with the correction produces a flat phi spectrum
    hPhiCorrection[iCentrality] = (TH1D*) hJetPhi[iCentrality]->Clone(Form("phiCorrection_C%d", iCentrality));
    hPhiCorrection[iCentrality]->Reset();

    for(int iBin = 1; iBin <= hJetPhi[iCentrality]->GetNbinsX(); iBin++){
      hPhiCorrection[iCentrality]->SetBinContent(iBin, averageLevel[iCentrality] / hJetPhi[iCentrality]->GetBinContent(iBin));
    }

    // Check that multiplying with correction produces flat distribution
    hCorrectedPhi[iCentrality] = (TH1D*) hJetPhi[iCentrality]->Clone(Form("phiCorrected_C%d", iCentrality));
    hCorrectedPhi[iCentrality]->Multiply(hPhiCorrection[iCentrality]);
  }

  
  // ==========================================================================
  //           Draw the phi histograms to check if correction works
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  TLegend* legend;
  
  // ===============================
  // Draw the jet phi distributions
  // ===============================

  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
  
    // Create a legend for the figure
    legend = new TLegend(0.45,0.02,0.67,0.31);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

    // Set the style for histograms
    hJetPhi[iCentrality]->SetLineColor(kBlue);
    hCorrectedPhi[iCentrality]->SetLineColor(kRed);
  
    // Draw the histograms to the upper canvas
    drawer->DrawHistogram(hJetPhi[iCentrality], "Jet #phi", "A.U.", " ");
    legend->AddEntry(hJetPhi[iCentrality], Form("Jet #phi, Cent %d", iCentrality), "l");
  
    hCorrectedPhi[iCentrality]->Draw("same");
    legend->AddEntry(hCorrectedPhi[iCentrality], "Correction in #phi", "l");
  
    // Draw the legend
    legend->Draw();
  
    // Save the figures to a file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/jetPhiWeightCheck_C%d%s.%s", iCentrality, saveComment, figureFormat));
    }
  }

  // Save the correction histograms to a file
  TFile* outputFile = new TFile("jetPhiWeights.root", "RECREATE");
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    hPhiCorrection[iCentrality]->Write();
  }
  outputFile->Close();
}
