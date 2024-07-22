#include "AlgorithmLibrary.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

void fitJetEventPlaneVn(){

  // Open the data files
  const int nFiles = 2;
  TFile *inputFile[nFiles];
  TString jetLegendString[nFiles];
  inputFile[0] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_smallTest_genJets_2024-06-26.root");
  jetLegendString[0] = "Gen jets";
  inputFile[1] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_smallTest_defaultFlow_2024-06-26.root");
  jetLegendString[1] = "Default flow";

  // Selection for jet type
  const int nJetType = 2;
  TString jetTypeString[nJetType] = {"inclusiveJet", "leadingJet"};
  TString looseJetTypeString[nJetType] = {"Inclusive jet", "Leading jet"};
  int iJetType = 0; // Select here 0 for inclusive jet, 1 for leading jet
  
  double centralityBinBorders[] = {0,10,30,50,90};
  
  bool matchYields = true;  // For comparison purposes, match the average yields of different jet collections
  int referenceYield = 0;   // Choose which jet collection to use as basis for yield matching
    
  bool drawAllInSamePlot = true;  // True: Draw all three jet collection to the same plot. False: Use separate plots
  bool hideFit = false;
  
  bool printVs = false; // True = Print the jet vn values to the console. False = Do not do that
  
  bool saveFigures = false;
  TString saveComment = "_manualJECpfCs";
  
  // Read the histograms from the data files
  const int nCentralityBins = 4;
  const int eventPlaneOrder = 2;
  TH1D *hJetEventPlane[nFiles][nCentralityBins];
  TF1 *fitFunctionJetEventPlane[nFiles][nCentralityBins];
  double averageYield[nFiles][nCentralityBins];
  
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      hJetEventPlane[iFile][iCentrality] = (TH1D*) inputFile[iFile]->Get(Form("%sEventPlaneOrder%d_C%d", jetTypeString[iJetType].Data(), eventPlaneOrder, iCentrality));    
    } // Centrality loop
  } // File loop
  
  // Scale the histograms such that yields between different jet collections match
  if(matchYields){
    
    // Find the average yield from each histogram and use it to scale the distributions
    for(int iFile = 0; iFile < nFiles; iFile++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        hJetEventPlane[iFile][iCentrality]->Fit("pol0");
        averageYield[iFile][iCentrality] = hJetEventPlane[iFile][iCentrality]->GetFunction("pol0")->GetParameter(0);
        hJetEventPlane[iFile][iCentrality]->Scale(1 / averageYield[iFile][iCentrality]);
      } // Centrality loop
    } // File loop
    
  } // Matching yields between different files
  
  
  // Do a fourier fit up to v4 to all histograms
  AlgorithmLibrary *fitter = new AlgorithmLibrary();
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){ 
      fitter->FourierFit(hJetEventPlane[iFile][iCentrality], 4);
      fitFunctionJetEventPlane[iFile][iCentrality] = hJetEventPlane[iFile][iCentrality]->GetFunction("fourier");
    } // Centrality loop
  } // File loop
  
  // Print values of vn components to console
  if(printVs){
    for(int iFile = 0; iFile < nFiles; iFile++){
      for(int iFlow = 2; iFlow < 5; iFlow++){
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          cout << Form("%s-event plane v%d. ", looseJetTypeString[iJetType].Data(), iFlow) << jetLegendString[iFile].Data() << Form(". Cent %.0f-%.0f: ", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]) << fitFunctionJetEventPlane[iFile][iCentrality]->GetParameter(iFlow) << endl;
        } // Centrality loop
      } // Flow component loop
    } // Jet type loop
  } // Printing values of vn components
    
  JDrawer *drawer = new JDrawer();
  
  TString centralityString;
  TString compactCentralityString;
  TLegend *legend;
  int colors[] = {kBlack, kRed, kBlue, kGreen+3};
  double maxYscale, minYscale;
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
    
    if(drawAllInSamePlot){
      legend = new TLegend(0.2,0.68,0.4,1);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, Form("%s",centralityString.Data()), "");
    }
    
    maxYscale = hJetEventPlane[referenceYield][iCentrality]->GetMaximum();
    maxYscale = maxYscale + 0.2*maxYscale; // 0.1
    minYscale = hJetEventPlane[referenceYield][iCentrality]->GetMinimum();
    minYscale = minYscale - 0.02*minYscale;
    
    for(int iFile = 0; iFile < nFiles; iFile++){
      
      if(hideFit) hJetEventPlane[iFile][iCentrality]->RecursiveRemove(fitFunctionJetEventPlane[iFile][iCentrality]);
      
      hJetEventPlane[iFile][iCentrality]->Rebin(2);
      hJetEventPlane[iFile][iCentrality]->Scale(0.5);
      hJetEventPlane[iFile][iCentrality]->GetYaxis()->SetRangeUser(minYscale, maxYscale);
      if(drawAllInSamePlot) {
        hJetEventPlane[iFile][iCentrality]->SetLineColor(colors[iFile]);
        fitFunctionJetEventPlane[iFile][iCentrality]->SetLineColor(colors[iFile]);
        
      }
      
      if(!drawAllInSamePlot || iFile == 0){
      drawer->DrawHistogram(hJetEventPlane[iFile][iCentrality], "#Delta#varphi", "A.U.", " ");
      } else {
        hJetEventPlane[iFile][iCentrality]->Draw("same");
      }
      
      if(!drawAllInSamePlot){
        legend = new TLegend(0.2,0.7,0.4,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, Form("%s, %s > 80 GeV", centralityString.Data(), looseJetTypeString[iJetType].Data()), "");
        legend->AddEntry((TObject*)0, jetLegendString[iFile], "");
        legend->Draw();
        
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/%sEventPlaneDeltaPhi%s_%s%s.pdf", jetTypeString[iJetType].Data(), saveComment.Data(), jetLegendString[iFile].Data(), compactCentralityString.Data()));
        }
        
      } else {
        legend->AddEntry(hJetEventPlane[iFile][iCentrality], Form("%s, v_{%d} = %.3f", jetLegendString[iFile].Data(), eventPlaneOrder, fitFunctionJetEventPlane[iFile][iCentrality]->GetParameter(eventPlaneOrder)) ,"l");
      }
      
    } // Jet type loop
    if(drawAllInSamePlot){
      legend->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/%sEventPlaneDeltaPhi%s_jetTypeComparison%s.pdf", jetTypeString[iJetType].Data(), saveComment.Data(), compactCentralityString.Data()));
      }
    }
  } // Centrality loop
}
