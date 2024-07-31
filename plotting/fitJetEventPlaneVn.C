#include "AlgorithmLibrary.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

void fitJetEventPlaneVn(){

  // Open the data files
  std::vector<TFile*> inputFile;
  std::vector<TString> jetLegendString;
  inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_generatorJets_fullRun2MC_bugFix_2024-07-29.root"));
  jetLegendString.push_back("Gen jets");
  inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_defaultFlowJets_fullRun2MC_bugFix_2024-07-29.root"));
  jetLegendString.push_back("Default flow jets");
  inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_iterativeFlowSubtraction_minJetPt15_2024-07-30.root"));
  jetLegendString.push_back("Iterative flow, 15 GeV");
  inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_iterativeFlowSubtraction_minJetPt20_2024-07-30.root"));
  jetLegendString.push_back("Iterative flow, 20 GeV");
  inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_iterativeFlowSubtraction_minJetPt30_2024-07-30.root"));
  jetLegendString.push_back("Iterative flow, 30 GeV");

  // Find the number of files
  const int nFiles = inputFile.size();

  // Selection for jet type
  const int nJetType = 2;
  TString jetTypeString[nJetType] = {"inclusiveJet", "leadingJet"};
  TString looseJetTypeString[nJetType] = {"Inclusive jet", "Leading jet"};
  int iJetType = 1; // Select here 0 for inclusive jet, 1 for leading jet
  
  double centralityBinBorders[] = {0,10,30,50,90};
  
  bool matchYields = true;  // For comparison purposes, match the average yields of different jet collections
    
  int nRebin = 2; // Option to rebin the histograms
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
      hJetEventPlane[iFile][iCentrality] = (TH1D*) inputFile.at(iFile)->Get(Form("%sEventPlaneOrder%d_C%d", jetTypeString[iJetType].Data(), eventPlaneOrder, iCentrality));    
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
          cout << Form("%s-event plane v%d. ", looseJetTypeString[iJetType].Data(), iFlow) << jetLegendString.at(iFile).Data() << Form(". Cent %.0f-%.0f: ", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]) << fitFunctionJetEventPlane[iFile][iCentrality]->GetParameter(iFlow) << endl;
        } // Centrality loop
      } // Flow component loop
    } // Jet type loop
  } // Printing values of vn components
    
  JDrawer *drawer = new JDrawer();
  
  TString centralityString;
  TString compactCentralityString;
  TLegend *legend;
  int colors[] = {kBlack, kRed, kBlue, kGreen+3, kMagenta, kCyan};
  std::pair<double,double> yAxisZoom;
  double spread;
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
    
    // Setup a legend for the plot
    legend = new TLegend(0.2,0.68,0.4,1);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*)0, Form("%s",centralityString.Data()), "");

    // Find the minimum and maximum scales from all possible histograms
    yAxisZoom.first = 1e10;
    yAxisZoom.second = -1e10;

    // Determine a good y-axis zoom range
    for(int iFile = 0; iFile < nFiles; iFile++){
      yAxisZoom = fitter->FindHistogramMinMax(hJetEventPlane[iFile][iCentrality], yAxisZoom);
    }
    spread = yAxisZoom.second - yAxisZoom.first;
    yAxisZoom.first = yAxisZoom.first - spread * 0.1;
    yAxisZoom.second = yAxisZoom.second + spread * 0.1;
    
    // Draw all the distributions to the same canvas
    for(int iFile = 0; iFile < nFiles; iFile++){
      
      // Option to hide the Fourier fits from the histograms
      if(hideFit) hJetEventPlane[iFile][iCentrality]->RecursiveRemove(fitFunctionJetEventPlane[iFile][iCentrality]);
      
      // Option to rebin the histograms
      if(nRebin > 1){
        hJetEventPlane[iFile][iCentrality]->Rebin(nRebin);
        hJetEventPlane[iFile][iCentrality]->Scale(1.0 / nRebin);
      }

      // Set drawing range and style for the histograms
      hJetEventPlane[iFile][iCentrality]->GetYaxis()->SetRangeUser(yAxisZoom.first, yAxisZoom.second);
      hJetEventPlane[iFile][iCentrality]->SetLineColor(colors[iFile]);
      fitFunctionJetEventPlane[iFile][iCentrality]->SetLineColor(colors[iFile]);
      
      // Draw the histograms to the canvas
      if(iFile == 0){
        drawer->DrawHistogram(hJetEventPlane[iFile][iCentrality], "#Delta#varphi", "A.U.", " ");
      } else {
        hJetEventPlane[iFile][iCentrality]->Draw("same");
      }
      
      // Add a legend entry for this histogram
      legend->AddEntry(hJetEventPlane[iFile][iCentrality], Form("%s, v_{%d} = %.3f", jetLegendString.at(iFile).Data(), eventPlaneOrder, fitFunctionJetEventPlane[iFile][iCentrality]->GetParameter(eventPlaneOrder)) ,"l");
      
    } // Jet type loop

    // Draw the legend
    legend->Draw();
     
    // Option to save the figures 
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/%sEventPlaneDeltaPhi%s_jetTypeComparison%s.pdf", jetTypeString[iJetType].Data(), saveComment.Data(), compactCentralityString.Data()));
    }
  } // Centrality loop
}
