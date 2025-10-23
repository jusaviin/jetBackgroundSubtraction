#include "AlgorithmLibrary.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "JetBackgroundCard.h"
#include "JetBackgroundHistogramManager.h"

/*
 * Macro for drawing jet-event plane correlation histograms and doing a flow fit to them
 *
 *  TString inputFileList = If defined, read the input files and legend strings from this file. If not, use the manually defined file names and legend strings 
 */
void debugJetEventPlaneVn(TString inputFileList = ""){

  // Define vectors for input files and legend string corresponding to said files
  std::vector<TFile*> inputFile;
  std::vector<TString> jetLegendString;
  std::vector<TString> saveNameString;
  std::vector<std::vector<bool>> drawFlowFitDebug; // Option to draw histograms with and without flow fit flag.
  TString saveComment;
  AlgorithmLibrary *fitter = new AlgorithmLibrary();

  // If a text file is provided as input, read the input files and legend string from there. Otherwise use manually defined ones
  if(inputFileList.EndsWith(".txt")){

    std::tie(inputFile, jetLegendString, saveNameString, drawFlowFitDebug, saveComment) = fitter->ReadFileListWithDebug(inputFileList);

    // If there was an error in loading the files, exit the program
    if(inputFile.at(0) == NULL){
      cout << "File loading failed! Cannot execute the code." << endl;
      return;
    }

  } else {

    std::vector<bool> sampleVector = {false, false, true};

    inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetBackgroundHistograms_genJets_2024-08-10.root"));
    jetLegendString.push_back("Gen jets");
    drawFlowFitDebug.push_back(sampleVector);
    inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetBackgroundAnalysis_optimizedCutFlowJets_ptFromMatchedGenJet_eventPlaneProjection_2025-07-23.root"));
    jetLegendString.push_back("Optimized flow with gen p_{T}");
    drawFlowFitDebug.push_back(sampleVector);
    inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetBackgroundAnalysis_optimizedCutFlowJetsNoIter_ptFromMatchedGenJet_eventPlaneProjection_2025-07-23.root"));
    jetLegendString.push_back("No iter optimized jets with gen p_{T}");
    drawFlowFitDebug.push_back(sampleVector);
    inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetBackgroundAnalysis_defaultFlowJets_ptFromMatchedGenJet_eventPlaneProjected_2025-07-22.root"));
    jetLegendString.push_back("Default flow with gen p_{T}");
    drawFlowFitDebug.push_back(sampleVector);
  }

  // Create a vector of cards from all input files
  std::vector<JetBackgroundCard*> cardVector;
  for(auto thisFile : inputFile){
    cardVector.push_back(new JetBackgroundCard(thisFile));
  }

  TString flowFlagString[3] = {" no fit", " flow fit", ""};

  // Find the number of files
  const int nFiles = inputFile.size();

  // Selection for jet type
  const int nJetType = 2;
  TString jetTypeString[nJetType] = {"inclusiveJet", "leadingJet"};
  TString looseJetTypeString[nJetType] = {"inclusive jet", "leading jet"};
  int iJetType = 0; // Select here 0 for inclusive jet, 1 for leading jet

  // If legend says were are doing calo jets, read the calo jets from the file
  std::vector<int> jetTypeVector;
  for(TString legendString : jetLegendString){
    if(legendString.Contains("alo")){
      jetTypeVector.push_back(2);
    } else {
      jetTypeVector.push_back(iJetType);
    }
  }
  
  // Select the bin ranges that are analyzed
  // Default centrality bins: 4, 14, 34, 54, 94
  std::vector<std::pair<int,int>> analyzedCentralityBin;
  analyzedCentralityBin.push_back(std::make_pair(4,14));
  analyzedCentralityBin.push_back(std::make_pair(14,34));
  analyzedCentralityBin.push_back(std::make_pair(34,54));
  analyzedCentralityBin.push_back(std::make_pair(54,94));

  // To use the whole region defined in the files, set the second bin to 0
  // Default jet pT bins: 80, 100, 120, 140, 160, 180, 200, 300, 500, 5020
  std::vector<std::pair<int,int>> analyzedJetPtBin;
  analyzedJetPtBin.push_back(std::make_pair(80, 0));
  //analyzedJetPtBin.push_back(std::make_pair(80, 100));
  //analyzedJetPtBin.push_back(std::make_pair(120, 140));
  
  bool matchYields = true;  // For comparison purposes, match the average yields of different jet collections
    
  int nRebin = 2; // Option to rebin the histograms
  bool hideFit = false;
  
  bool printVs = false; // True = Print the jet vn values to the console. False = Do not do that
  
  bool saveFigures = true;
  if(!inputFileList.EndsWith(".txt")) saveComment = "manualJECpfCs";

  if(saveComment.CompareTo("", TString::kExact) != 0) saveComment.Prepend("_");
  
  // Fing the number of bins in the files
  const int nCentralityBins = cardVector.at(0)->GetNCentralityBins();
  const int nJetPtBins = cardVector.at(0)->GetNJetPtBins();
  const int eventPlaneOrder = 2;

  // Initialize objects to null and numbers negative
  TH1D* hJetEventPlane[nFiles][nCentralityBins][nJetPtBins+1][3];
  TF1* fitFunctionJetEventPlane[nFiles][nCentralityBins][nJetPtBins+1][3];
  double averageYield[nFiles][nCentralityBins][nJetPtBins+1][3];

  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBins+1; iJetPt++){
        for(int iFlowFlag = 0; iFlowFlag < 3; iFlowFlag++){
          hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag] = NULL;
          fitFunctionJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag] = NULL;
          averageYield[iFile][iCentrality][iJetPt][iFlowFlag] = -999;
        }
      } // Jet pT loop
    } // Centrality loop
  } // File loop

  // Initialize histogram managers from each input file
  std::vector<JetBackgroundHistogramManager*> histograms;
  for(auto thisFile : inputFile){

    // Add the histogram manager with properly loaded histograms to the manager of histogram managers
    histograms.push_back(new JetBackgroundHistogramManager(thisFile));
  }

  

  // Read the histograms from the histogram managers
  int iCentrality, iCentralityMatched;
  int iJetPt, iJetPtMatched;
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iFlowFlag = 0; iFlowFlag < 3; iFlowFlag++){
      if(!drawFlowFitDebug.at(iFile).at(iFlowFlag)) continue;
      for(auto centralityBin : analyzedCentralityBin){
        iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
        iCentralityMatched = cardVector.at(iFile)->FindBinIndexCentrality(centralityBin);
        for(auto jetPtBin : analyzedJetPtBin){
          if(jetPtBin.second == 0){
            hJetEventPlane[iFile][iCentrality][nJetPtBins][iFlowFlag] = histograms.at(iFile)->GetHistogramJetEventPlane(eventPlaneOrder, jetTypeVector.at(iFile), iCentralityMatched, -1, iFlowFlag);
          } else {
            iJetPt = cardVector.at(0)->FindBinIndexJetPt(jetPtBin);
            iJetPtMatched = cardVector.at(iFile)->FindBinIndexJetPt(jetPtBin);
            hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag] = histograms.at(iFile)->GetHistogramJetEventPlane(eventPlaneOrder, jetTypeVector.at(iFile), iCentralityMatched, iJetPtMatched, iFlowFlag);
          }
        } // Jet pT loop 
      } // Centrality loop
    } // Flow debug flag
  } // File loop
  
  // Scale the histograms such that yields between different jet collections match
  if(matchYields){
    
    // Find the average yield from each histogram and use it to scale the distributions
    for(int iFile = 0; iFile < nFiles; iFile++){
      for(int iFlowFlag = 0; iFlowFlag < 3; iFlowFlag++){
        if(!drawFlowFitDebug.at(iFile).at(iFlowFlag)) continue;
        for(auto centralityBin : analyzedCentralityBin){
          for(auto jetPtBin : analyzedJetPtBin){
            iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
            iJetPt = jetPtBin.second == 0 ? nJetPtBins : cardVector.at(0)->FindBinIndexJetPt(jetPtBin);
            hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->Fit("pol0","0");
            averageYield[iFile][iCentrality][iJetPt][iFlowFlag] = hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->GetFunction("pol0")->GetParameter(0);
            hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->Scale(1 / averageYield[iFile][iCentrality][iJetPt][iFlowFlag]);
          } // Jet pT loop
        } // Centrality loop
      } // Flow debug flag
    } // File loop
    
  } // Matching yields between different files
  
  
  // Do a fourier fit up to v4 to all histograms
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iFlowFlag = 0; iFlowFlag < 3; iFlowFlag++){
      if(!drawFlowFitDebug.at(iFile).at(iFlowFlag)) continue;
      for(auto centralityBin : analyzedCentralityBin){
        for(auto jetPtBin : analyzedJetPtBin){
          iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
          iJetPt = jetPtBin.second == 0 ? nJetPtBins : cardVector.at(0)->FindBinIndexJetPt(jetPtBin);
          fitter->FourierFit(hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag], 4, false, "0");
          fitFunctionJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag] = hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->GetFunction("fourier");
        }
      } // Centrality loop
    } // Flow debug flag
  } // File loop
  
  // Print values of vn components to console
  if(printVs){
    for(int iFile = 0; iFile < nFiles; iFile++){
      for(int iFlowFlag = 0; iFlowFlag < 3; iFlowFlag++){
        if(!drawFlowFitDebug.at(iFile).at(iFlowFlag)) continue;
        for(int iFlow = 2; iFlow < 5; iFlow++){
          for(auto centralityBin : analyzedCentralityBin){
            for(auto jetPtBin : analyzedJetPtBin){
              iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
              iJetPt = jetPtBin.second == 0 ? nJetPtBins : cardVector.at(0)->FindBinIndexJetPt(jetPtBin);
              cout << Form("%s-event plane v%d. ", looseJetTypeString[iJetType].Data(), iFlow) << jetLegendString.at(iFile).Data() << Form(". Cent %d-%d: ", centralityBin.first, centralityBin.second) << fitFunctionJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->GetParameter(iFlow) << endl;
            } // Jet pT loop
          } // Centrality loop
        } // Flow component loop
      } // Jet type loop
    } // Flow debug flag
  } // Printing values of vn components
    
  JDrawer* drawer = new JDrawer();
  
  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;
  TLegend* legend;
  int colors[] = {kBlack, kRed, kBlue, kGreen+3, kMagenta, kCyan};
  std::pair<double,double> yAxisZoom;
  double spread;
  int colorIndex = 0;
  bool firstDistribution = true;
  
  for(auto centralityBin : analyzedCentralityBin){
    iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
    
    centralityString = Form("Cent: %d-%d%%", centralityBin.first, centralityBin.second);
    compactCentralityString = Form("_C=%d-%d", centralityBin.first, centralityBin.second);

    for(auto jetPtBin : analyzedJetPtBin){

      // Clear color index
      colorIndex = 0;

      if(jetPtBin.second == 0){
        iJetPt = nJetPtBins;
        jetPtString = Form("%s p_{T} > %.0f GeV", looseJetTypeString[iJetType].Data(), cardVector.at(0)->GetLowBinBorderJetPt(0));
        compactJetPtString = Form("_J>%.0f", cardVector.at(0)->GetLowBinBorderJetPt(0));
      } else {
        iJetPt = cardVector.at(0)->FindBinIndexJetPt(jetPtBin);
        jetPtString = Form("%d < %s p_{T} < %d GeV", jetPtBin.first, looseJetTypeString[iJetType].Data(), jetPtBin.second);
        compactJetPtString = Form("_J=%d-%d", jetPtBin.first, jetPtBin.second);
      }

    
      // Setup a legend for the plot
      legend = new TLegend(0.2,0.68,0.4,1);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, Form("%s, %s",centralityString.Data(), jetPtString.Data()), "");

      // Find the minimum and maximum scales from all possible histograms
      yAxisZoom.first = 1e10;
      yAxisZoom.second = -1e10;

      // Determine a good y-axis zoom range
      for(int iFile = 0; iFile < nFiles; iFile++){
        for(int iFlowFlag = 0; iFlowFlag < 3; iFlowFlag++){
          if(!drawFlowFitDebug.at(iFile).at(iFlowFlag)) continue;
          yAxisZoom = fitter->FindHistogramMinMax(hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag], yAxisZoom);
        } // Flow debug flag
      } // File loop
      spread = yAxisZoom.second - yAxisZoom.first;
      yAxisZoom.first = yAxisZoom.first - spread * 0.1;
      yAxisZoom.second = yAxisZoom.second + spread * 0.1;
    
      // Draw all the distributions to the same canvas
      firstDistribution = true;
      for(int iFile = 0; iFile < nFiles; iFile++){
        for(int iFlowFlag = 2; iFlowFlag >= 0; iFlowFlag--){
          if(!drawFlowFitDebug.at(iFile).at(iFlowFlag)) continue;
            
          // Option to rebin the histograms
          if(nRebin > 1){
            hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->Rebin(nRebin);
            hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->Scale(1.0 / nRebin);
          }

          // Set drawing range and style for the histograms
          hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->GetYaxis()->SetRangeUser(yAxisZoom.first, yAxisZoom.second);
          hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->SetLineColor(colors[colorIndex]);
          fitFunctionJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->SetLineColor(colors[colorIndex++]);
      
          // Draw the histograms to the canvas
          if(firstDistribution){
            drawer->DrawHistogram(hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag], "#Delta#varphi", "A.U.", " ");
            firstDistribution = false;
          } else {
            hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->Draw("same");
          }
      
          // Add a legend entry for this histogram
          legend->AddEntry(hJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag], Form("%s%s, v_{%d} = %.3f", jetLegendString.at(iFile).Data(), flowFlagString[iFlowFlag].Data(), eventPlaneOrder, fitFunctionJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->GetParameter(eventPlaneOrder)) ,"l");

          // Draw the fit if not explicitly required to hide it
          if(!hideFit){
            fitFunctionJetEventPlane[iFile][iCentrality][iJetPt][iFlowFlag]->Draw("same");
          }

        } // Flow flag debug 
      } // File loop

      // Draw the legend
      legend->Draw();
     
      // Option to save the figures 
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/%sEventPlaneDeltaPhi%s_jetTypeComparison%s%s.pdf", jetTypeString[iJetType].Data(), saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data()));
      }
    } // Jet pT loop 
  } // Centrality loop
}
