#include "AlgorithmLibrary.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "JetBackgroundCard.h"
#include "JetBackgroundHistogramManager.h"

/*
 * Macro for drawing jet-event plane correlation histograms and doing a flow fit to them
 *
 *  TString inputFileList = If defined, read the input files and legend strings from this file. If not, use the manually defined file names and legend strings 
 */
void fitJetEventPlaneVn(TString inputFileList = ""){

  // Define vectors for input files and legend string corresponding to said files
  std::vector<TFile*> inputFile;
  std::vector<TString> jetLegendString;
  TString saveComment;

  // If a text file is provided as input, read the input files and legend string from there. Otherwise use manually defined ones
  if(inputFileList.EndsWith(".txt")){

    // Set up the input file list for reading
    std::ifstream file_stream(inputFileList);
    std::string line;
    TObjArray* lineContents;
    TObjString* lineItem;

    bool firstLine = true;

    // Open the input file list for reading
    if( file_stream.is_open() ) {
    
      // Loop over the lines in the input file list
      while( !file_stream.eof() ) {
        getline(file_stream, line);
        TString lineString(line);
      
        // If the line is non-empty, extract the file name and legend comment from it
        if( lineString.CompareTo("", TString::kExact) != 0 ) {

          if(firstLine){
            saveComment = lineString;
            firstLine = false;
          } else {
            // Other lines define the files that are compared, and a comment given to them in legend
            lineContents = lineString.Tokenize("&"); // Tokenize the string from '&' character

            // It is assumed that the line content before '|' character gives the file name
            lineItem = (TObjString*)lineContents->At(0);
            inputFile.push_back(TFile::Open(lineItem->String()));

            // It is assumed taht the line content after '|' character gives the legend comment
            lineItem = (TObjString*)lineContents->At(1);
            jetLegendString.push_back(lineItem->String());
          }
        } // Empty line if
      
      } // Loop over lines in the file
    
    // If cannot read the file, give error and end program
    } else {
      std::cout << "Error, could not open " << inputFileList.Data() << " for reading" << std::endl;
      std::cout << "Please check the file name! Will not run the code!" << std::endl;
      return;
    }

  } else {

    inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_generatorJets_fullRun2MC_2024-07-29.root"));
    jetLegendString.push_back("Gen jets");
    inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_defaultFlowJets_fullRun2MC_2024-07-29.root"));
    jetLegendString.push_back("Default flow jets");
    //inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_iterativeFlowSubtraction_minJetPt15_2024-07-30.root"));
    //jetLegendString.push_back("Iterative flow, 15 GeV");
    //inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_iterativeFlowSubtraction_minJetPt20_2024-07-30.root"));
    //jetLegendString.push_back("Iterative flow, 20 GeV");
    //inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_iterativeFlowSubtraction_minJetPt30_2024-07-30.root"));
    //jetLegendString.push_back("Iterative flow, 30 GeV");
    inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_iterativeFlowSubtraction_minJetPt40_2024-08-01.root"));
    jetLegendString.push_back("Iterative flow, 40 GeV");
    //inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_iterativeFlowSubtraction_minJetPt50_2024-08-01.root"));
    //jetLegendString.push_back("Iterative flow, 50 GeV");
    //inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_iterativeFlowSubtraction_minJetPt60_2024-08-01.root"));
    //jetLegendString.push_back("Iterative flow, 60 GeV");
  }

  // Create a vector of cards from all input files
  std::vector<JetBackgroundCard*> cardVector;
  for(auto thisFile : inputFile){
    cardVector.push_back(new JetBackgroundCard(thisFile));
  }

  // Find the number of files
  const int nFiles = inputFile.size();

  // Selection for jet type
  const int nJetType = 2;
  TString jetTypeString[nJetType] = {"inclusiveJet", "leadingJet"};
  TString looseJetTypeString[nJetType] = {"inclusive jet", "leading jet"};
  int iJetType = 0; // Select here 0 for inclusive jet, 1 for leading jet
  
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
  TH1D* hJetEventPlane[nFiles][nCentralityBins][nJetPtBins+1];
  TF1* fitFunctionJetEventPlane[nFiles][nCentralityBins][nJetPtBins+1];
  double averageYield[nFiles][nCentralityBins][nJetPtBins+1];

  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBins+1; iJetPt++){
        hJetEventPlane[iFile][iCentrality][iJetPt] = NULL;
        fitFunctionJetEventPlane[iFile][iCentrality][iJetPt] = NULL;
        averageYield[iFile][iCentrality][iJetPt] = -999;
      } // Jet pT loop
    } // Centrality loop
  } // File loop

  // Initialize histogram managers from each input file
  std::vector<JetBackgroundHistogramManager*> histograms;
  JetBackgroundHistogramManager* manager;

  for(auto thisFile : inputFile){
    manager = new JetBackgroundHistogramManager(thisFile);

    // Load the jet-event plane correlation histograms
    manager->SetLoadJetEventPlaneHistograms(true);
    manager->LoadProcessedHistograms();

    // Add the histogram manager with properly loaded histograms to the manager of histogram managers
    histograms.push_back(manager);
  }

  

  // Read the histograms from the histogram managers
  int iCentrality, iCentralityMatched;
  int iJetPt, iJetPtMatched;
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(auto centralityBin : analyzedCentralityBin){
      iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
      iCentralityMatched = cardVector.at(iFile)->FindBinIndexCentrality(centralityBin);
      for(auto jetPtBin : analyzedJetPtBin){
        if(jetPtBin.second == 0){
          hJetEventPlane[iFile][iCentrality][nJetPtBins] = histograms.at(iFile)->GetHistogramJetEventPlane(eventPlaneOrder, iJetType, iCentralityMatched);
        } else {
          iJetPt = cardVector.at(0)->FindBinIndexJetPt(jetPtBin);
          iJetPtMatched = cardVector.at(0)->FindBinIndexJetPt(jetPtBin);
          hJetEventPlane[iFile][iCentrality][iJetPt] = histograms.at(iFile)->GetHistogramJetEventPlane(eventPlaneOrder, iJetType, iCentralityMatched, iJetPtMatched);
        }
      } // Jet pT loop 
    } // Centrality loop
  } // File loop
  
  // Scale the histograms such that yields between different jet collections match
  if(matchYields){
    
    // Find the average yield from each histogram and use it to scale the distributions
    for(int iFile = 0; iFile < nFiles; iFile++){
      for(auto centralityBin : analyzedCentralityBin){
        for(auto jetPtBin : analyzedJetPtBin){
          iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
          iJetPt = jetPtBin.second == 0 ? nJetPtBins : cardVector.at(0)->FindBinIndexJetPt(jetPtBin);
          hJetEventPlane[iFile][iCentrality][iJetPt]->Fit("pol0","0");
          averageYield[iFile][iCentrality][iJetPt] = hJetEventPlane[iFile][iCentrality][iJetPt]->GetFunction("pol0")->GetParameter(0);
          hJetEventPlane[iFile][iCentrality][iJetPt]->Scale(1 / averageYield[iFile][iCentrality][iJetPt]);
        } // Jet pT loop
      } // Centrality loop
    } // File loop
    
  } // Matching yields between different files
  
  
  // Do a fourier fit up to v4 to all histograms
  AlgorithmLibrary *fitter = new AlgorithmLibrary();
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(auto centralityBin : analyzedCentralityBin){
      for(auto jetPtBin : analyzedJetPtBin){
        iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
        iJetPt = jetPtBin.second == 0 ? nJetPtBins : cardVector.at(0)->FindBinIndexJetPt(jetPtBin);
        fitter->FourierFit(hJetEventPlane[iFile][iCentrality][iJetPt], 4, false, "0");
        fitFunctionJetEventPlane[iFile][iCentrality][iJetPt] = hJetEventPlane[iFile][iCentrality][iJetPt]->GetFunction("fourier");
      }
    } // Centrality loop
  } // File loop
  
  // Print values of vn components to console
  if(printVs){
    for(int iFile = 0; iFile < nFiles; iFile++){
      for(int iFlow = 2; iFlow < 5; iFlow++){
        for(auto centralityBin : analyzedCentralityBin){
          for(auto jetPtBin : analyzedJetPtBin){
            iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
            iJetPt = jetPtBin.second == 0 ? nJetPtBins : cardVector.at(0)->FindBinIndexJetPt(jetPtBin);
            cout << Form("%s-event plane v%d. ", looseJetTypeString[iJetType].Data(), iFlow) << jetLegendString.at(iFile).Data() << Form(". Cent %d-%d: ", centralityBin.first, centralityBin.second) << fitFunctionJetEventPlane[iFile][iCentrality][iJetPt]->GetParameter(iFlow) << endl;
          } // Jet pT loop
        } // Centrality loop
      } // Flow component loop
    } // Jet type loop
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
  
  for(auto centralityBin : analyzedCentralityBin){
    iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
    
    centralityString = Form("Cent: %d-%d%%", centralityBin.first, centralityBin.second);
    compactCentralityString = Form("_C=%d-%d", centralityBin.first, centralityBin.second);

    for(auto jetPtBin : analyzedJetPtBin){

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
        yAxisZoom = fitter->FindHistogramMinMax(hJetEventPlane[iFile][iCentrality][iJetPt], yAxisZoom);
      }
      spread = yAxisZoom.second - yAxisZoom.first;
      yAxisZoom.first = yAxisZoom.first - spread * 0.1;
      yAxisZoom.second = yAxisZoom.second + spread * 0.1;
    
      // Draw all the distributions to the same canvas
      for(int iFile = 0; iFile < nFiles; iFile++){
            
        // Option to rebin the histograms
        if(nRebin > 1){
          hJetEventPlane[iFile][iCentrality][iJetPt]->Rebin(nRebin);
          hJetEventPlane[iFile][iCentrality][iJetPt]->Scale(1.0 / nRebin);
        }

        // Set drawing range and style for the histograms
        hJetEventPlane[iFile][iCentrality][iJetPt]->GetYaxis()->SetRangeUser(yAxisZoom.first, yAxisZoom.second);
        hJetEventPlane[iFile][iCentrality][iJetPt]->SetLineColor(colors[iFile]);
        fitFunctionJetEventPlane[iFile][iCentrality][iJetPt]->SetLineColor(colors[iFile]);
      
        // Draw the histograms to the canvas
        if(iFile == 0){
          drawer->DrawHistogram(hJetEventPlane[iFile][iCentrality][iJetPt], "#Delta#varphi", "A.U.", " ");
        } else {
          hJetEventPlane[iFile][iCentrality][iJetPt]->Draw("same");
        }
      
        // Add a legend entry for this histogram
        legend->AddEntry(hJetEventPlane[iFile][iCentrality][iJetPt], Form("%s, v_{%d} = %.3f", jetLegendString.at(iFile).Data(), eventPlaneOrder, fitFunctionJetEventPlane[iFile][iCentrality][iJetPt]->GetParameter(eventPlaneOrder)) ,"l");

        // Draw the fit if not explicitly required to hide it
        if(!hideFit){
          fitFunctionJetEventPlane[iFile][iCentrality][iJetPt]->Draw("same");
        }
      
      } // Jet type loop

      // Draw the legend
      legend->Draw();
     
      // Option to save the figures 
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/%sEventPlaneDeltaPhi%s_jetTypeComparison%s%s.pdf", jetTypeString[iJetType].Data(), saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data()));
      }
    } // Jet pT loop 
  } // Centrality loop
}
