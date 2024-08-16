#include "JetBackgroundHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "../src/JetBackgroundHistograms.h"
#include "JetBackgroundCard.h"
#include "JDrawer.h"
#include <tuple>

/*
 * Fit a Gauss function to a histogram and extract parameters from that
 *
 *  TH1* histogram = Histogram from which drawing range in searched
 *
 *  return: Gauss mean, Gauss sigma, Error for Gauss mean, Error for Gauss sigma
 */
std::tuple<double,double,double,double> fitGauss(TH1* histogram, TString title = "", TString jetTypeString = "", TString centralityBin = "", TString ptBin = "",  TString saveName = ""){
  histogram->Fit("gaus","","",0.5,1.5);
  TF1* gaussFit = histogram->GetFunction("gaus");
  
  double gaussMean = 0;
  double gaussSigma = 0;
  double gaussMeanError = 0;
  double gaussSigmaError = 0;
  
  if(gaussFit){
    gaussMean = gaussFit->GetParameter(1);
    gaussSigma = gaussFit->GetParameter(2);
    gaussMeanError = gaussFit->GetParError(1);
    gaussSigmaError = gaussFit->GetParError(2);
  }
  
  // If title is given, print the fit
  if(!title.EqualTo("")){
    JDrawer *temporaryDrawer = new JDrawer();
    temporaryDrawer->DrawHistogram(histogram,"Reco p_{T} / Gen p_{T}","Counts", " ");
    TLegend *legend = new TLegend(0.57,0.68,0.8,0.93);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->SetHeader(title);
    legend->AddEntry((TObject*)0, jetTypeString, "");
    legend->AddEntry((TObject*)0, centralityBin, "");
    legend->AddEntry((TObject*)0, ptBin, "");
    legend->Draw();
    
    if(!saveName.EqualTo("")){
      gPad->GetCanvas()->SaveAs(saveName);
    }
    
  }
  
  return std::make_tuple(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError);
}

/*
 * Draw a closure histogram with quark/gluon discrimination to the plot
 *
 *  std::vector<TH1D*> histogram = Vector of histograms to be drawn
 *  JetBackgroundCard* card = Card with binning information
 *  const char* xTitle = Title given to the x-axis
 *  const char* yTitle = Title given to the y-axis
 *  int iCentrality = Index of the centrality bin
 *  int legendNzoom = Define the y-axis zoom and the legend position in the plot
 *  TString legendTitle = Title given to the legend
 *  std::vector<TString> legendString = Vector of strings to be put into legend for each histogram
 *  const char* saveComment = Comment given to the save name file
 *  bool saveFigures = Choose whether to save the figures or not
 */
void drawClosureHistogram(std::vector<TH1D*> histogram, JetBackgroundCard* card, const char* xTitle, const char* yTitle, int iCentrality, int legendNzoom, TString legendTitle, std::vector<TString> legendString, const char* saveComment, bool saveFigures){
  
  // Create a new drawer and define bin borders and drawing style
  JDrawer* drawer = new JDrawer();
  drawer->SetCanvasSize(700,700);
  const char* centralityString;
  const char* centralitySaveName;
  int lineColors[] = {kBlack, kRed, kBlue, kGreen+3, kMagenta, kCyan};
  
  // Zooming and legend position options
  double yZoomLow = 0.9;
  double yZoomHigh = 1.1;
  double legendX1 = 0.51;
  double legendX2 = 0.93;
  double legendY1 = 0.65;
  double legendY2 = 0.92;
  
  if(legendNzoom == 1){
    yZoomLow = 0;
    yZoomHigh = 0.34;
    legendX1 = 0.51;
    legendX2 = 0.93;
    legendY1 = 0.65;
    legendY2 = 0.92;
  } else if(legendNzoom == 2){
    legendX1 = 0.16;
    legendX2 = 0.58;
    legendY1 = 0.78;
    legendY2 = 0.99;
  }

  const int nHistograms = histogram.size();

  // Draw all the histograms to the same canvas
  TLegend* legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);

  centralityString = Form(", Cent:%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
  centralitySaveName = Form("_C=%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));

  legend->SetHeader(Form("%s%s", legendTitle.Data(), centralityString));

  for(int iHistogram = 0; iHistogram < nHistograms; iHistogram++){
  
    // Set a good style for histograms and draw them
    histogram.at(iHistogram)->SetLineColor(lineColors[iHistogram]);
    histogram.at(iHistogram)->GetYaxis()->SetRangeUser(yZoomLow,yZoomHigh);
    if(iHistogram == 0){
      drawer->DrawHistogram(histogram.at(iHistogram), xTitle, yTitle, " ");
    } else {
      histogram.at(iHistogram)->Draw("same");
    }
  
    // Add a legend for the histograms
    legend->AddEntry(histogram.at(iHistogram), legendString.at(iHistogram), "l");
  
  } // Histogram loop
  
  // Draw the legend
  legend->Draw();
  
  // Draw a dashed line at one
  double startPoint = histogram.at(0)->GetXaxis()->GetBinLowEdge(1);
  double endPoint = histogram.at(0)->GetXaxis()->GetBinUpEdge(histogram.at(0)->GetNbinsX());
  TLine* oneLine = new TLine(startPoint, 1, endPoint, 1);
  oneLine->SetLineStyle(2);
  oneLine->SetLineColor(kBlack);
  oneLine->Draw("same");
  
  // Save the figures if selected to do so
  if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/jet%s%s.pdf",saveComment,centralitySaveName));
  }
  
}

/*
 * Macro for constructing jet pT closures and comparing them between different files
 */
void constructJetPtClosures(TString inputFileList = ""){

  // ==================================================================
  // =========================  Input files  ==========================
  // ==================================================================

  // Define vectors for input files and legend string corresponding to said files
  std::vector<TFile*> inputFile;
  std::vector<TString> jetLegendString;
  std::vector<TString> saveNameString;
  TString saveComment;
  AlgorithmLibrary *jackOfAllTrades = new AlgorithmLibrary();

  // If a text file is provided as input, read the input files and legend string from there. Otherwise use manually defined ones
  if(inputFileList.EndsWith(".txt")){

    std::tie(inputFile, jetLegendString, saveNameString, saveComment) = jackOfAllTrades->ReadFileList(inputFileList);

    // If there was an error in loading the files, exit the program
    if(inputFile.at(0) == NULL){
      cout << "File loading failed! Cannot execute the code." << endl;
      return;
    }

  } else {

    inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetBackgroundHistograms_defaultFlow_2024-08-10.root"));
    jetLegendString.push_back("Default flow");
    //inputFile.push_back(TFile::Open("eventPlaneCorrelation/jetBackgroundHistograms_noProbabilityCuts_2024-08-09.root"));
    //jetLegendString.push_back("No probability cuts");
    saveComment = "defaultFlow";
  }

  // Create a vector of cards from all input files
  std::vector<JetBackgroundCard*> cardVector;
  for(auto thisFile : inputFile){
    cardVector.push_back(new JetBackgroundCard(thisFile));
  }

  // Find the number of files
  const int nFiles = inputFile.size();

  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Select the analyzed centrality bin range
  // Default centrality bins: 4, 14, 34, 54, 94
  std::vector<std::pair<int,int>> analyzedCentralityBin;
  analyzedCentralityBin.push_back(std::make_pair(4,14));
  //analyzedCentralityBin.push_back(std::make_pair(14,34));
  //analyzedCentralityBin.push_back(std::make_pair(34,54));
  //analyzedCentralityBin.push_back(std::make_pair(54,94));

  
  bool drawPtClosure = true;
  bool drawEtaClosure = false;
  bool drawPhiClosure = false;
  
  bool includeQuarkGluon = (nFiles == 1); // Include only quark and only gluon jet curves is only one file is provided
  bool drawGaussFitsPt = false;
    
  bool fitResolution = false;  // Fit the jet pT resolution histograms
  
  bool saveFigures = false;  // Save the figures to file
  
  // ==================================================================
  // =================== Configuration ready ==========================
  // ==================================================================

  // Initialize histogram managers from each input file
  std::vector<JetBackgroundHistogramManager*> closureHistograms;
  JetBackgroundHistogramManager* manager;

  for(auto thisFile : inputFile){
    manager = new JetBackgroundHistogramManager(thisFile);

    // Load the jet-event plane correlation histograms
    manager->SetLoadJetPtClosureHistograms(true);
    manager->LoadProcessedHistograms();

    // Add the histogram manager with properly loaded histograms to the manager of histogram managers
    closureHistograms.push_back(manager);
  }
  
  // Find the correct number of centrality bins
  const int nCentralityBins = cardVector.at(0)->GetNCentralityBins();
  
  // Initialize reco/gen ratio and closure histograms
  TH1D* hRecoGenRatio[nFiles][JetBackgroundHistogramManager::knGenJetPtBins][nCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1];
  TH1D* hRecoGenRatioEta[nFiles][JetBackgroundHistogramManager::knJetEtaBins][nCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1];
  TH1D* hRecoGenRatioPhi[nFiles][JetBackgroundHistogramManager::knJetPhiBins][nCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1];
  TH1D* hJetPtClosure[nFiles][nCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1];
  TH1D* hJetPtClosureSigma[nFiles][nCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1];
  TH1D* hJetPtClosureEta[nFiles][nCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1];
  TH1D* hJetPtClosureSigmaEta[nFiles][nCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1];
  TH1D* hJetPtClosurePhi[nFiles][nCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1];
  TH1D* hJetPtClosureSigmaPhi[nFiles][nCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1];
  const char* histogramNamer;

    // Prepare a vector with indices for different jet types
  std::vector<int> partonIndices;
  partonIndices.push_back(JetBackgroundHistograms::knInitialPartonTypes);
  if(includeQuarkGluon){
    partonIndices.push_back(JetBackgroundHistograms::kQuark);
    partonIndices.push_back(JetBackgroundHistograms::kGluon);
  }
  
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes+1; iParton++){
        for(int iGenJetPt = 0; iGenJetPt < JetBackgroundHistogramManager::knGenJetPtBins; iGenJetPt++){
          hRecoGenRatio[iFile][iGenJetPt][iCentrality][iParton] = NULL;
        } // Generator level jet pT loop
        for(int iJetEta = 0; iJetEta < JetBackgroundHistogramManager::knJetEtaBins; iJetEta++){
          hRecoGenRatioEta[iFile][iJetEta][iCentrality][iParton] = NULL;
        }
        for(int iJetPhi = 0; iJetPhi < JetBackgroundHistogramManager::knJetPhiBins; iJetPhi++){
          hRecoGenRatioPhi[iFile][iJetPhi][iCentrality][iParton] = NULL;
        }
        histogramNamer = Form("jetPtClosure%d_Cent%d_Part%d", iFile, iCentrality, iParton);
        hJetPtClosure[iFile][iCentrality][iParton] = new TH1D(histogramNamer,histogramNamer,45,50,500);
        histogramNamer = Form("jetPtClosureSigma%d_Cent%d_Part%d", iFile, iCentrality, iParton);
        hJetPtClosureSigma[iFile][iCentrality][iParton] = new TH1D(histogramNamer,histogramNamer,45,50,500);
        histogramNamer = Form("jetPtClosureEta%d_Cent%d_Part%d", iFile, iCentrality, iParton);
        hJetPtClosureEta[iFile][iCentrality][iParton] = new TH1D(histogramNamer,histogramNamer,50,-2.5,2.5);
        histogramNamer = Form("jetPtClosureSigmaEta%d_Cent%d_Part%d", iFile, iCentrality, iParton);
        hJetPtClosureSigmaEta[iFile][iCentrality][iParton] = new TH1D(histogramNamer,histogramNamer,50,-2.5,2.5);
        histogramNamer = Form("jetPtClosurePhi%d_Cent%d_Part%d", iFile, iCentrality, iParton);
        hJetPtClosurePhi[iFile][iCentrality][iParton] = new TH1D(histogramNamer,histogramNamer,64,-TMath::Pi(),TMath::Pi());
        histogramNamer = Form("jetPtClosureSigmaPhi%d_Cent%d_Part%d", iFile, iCentrality, iParton);
        hJetPtClosureSigmaPhi[iFile][iCentrality][iParton] = new TH1D(histogramNamer,histogramNamer,64,-TMath::Pi(),TMath::Pi());
      } // Closure particle loop (quark/gluon/no selection)
    } // Centrality loop
  } // File loop
  
  JDrawer* drawer = new JDrawer();
  TF1* gaussFit;
  double gaussMean = 0;
  double gaussSigma = 0;
  double gaussMeanError = 0;
  double gaussSigmaError = 0;
  int minGenPt = 7;  // Set this to 7 to skip bins below 120 GeV
  TString genPtString;
  TString centralityString;
  TString jetTypeName[JetBackgroundHistograms::knInitialPartonTypes+1] = {"Quark", "Gluon", "Undetermined", "All"};
  TString jetTypeString;
  TString gaussFitSaveString;
  int iCentrality, iCentralityMatched;
  
  // Read the reco/gen histograms from the file and fit them to construct the closure plots
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(auto centralityBin : analyzedCentralityBin){
      iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
      iCentralityMatched = cardVector.at(iFile)->FindBinIndexCentrality(centralityBin);
      for(int iParton : partonIndices){
        if(drawPtClosure){
          for(int iGenJetPt = minGenPt; iGenJetPt < JetBackgroundHistogramManager::knGenJetPtBins; iGenJetPt++){

            // Read the reco/gen histogram from the file
            hRecoGenRatio[iFile][iGenJetPt][iCentrality][iParton] = closureHistograms.at(iFile)->GetHistogramJetPtClosure(iGenJetPt, JetBackgroundHistogramManager::knJetEtaBins, JetBackgroundHistogramManager::knJetPhiBins, iCentralityMatched, iParton);
          
            // Fit a gauss to the histogram
            if(drawGaussFitsPt){
              genPtString = Form("%d < Gen p_{T} < %d", 50+10*iGenJetPt, 60+10*iGenJetPt);
              centralityString = Form("Cent: %d-%d", centralityBin.first, centralityBin.second);
              jetTypeString = Form("%s jets", jetTypeName[iParton].Data());
              gaussFitSaveString = Form("figures/jetPtClosureGaussFit%sJets_T%dC%d.pdf", jetTypeName[iParton].Data(), iGenJetPt, iCentrality);
              std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatio[iFile][iGenJetPt][iCentrality][iParton], " ", jetTypeString, centralityString, genPtString, gaussFitSaveString);
            } else {
              std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatio[iFile][iGenJetPt][iCentrality][iParton]);
            }
          
            // Fill the histogram with the fit parameters
            hJetPtClosure[iFile][iCentrality][iParton]->SetBinContent(iGenJetPt+1,gaussMean);
            hJetPtClosure[iFile][iCentrality][iParton]->SetBinError(iGenJetPt+1,gaussMeanError);
            hJetPtClosureSigma[iFile][iCentrality][iParton]->SetBinContent(iGenJetPt+1,gaussSigma);
            hJetPtClosureSigma[iFile][iCentrality][iParton]->SetBinError(iGenJetPt+1,gaussSigmaError);
          
          } // Generator level jet pT loop
        } // pT closure if
      
        if(drawEtaClosure){
          // For eta, bins from 9 to nBins-9 cover the region -1.6 < eta < 1.6
          for(int iJetEta = 9; iJetEta < JetBackgroundHistogramManager::knJetEtaBins-9; iJetEta++){
          
            // Read the reco/gen histogram from the file
            hRecoGenRatioEta[iFile][iJetEta][iCentrality][iParton] = closureHistograms.at(iFile)->GetHistogramJetPtClosure(JetBackgroundHistogramManager::knGenJetPtBins, iJetEta, JetBackgroundHistogramManager::knJetPhiBins, iCentralityMatched, iParton);
          
            // Fit a gauss to the histogram
            std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatioEta[iFile][iJetEta][iCentrality][iParton]);
          
            // Fill the histogram with the fit parameters
            hJetPtClosureEta[iFile][iCentrality][iParton]->SetBinContent(iJetEta+1,gaussMean);
            hJetPtClosureEta[iFile][iCentrality][iParton]->SetBinError(iJetEta+1,gaussMeanError);
            hJetPtClosureSigmaEta[iFile][iCentrality][iParton]->SetBinContent(iJetEta+1,gaussSigma);
            hJetPtClosureSigmaEta[iFile][iCentrality][iParton]->SetBinError(iJetEta+1,gaussSigmaError);
          
          } // Jet eta loop
        } // eta closure if

        if(drawPhiClosure){

          // Loop over all jet phi bins
          for(int iJetPhi = 0; iJetPhi < JetBackgroundHistogramManager::knJetPhiBins; iJetPhi++){

            // In PbPb, due to detector inefficiencies, we cut the region -0.1 < phi < 1.2
            if((iJetPhi > 30) && (iJetPhi <46)) continue;
          
            // Read the reco/gen histogram from the file
            hRecoGenRatioPhi[iFile][iJetPhi][iCentrality][iParton] = closureHistograms.at(iFile)->GetHistogramJetPtClosure(JetBackgroundHistogramManager::knGenJetPtBins, JetBackgroundHistogramManager::knJetEtaBins, iJetPhi, iCentralityMatched, iParton);
          
            // Fit a gauss to the histogram
            std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatioPhi[iFile][iJetPhi][iCentrality][iParton]);
          
            // Fill the histogram with the fit parameters
            hJetPtClosurePhi[iFile][iCentrality][iParton]->SetBinContent(iJetPhi+1,gaussMean);
            hJetPtClosurePhi[iFile][iCentrality][iParton]->SetBinError(iJetPhi+1,gaussMeanError);
            hJetPtClosureSigmaPhi[iFile][iCentrality][iParton]->SetBinContent(iJetPhi+1,gaussSigma);
            hJetPtClosureSigmaPhi[iFile][iCentrality][iParton]->SetBinError(iJetPhi+1,gaussSigmaError);
          
          } // Jet eta loop
        } // eta closure if
      
      } // Closure particle loop (quark/gluon/no selection)
    } // Centrality loop
  } // File loop
  
  double minFitPt = 50+10*minGenPt;
  double maxFitPt = 500;
  
  // Fit the resolution plots with a polynomial function
  if(fitResolution){
    for(int iFile = 0; iFile < nFiles; iFile++){
      for(auto centralityBin : analyzedCentralityBin){
        iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);
        hJetPtClosureSigma[iFile][iCentrality][JetBackgroundHistograms::knInitialPartonTypes]->Fit("pol4","","",minFitPt,maxFitPt);
      } // Centrality loop
    } // File loop
  } // Fitting the resolution
  
  // Setup the legends for the closure plots
  std::vector<TH1D*> drawnHistograms;
  std::vector<TString> histogramLegend;
  TString legendTitle;

  // The include quark-gluon flag is automatically set to true if only one input file is provided
  if(includeQuarkGluon){
    // Setup legend for separating all jets, quark jet, and gluon jet components
    histogramLegend.push_back("Inclusive jets");
    histogramLegend.push_back("Quark jets");
    histogramLegend.push_back("Gluon jets");
    legendTitle = jetLegendString.at(0);
  } else {
    // If there are several files, draw only the all jets histogram from each file
    for(int iFile = 0; iFile < nFiles; iFile++){
      histogramLegend.push_back(jetLegendString.at(iFile));
    }
    legendTitle = "Inclusive jet";
  }

  // Draw the closure plots
  for(auto centralityBin : analyzedCentralityBin){
    iCentrality = cardVector.at(0)->FindBinIndexCentrality(centralityBin);

    if(drawPtClosure){

      // Set up the drawn histograms jet energy scale as a function of pT
      drawnHistograms.clear();

      for(int iFile = 0; iFile < nFiles; iFile++){
        for(int iParton : partonIndices){
          drawnHistograms.push_back(hJetPtClosure[iFile][iCentrality][iParton]);
        }
      }
      drawClosureHistogram(drawnHistograms, cardVector.at(0), "Gen p_{T} (GeV)", "#mu(reco p_{T} / gen p_{T})", iCentrality, 0, legendTitle, histogramLegend, Form("PtClosure_%s", saveComment.Data()), saveFigures);

      // Set up the drawn histograms jet energy resolution as a function of pT
      drawnHistograms.clear();
      for(int iFile = 0; iFile < nFiles; iFile++){
        for(int iParton : partonIndices){
          drawnHistograms.push_back(hJetPtClosureSigma[iFile][iCentrality][iParton]);
        }
      }
      drawClosureHistogram(drawnHistograms, cardVector.at(0), "Gen p_{T} (GeV)", "#sigma(reco p_{T} / gen p_{T})", iCentrality, 1, legendTitle, histogramLegend, Form("PtResolution_%s", saveComment.Data()), saveFigures);
    }
    
    if(drawEtaClosure){

      // Set up the drawn histograms jet energy scale as a function of eta
      drawnHistograms.clear();
      for(int iFile = 0; iFile < nFiles; iFile++){
        for(int iParton : partonIndices){
          drawnHistograms.push_back(hJetPtClosureEta[iFile][iCentrality][iParton]);
        }
      }
      drawClosureHistogram(drawnHistograms, cardVector.at(0), "#eta", "#mu(reco p_{T} / gen p_{T})", iCentrality, 0, legendTitle, histogramLegend, Form("EtaClosure_%s", saveComment.Data()), saveFigures);
      
      // Set up the drawn histograms jet energy resolution as a function of eta
      drawnHistograms.clear();
      for(int iFile = 0; iFile < nFiles; iFile++){
        for(int iParton : partonIndices){
          drawnHistograms.push_back(hJetPtClosureSigmaEta[iFile][iCentrality][iParton]);
        }
      }
      drawClosureHistogram(drawnHistograms, cardVector.at(0), "#eta", "#sigma(reco p_{T} / gen p_{T})", iCentrality, 1, legendTitle, histogramLegend, Form("EtaResolution_%s", saveComment.Data()), saveFigures);
      
    }

    if(drawPhiClosure){

      // Set up the drawn histograms jet energy scale as a function of phi
      drawnHistograms.clear();
      for(int iFile = 0; iFile < nFiles; iFile++){
        for(int iParton : partonIndices){
          drawnHistograms.push_back(hJetPtClosurePhi[iFile][iCentrality][iParton]);
        }
      }
      drawClosureHistogram(drawnHistograms, cardVector.at(0), "#varphi", "#mu(reco p_{T} / gen p_{T})", iCentrality, 0, legendTitle, histogramLegend, Form("PhiClosure_%s", saveComment.Data()), saveFigures);
      
      // Set up the drawn histograms jet energy resolution as a function of phi
      drawnHistograms.clear();
      for(int iFile = 0; iFile < nFiles; iFile++){
        for(int iParton : partonIndices){
          drawnHistograms.push_back(hJetPtClosureSigmaPhi[iFile][iCentrality][iParton]);
        }
      }
      drawClosureHistogram(drawnHistograms, cardVector.at(0), "#varphi", "#sigma(reco p_{T} / gen p_{T})", iCentrality, 1, legendTitle, histogramLegend, Form("PhiResolution_%s", saveComment.Data()), saveFigures);
    }
    
  } // Centrality loop
  
}
