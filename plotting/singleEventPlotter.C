#include "AlgorithmLibrary.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Plot particle phi distribution with reconstructed and generator level jet information from single events
 *
 *  Arguments:
 *    TString inputFileName = Input root file
 *    TString logFileName = Text log for running the analysis
 */
void singleEventPlotter(TString inputFileName, TString logFileName){
  
  // Extract interesting information from the log file:
  std::ifstream file_stream(logFileName);
  std::string line;
  TObjArray* lineContents;
  TObjString* lineItem;
  TString itemString;

  double eventPlaneAngle;
  std::vector<double> recoJetPhi;
  std::vector<double> recoJetPt;
  std::vector<double> genJetPhi;
  std::vector<double> genJetPt;

  double flowFitV1 = 0;
  double flowFitEventPlane1 = 0;
  double flowFitV2 = 0;
  double flowFitEventPlane2 = 0;
  double flowFitV3 = 0;
  double flowFitEventPlane3 = 0;
  double flowFitV4 = 0;
  double flowFitEventPlane4 = 0;
  double flowFitQuality = 0;
  double flowFitAmplitude = 0;
  int iEvent = 0;

  // Open the input file list for reading
  if( file_stream.is_open() ) {
    
    // Loop over the lines in the input file list
    while( !file_stream.eof() ) {
      getline(file_stream, line);
      TString lineString(line);
      
      // If the line is non-empty, check if it contains a coded message and interpret it
      if( lineString.CompareTo("", TString::kExact) != 0 ) {

        // Tokenize the line and search for coded messages from the first part
        lineContents = lineString.Tokenize(" "); // Tokenize the string from ' ' character
        lineItem = (TObjString*)lineContents->At(0);
        itemString = lineItem->String().Strip(TString::kBoth, ' ');

        // If the first part of the string is EVENTPLANE, the second part gives the second order event plane angle
        if(itemString.CompareTo("EVENTPLANE", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          eventPlaneAngle = lineItem->String().Atof();
        }

        // Of the first part of the string is RECOJET, read the phi and pT of the reconstructed jet
        if(itemString.CompareTo("RECOJET", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          recoJetPhi.push_back(lineItem->String().Atof());
          lineItem = (TObjString*)lineContents->At(2);
          recoJetPt.push_back(lineItem->String().Atof());
        }
        
        // Of the first part of the string is GENJET, read the phi and pT of the generator level jet
        if(itemString.CompareTo("GENJET", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          genJetPhi.push_back(lineItem->String().Atof());
          lineItem = (TObjString*)lineContents->At(2);
          genJetPt.push_back(lineItem->String().Atof());
        }

        // If the first part of the string is FLOWFITV1, the second part gives the v1 component in the flow fit
        if(itemString.CompareTo("FLOWFITV1", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          flowFitV1 = lineItem->String().Atof();
        }

        // If the first part of the string is FLOWFITV2, the second part gives the v2 component in the flow fit
        if(itemString.CompareTo("FLOWFITV2", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          flowFitV2 = lineItem->String().Atof();
        }

        // If the first part of the string is FLOWFITV3, the second part gives the v3 component in the flow fit
        if(itemString.CompareTo("FLOWFITV3", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          flowFitV3 = lineItem->String().Atof();
        }

        // If the first part of the string is FLOWFITV4, the second part gives the v4 component in the flow fit
        if(itemString.CompareTo("FLOWFITV4", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          flowFitV4 = lineItem->String().Atof();
        }

        // If the first part of the string is FLOWFITEVENTPLANE1, the second part gives the order 1 event plane in the flow fit
        if(itemString.CompareTo("FLOWFITEVENTPLANE1", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          flowFitEventPlane1 = lineItem->String().Atof();
        }

        // If the first part of the string is FLOWFITEVENTPLANE2, the second part gives the order 2 event plane in the flow fit
        if(itemString.CompareTo("FLOWFITEVENTPLANE2", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          flowFitEventPlane2 = lineItem->String().Atof();
        }

        // If the first part of the string is FLOWFITEVENTPLANE3, the second part gives the order 3 event plane in the flow fit
        if(itemString.CompareTo("FLOWFITEVENTPLANE3", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          flowFitEventPlane3 = lineItem->String().Atof();
        }

        // If the first part of the string is FLOWFITEVENTPLANE4, the second part gives the order 4 event plane in the flow fit
        if(itemString.CompareTo("FLOWFITEVENTPLANE4", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          flowFitEventPlane4 = lineItem->String().Atof();
        }


        // If the first part of the string is FLOWFITQUALITY, the second part gives the quality measure of the flow fit
        if(itemString.CompareTo("FLOWFITQUALITY", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          flowFitQuality = lineItem->String().Atof();
        }

        // If the first part of the string is FLOWFITAMPLITUDE, the second part gives the amplitude of the flow fit
        if(itemString.CompareTo("FLOWFITAMPLITUDE", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(1);
          flowFitAmplitude = lineItem->String().Atof();
        }

        // If the first part of the string is Analyzing, the third part gives the event number
        if(itemString.CompareTo("Analyzing", TString::kExact) == 0){
          lineItem = (TObjString*)lineContents->At(2);
          iEvent = lineItem->String().Atoi();
        }

      } // Empty line if
      
    } // Loop over lines in the file
    
  // If cannot read the file, give error and end program
  } else {
  	cout << "ERROR! Could not open the file " << logFileName.Data() << endl;
  	return;
  }

  // AlgorithmLibrary for histogram min-maxing
  AlgorithmLibrary* minMaxer = new AlgorithmLibrary();
  std::pair<double,double> currentMinMax = std::make_pair(10e10,-10e10);
  
  // Open the input file
  TFile* inputFile = TFile::Open(inputFileName);

  // Read the histograms for phi and pT sums as a function of phi from Pythia and Hydjet
  TH1D* hPhiPythia = (TH1D*) inputFile->Get("phiPythia");
  TH1D* hPhiHydjet = (TH1D*) inputFile->Get("phiHydjet");
  TH1D* hPtSumPythia = (TH1D*) inputFile->Get("phiPythiaPt");
  TH1D* hPtSumHydjet = (TH1D*) inputFile->Get("phiHydjetPt");
  TH1D* hPhiPfCandidate = (TH1D*) inputFile->Get("phiPfCandidate");
  TH1D* hPtSumPfCandidate = (TH1D*) inputFile->Get("phiPfCandidatePt");

  // Determine average level in the hydjet phi distribution
  double averageAmplitude = hPhiHydjet->Integral("width") / TMath::Pi() / 2.0;

  // Create a function that is used in the flow fit from the extracted input
  TF1* flowFit = new TF1("flowFit", "[0] * (1 + 2*[1]*cos((x-[2])) + 2*[3]*cos(2*(x-[4])) + 2*[5]*cos(3*(x-[6])) + 2*[7]*cos(4*(x-[8])))", -TMath::Pi(), TMath::Pi());
  flowFit->SetParameters(averageAmplitude, flowFitV1, flowFitEventPlane1, flowFitV2, flowFitEventPlane2, flowFitV3, flowFitEventPlane3, flowFitV4, flowFitEventPlane4);

  TF1* flowFitPfCand = new TF1("flowFitPfCand", "[0] * (1 + 2*[1]*cos((x-[2])) + 2*[3]*cos(2*(x-[4])) + 2*[5]*cos(3*(x-[6])) + 2*[7]*cos(4*(x-[8])))", -TMath::Pi(), TMath::Pi());
  flowFitPfCand->SetParameters(flowFitAmplitude, flowFitV1, flowFitEventPlane1, flowFitV2, flowFitEventPlane2, flowFitV3, flowFitEventPlane3, flowFitV4, flowFitEventPlane4);
  flowFitPfCand->SetLineColor(kGreen+3);

  // Change line colors for histograms
  hPhiPythia->SetLineColor(kBlue);
  hPhiHydjet->SetLineColor(kRed);
  hPhiPfCandidate->SetLineColor(kGreen+3);
  hPtSumPythia->SetLineColor(kBlue);
  hPtSumHydjet->SetLineColor(kRed);
  hPtSumPfCandidate->SetLineColor(kGreen+3);

  // Draw the histograms to the same canvas
  JDrawer* drawer = new JDrawer();
  drawer->DrawHistogram(hPtSumPythia, "#phi", "#Sigma p_{T}", " ");
  hPtSumHydjet->Draw("same");

  // Draw lines for all reconstructed and generator level jets
  TLine* recoDrawer = new TLine();
  recoDrawer->SetLineColor(kMagenta);
  recoDrawer->SetLineWidth(2);
  TLine* genDrawer = new TLine();
  genDrawer->SetLineStyle(2);
  genDrawer->SetLineWidth(2);
  genDrawer->SetLineColor(kCyan);
  TLine* eventPlaneDrawer = new TLine();
  eventPlaneDrawer->SetLineWidth(2);
  eventPlaneDrawer->SetLineColor(kGreen+3);
  TLine* fitEventPlaneDrawer = new TLine();
  fitEventPlaneDrawer->SetLineWidth(2);
  fitEventPlaneDrawer->SetLineColor(9);

  double maxValue = hPtSumPythia->GetMaximum();

  eventPlaneDrawer->DrawLine(eventPlaneAngle, 0, eventPlaneAngle, maxValue);
  for(double phi : recoJetPhi){
  	recoDrawer->DrawLine(phi, 0, phi, maxValue);
  }
  for(double phi : genJetPhi){
  	genDrawer->DrawLine(phi, 0, phi, maxValue);
  }

  // Add a legend to the top right corner
  TLegend* legend = new TLegend(0.2,0.56,0.4,0.9);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry(hPtSumPythia, "Pythia", "l");
  legend->AddEntry(hPtSumHydjet, "Hydjet", "l");
  legend->AddEntry(eventPlaneDrawer, "#Psi_{2}", "l");
  legend->AddEntry(recoDrawer, "Reco jets", "l");
  legend->AddEntry(genDrawer, "Gen jets", "l");
  legend->Draw();

  // Determine a good drawing range for the histograms
  currentMinMax = std::make_pair(1e10,-1e10);
  currentMinMax = minMaxer->FindHistogramMinMax(hPhiPythia, currentMinMax);
  currentMinMax = minMaxer->FindHistogramMinMax(hPhiHydjet, currentMinMax);

  hPhiPythia->GetYaxis()->SetRangeUser(0, currentMinMax.second + 0.15*currentMinMax.second);

  // Make another canvas for histograms without pT weights
  drawer->DrawHistogram(hPhiPythia, "#phi", "Counts", Form("Event %d, fit quality: %.5f", iEvent, flowFitQuality));
  hPhiHydjet->Draw("same");
  hPhiPfCandidate->Draw("same");

  // Draw the flow fit to the same canvas
  flowFit->Draw("same");
  flowFitPfCand->Draw("same");

  // Draw the jet location and determined event plane to the same figure
  maxValue = hPhiPythia->GetMaximum();

  eventPlaneDrawer->DrawLine(eventPlaneAngle, 0, eventPlaneAngle, maxValue);
  fitEventPlaneDrawer->DrawLine(flowFitEventPlane2, 0, flowFitEventPlane2, maxValue);
  //fitEventPlaneDrawer->DrawLine(flowFitEventPlane2+TMath::Pi(), 0, flowFitEventPlane2+TMath::Pi(), maxValue);
  for(double phi : recoJetPhi){
  	recoDrawer->DrawLine(phi, 0, phi, maxValue);
  }
  for(double phi : genJetPhi){
  	genDrawer->DrawLine(phi, 0, phi, maxValue);
  }

  // Add a legend to the figure
  legend = new TLegend(0.2,0.56,0.4,0.9);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry(hPhiPythia, "Pythia", "l");
  legend->AddEntry(hPhiHydjet, "Hydjet", "l");
  legend->AddEntry(hPhiPfCandidate, "PF candidates", "l");
  legend->AddEntry(eventPlaneDrawer, "#Psi_{2} (gen)", "l");
  legend->AddEntry(fitEventPlaneDrawer, "#Psi_{2} (PF cand)", "l");
  legend->AddEntry(recoDrawer, "Reco jets", "l");
  legend->AddEntry(genDrawer, "Gen jets", "l");

  legend->Draw();

}
