#include "JDrawer.h"

void plotHistogram(){

  // Read the input file
  TFile* inputFile = TFile::Open("eventPlaneCorrelation/jetBackgroundHistograms_optimizedCutSet1_2024-08-09.root");

  // Read the desired histogram from the input file
  //TString histogramName = "flowFitParameters/flowFitProbability_C0";
  TString histogramName = "inclusiveJet/inclusiveJetPhi_C0";
  TH1D* myHistogram = (TH1D*) inputFile->Get(histogramName);

  // Draw the histogram
  JDrawer* drawer = new JDrawer();
  TLegend* legend;

  TString xAxisTitle = "Inclusive jet #varphi";
  TString yAxisTitle = "A.U.";

  myHistogram->SetMarkerStyle(kFullCircle);
  drawer->DrawHistogram(myHistogram, xAxisTitle, yAxisTitle, " ");

  legend = new TLegend(0.2,0.65,0.4,0.84);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*)0, "Run2 MC (2018)", "");
  legend->AddEntry((TObject*)0, "Cent: 4-14%", "");
  legend->AddEntry((TObject*)0, "Jet p_{T} > 80 GeV", "");
  legend->Draw();

}
