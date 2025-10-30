#include "JDrawer.h"

/*
 *  Implementation of the flow function
 */
double flowFunction(double* x, double* par) {
  unsigned int nFlow = par[0];          // Number of fitted flow components is defined in the first parameter
  unsigned int firstFittedVn = par[1];  // The first fitted flow component is defined in the second parameter

  // Add each component separately to the total fit value
  double flowModulation = par[2];
  for (unsigned int iFlow = 0; iFlow < nFlow; iFlow++) {
    flowModulation += par[2] * 2.0 * par[2 * iFlow + 3] * std::cos((iFlow + firstFittedVn) * (x[0] - par[2 * iFlow + 4]));
  }
  return flowModulation;
}

/*
 * Macro for plotting example flow histogram and fit to it
 */
void plotFlowFitExample(){

	// Configuration
	bool saveFigures = false;
	TString figureComment = "_central_lowProbability";
	int maxPlots = 20;
	double minCentrality = 0;
	double maxCentrality = 10;
	double minProbability = 0.4;
	double maxProbability = 0.6;

	// Open the file with example histograms
	TFile* inputFile = TFile::Open("data/HiForestMiniAOD_defaultFlow_withDebugVectors_42.root");

	// Read the tree containing the flow fit debug vectors
	TTree* jetTree = (TTree*) inputFile->Get("akFlowPuCs4PFJetAnalyzer/t");
	TTree* eventTree = (TTree*) inputFile->Get("hiEvtAnalyzer/HiTree");

	// By default, disable the branches in the trees. We are only interested in a small subset
	jetTree->SetBranchStatus("*",0);
	eventTree->SetBranchStatus("*",0);

	// From the event tree, we only need centrality information
	Int_t hiBin; TBranch* hiBinBranch;
	eventTree->SetBranchStatus("hiBin",1);
	eventTree->SetBranchAddress("hiBin", &hiBin, &hiBinBranch);

	// From the jet tree, we want to read the vectors containing flow fit debug information
	vector<double>* flowFitParameters = new vector<double>(); flowFitParameters->clear();
	TBranch* flowFitParametersBranch;
  jetTree->SetBranchStatus("flowFitParameters", 1);
  jetTree->SetBranchAddress("flowFitParameters", &flowFitParameters, &flowFitParametersBranch);

  vector<double>* flowFitDebugInfo = new vector<double>(); flowFitDebugInfo->clear();
  TBranch* flowDebugInfoBranch;
  jetTree->SetBranchStatus("flowFitDebugInfo", 1);
  jetTree->SetBranchAddress("flowFitDebugInfo", &flowFitDebugInfo, &flowDebugInfoBranch);

  Int_t nEvents = flowFitParametersBranch->GetEntries();

	// Variables needed to hold the event variables
	double centrality;
	int firstFittedVn, lastFittedVn, nFlow;
	int chi2index;
	double flowFitProbability;

  int nFlowFitBins;
	TH1D* flowFitHistogram = new TH1D("flowFitHistogram", "flowFitHistogram", 10, -TMath::Pi(), TMath::Pi());
  flowFitHistogram->Sumw2();

  // For the fit function, we need to define the number of parameters. These can be read from the first event
  // Also read all other variables that do not change event to event
  jetTree->GetEntry(0);
  
  firstFittedVn = flowFitParameters->back(); // The last entry in the flow parameter vector gives the first fitted flow component
  lastFittedVn = firstFittedVn + ((flowFitParameters->size() - 4) / 2) - 1; // The last fitted component can be calculated from the size of the fir parameter vector
  chi2index = flowFitParameters->size() - 3; // Index in the vector that holds fit chi2 value
  nFlow = lastFittedVn - firstFittedVn + 1;

  TF1* flowFitFunction = new TF1("flowFitFunction", flowFunction, -TMath::Pi(), TMath::Pi(), nFlow * 2 + 3);
  flowFitFunction->FixParameter(0, nFlow);           // The first parameter defines the number of fitted flow components
  flowFitFunction->FixParameter(1, firstFittedVn);   // The second parameter defines the first fitted flow component

  // Define helper variables for drawing
  JDrawer* drawer = new JDrawer();
  TLegend* legend;

  // Loop variables
  int iPlot = 0;

  // Event loop
  for(int iEvent = 0; iEvent < nEvents; iEvent++){

	  // Get the first event from the trees
	  eventTree->GetEntry(iEvent);
	  jetTree->GetEntry(iEvent);

    // Decode the information read from the trees
	  centrality = hiBin / 2.0;

	  if(centrality < minCentrality) continue;
	  if(centrality >= maxCentrality) continue;

    // The flow fit probability is calculated using the following formula
    flowFitProbability = ROOT::Math::chisquared_cdf_c(flowFitParameters->at(chi2index), flowFitParameters->at(chi2index + 1));

    if(flowFitProbability < minProbability) continue;
    if(flowFitProbability >= maxProbability) continue;

    // The histogram used to do the flow fit can be accessed through the flow fit debug info vector
    nFlowFitBins = flowFitDebugInfo->size()-1;
    flowFitHistogram->SetBins(nFlowFitBins, -TMath::Pi(), TMath::Pi());
    for(size_t iBin = 1; iBin <= nFlowFitBins; iBin++){
      flowFitHistogram->SetBinContent(iBin, flowFitDebugInfo->at(iBin));
      flowFitHistogram->SetBinError(iBin, TMath::Sqrt(flowFitDebugInfo->at(iBin)));
    }

    // The flow fit function can be reconstructed from the information available in the flow fit parameter vector
    flowFitFunction->SetParameter(2, flowFitParameters->at(0));
    for(Int_t iFlow = 0; iFlow < nFlow*2; iFlow++){
      flowFitFunction->SetParameter(iFlow+3, flowFitParameters->at(iFlow+1));
    }

	  // Draw the histogram and fit function to the same canvas
	  drawer->DrawHistogram(flowFitHistogram, "#varphi", "Counts", " ");
	  flowFitFunction->Draw("same");

	  legend = new TLegend(0.2,0.74,0.4,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*)0, Form("Cent: %.1f-%.1f%%", centrality, centrality+0.5), "");
    legend->AddEntry((TObject*)0, Form("Fit prob: %.4f", flowFitProbability), "");
    legend->Draw();

    // Option to save the figures 
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/pfCandidateFlowFit%s_event%d.pdf", figureComment.Data(), iEvent));
    }

    // If we have reached the maximum number of plots, break the event loop
    if(iPlot++ >= maxPlots) break;
  }

}
