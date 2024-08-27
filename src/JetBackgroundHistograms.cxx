// Class for histograms needed in the jet background analysis

// C++ includes
#include <assert.h>

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "JetBackgroundHistograms.h"

/*
 * Default constructor
 */
JetBackgroundHistograms::JetBackgroundHistograms() :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhEvents(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhInclusiveJet(0),
  fhLeadingJet(0),
  fhCalorimeterJet(0),
  fhJetPtClosure(0),
  fCard(0)
{
  // Default constructor

  for(int iEventPlane = 0; iEventPlane < knEventPlanes; iEventPlane++){
    fhInclusiveJetEventPlane[iEventPlane] = NULL;
    fhLeadingJetEventPlane[iEventPlane] = NULL;
    fhCalorimeterJetEventPlane[iEventPlane] = NULL;
  }
  
}

/*
 * Custom constructor
 */
JetBackgroundHistograms::JetBackgroundHistograms(ConfigurationCard* newCard) :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhEvents(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhInclusiveJet(0),
  fhLeadingJet(0),
  fhCalorimeterJet(0),
  fhJetPtClosure(0),
  fCard(newCard)
{
  // Custom constructor

  for(int iEventPlane = 0; iEventPlane < knEventPlanes; iEventPlane++){
    fhInclusiveJetEventPlane[iEventPlane] = NULL;
    fhLeadingJetEventPlane[iEventPlane] = NULL;
    fhCalorimeterJetEventPlane[iEventPlane] = NULL;
  }
}

/*
 * Copy constructor
 */
JetBackgroundHistograms::JetBackgroundHistograms(const JetBackgroundHistograms& in) :
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhEvents(in.fhEvents),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted),
  fhPtHat(in.fhPtHat),
  fhPtHatWeighted(in.fhPtHatWeighted),
  fhInclusiveJet(in.fhInclusiveJet),
  fhLeadingJet(in.fhLeadingJet),
  fhCalorimeterJet(in.fhCalorimeterJet),
  fhJetPtClosure(in.fhJetPtClosure),
  fCard(in.fCard)
{
  // Copy constructor

  for(int iEventPlane = 0; iEventPlane < knEventPlanes; iEventPlane++){
    fhInclusiveJetEventPlane[iEventPlane] = in.fhInclusiveJetEventPlane[iEventPlane];
    fhLeadingJetEventPlane[iEventPlane] = in.fhLeadingJetEventPlane[iEventPlane];
    fhCalorimeterJetEventPlane[iEventPlane] = in.fhCalorimeterJetEventPlane[iEventPlane];
  }

}

/*
 * Assingment operator
 */
JetBackgroundHistograms& JetBackgroundHistograms::operator=(const JetBackgroundHistograms& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fhVertexZ = in.fhVertexZ;
  fhVertexZWeighted = in.fhVertexZWeighted;
  fhEvents = in.fhEvents;
  fhCentrality = in.fhCentrality;
  fhCentralityWeighted = in.fhCentralityWeighted;
  fhPtHat = in.fhPtHat;
  fhPtHatWeighted = in.fhPtHatWeighted;
  fhInclusiveJet = in.fhInclusiveJet;
  fhLeadingJet = in.fhLeadingJet;
  fhCalorimeterJet = in.fhCalorimeterJet;
  fhJetPtClosure = in.fhJetPtClosure;
  fCard = in.fCard;

  for(int iEventPlane = 0; iEventPlane < knEventPlanes; iEventPlane++){
    fhInclusiveJetEventPlane[iEventPlane] = in.fhInclusiveJetEventPlane[iEventPlane];
    fhLeadingJetEventPlane[iEventPlane] = in.fhLeadingJetEventPlane[iEventPlane];
    fhCalorimeterJetEventPlane[iEventPlane] = in.fhCalorimeterJetEventPlane[iEventPlane];
  }
  
  return *this;
}

/*
 * Destructor
 */
JetBackgroundHistograms::~JetBackgroundHistograms(){
  // destructor
  delete fhVertexZ;
  delete fhVertexZWeighted;
  delete fhEvents;
  delete fhCentrality;
  delete fhCentralityWeighted;
  delete fhPtHat;
  delete fhPtHatWeighted;
  delete fhInclusiveJet;
  delete fhLeadingJet;
  delete fhCalorimeterJet;
  delete fhJetPtClosure;

  for(int iEventPlane = 0; iEventPlane < knEventPlanes; iEventPlane++){
    delete fhInclusiveJetEventPlane[iEventPlane];
    delete fhLeadingJetEventPlane[iEventPlane];
    delete fhCalorimeterJetEventPlane[iEventPlane];
  }
}

/*
 * Set the configuration card used for the histogram class
 */
void JetBackgroundHistograms::SetCard(ConfigurationCard* newCard){
  fCard = newCard;
}

/*
 * Create the necessary histograms
 */
void JetBackgroundHistograms::CreateHistograms(){
  
  // ======== Common binning information for histograms =========
  
  // Centrality
  const Double_t minCentrality = -1;   // Minimum centrality bin is negative since hiBin is -0.5 for pp
  const Double_t maxCentrality = 100;  // Maximum centrality bin
  const Int_t nCentralityBins = 202;      // Number of centrality bins
  
  // Jet pT
  const Double_t minPtJet = 0;     // Minimum jet pT
  const Double_t maxPtJet = 500;   // Maximum jet pT
  const Int_t nPtBinsJet = 100;    // Number of jet pT bins
  
  // Phi
  const Double_t minPhi = -TMath::Pi();  // Minimum phi
  const Double_t maxPhi = TMath::Pi();   // Maximum phi
  const Int_t nPhiBins = 64;             // Number of phi bins
  
  // Eta
  const Double_t minEta = -2.5;    // Minimum eta (current eta cut for tracks = 2.4)
  const Double_t maxEta = 2.5;     // Maximum eta (current eta cut for tracks = 2.4)
  const Int_t nEtaBins = 50;       // Number of eta bins
  
  // Vertex z-position
  const Double_t minVz = -20;   // Minimum vz
  const Double_t maxVz = 20;    // Maximum vz
  const Int_t nVzBins = 80;     // Number of vz bins
  
  // pT hat
  const Double_t minPtHat = 0;     // Minimum pT hat
  const Double_t maxPtHat = 460;   // Maximum pT hat
  const Int_t nFinePtHatBins = 230; // Number of fine pT hat bins
  
  // Generator level pT binning for closure histograms
  const Double_t minClosurePt = 50;                             // Minimum gen jet pT for closure plots
  const Double_t maxClosurePt = 500;                            // Maximum gen jet pT for closure plots
  const Int_t nClosurePtBins = (maxClosurePt-minClosurePt)/10;  // Bin width of 10 for the Gen pT in closure plots
  
  // Particle type for closure plots (0 = quark, 1 = gluon, 2 = undetermined)
  const Double_t minClosureParticleType = -0.5;                        // Closure particle type indexing starts from zero
  const Double_t maxClosureParticleType = knInitialPartonTypes-0.5;  // Maximum closure particle type index
  const Int_t nClosureParticleTypeBins = knInitialPartonTypes;       // Bin width for particle type is 1

  // Flag for existance of a matching jet
  const Double_t minMatchingJetFlag = -0.5;
  const Double_t maxMatchingJetFlag = knMatchingTypes-0.5;
  const Int_t nMatchingJetFlags = knMatchingTypes;
  
  // Binning for reco/gen ratio for closure histograms
  const Double_t minClosureRatio = 0.025;    // Minimum ratio for the closure plots
  const Double_t maxClosureRatio = 2.025;    // Maximum ratio for the closure plots
  const Int_t nClosureRatioBins = 40;    // Number of closure ratio bins

  // DeltaPhi in [-pi/2,3pi/2]
  const Double_t minDeltaPhiJetEventPlane = -TMath::Pi()/2.0;    // Minimum deltaPhi for jet-event plane correlations
  const Double_t maxDeltaPhiJetEventPlane = 3.0*TMath::Pi()/2.0; // Maximum deltaPhi for jet-event plane correlations
  const Int_t nDeltaPhiBinsJetEventPlane = 200;                  // Number of deltaPhi bins for jet-event plane correlations
  
  // Centrality bins for THnSparses (We run into memory issues, if have all the bins)
  const Int_t nWideCentralityBins = fCard->GetNBin("CentralityBinEdges");
  Double_t wideCentralityBins[nWideCentralityBins+1];
  for(Int_t iCentrality = 0; iCentrality < nWideCentralityBins+1; iCentrality++){
    wideCentralityBins[iCentrality] = fCard->Get("CentralityBinEdges",iCentrality);
  }
  
  // Bins for the pT hat histogram
  const Int_t nPtHatBins = fCard->GetNBin("PtHatBinEdges");
  Double_t ptHatBins[nPtHatBins+1];
  for(Int_t iPtHat = 0; iPtHat < nPtHatBins+1; iPtHat++){
    ptHatBins[iPtHat] = fCard->Get("PtHatBinEdges",iPtHat);
  }
  
  // Jet pT binning for event plane correlation histograms
  const Int_t nJetPtBinsEventPlane = fCard->GetNBin("JetPtBinEdges");
  Double_t jetPtBinsEventPlane[nJetPtBinsEventPlane+1];
  for(Int_t iJetPt = 0; iJetPt < nJetPtBinsEventPlane+1; iJetPt++){
    jetPtBinsEventPlane[iJetPt] = fCard->Get("JetPtBinEdges",iJetPt);
  }
  const Double_t minJetPtEventPlane = jetPtBinsEventPlane[0];
  const Double_t maxJetPtEventPlane = jetPtBinsEventPlane[nJetPtBinsEventPlane];

  // Arrays for creating THnSparses
  const Int_t nAxesJet = 6;
  Int_t nBinsJet[nAxesJet];
  Double_t lowBinBorderJet[nAxesJet];
  Double_t highBinBorderJet[nAxesJet];
  
  const Int_t nAxesJetClosure = 7;
  Int_t nBinsJetClosure[nAxesJetClosure];
  Double_t lowBinBorderJetClosure[nAxesJetClosure];
  Double_t highBinBorderJetClosure[nAxesJetClosure];

  const Int_t nAxesJetEventPlaneCorrelation = 3;
  Int_t nBinsJetPtEventPlaneCorrelation[nAxesJetEventPlaneCorrelation];
  Double_t lowBinBorderJetEventPlaneCorrelation[nAxesJetEventPlaneCorrelation];
  Double_t highBinBorderJetEventPlaneCorrelation[nAxesJetEventPlaneCorrelation];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1F("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhVertexZWeighted = new TH1F("vertexZweighted","vertexZweighted",nVzBins,minVz,maxVz); fhVertexZWeighted->Sumw2();
  fhEvents = new TH1F("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhCentrality = new TH1F("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  fhCentralityWeighted = new TH1F("centralityWeighted","centralityWeighted",nCentralityBins,minCentrality,maxCentrality); fhCentralityWeighted->Sumw2();
  fhPtHat = new TH1F("pthat","pthat",nPtHatBins,ptHatBins); fhPtHat->Sumw2();
  fhPtHatWeighted = new TH1F("pthatWeighted","pthatWeighted",nFinePtHatBins,minPtHat,maxPtHat); fhPtHatWeighted->Sumw2();
  
  // For the event histogram, label each bin corresponding to an event cut
  for(Int_t i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // ======== THnSparse for all jets ========
  
  // Axis 0 for the any jet histogram: jet pT
  nBinsJet[0] = nPtBinsJet;         // nBins for any jet pT
  lowBinBorderJet[0] = minPtJet;    // low bin border for any jet pT
  highBinBorderJet[0] = maxPtJet;   // high bin border for any jet pT
  
  // Axis 1 for the any jet histogram: jet phi
  nBinsJet[1] = nPhiBins;        // nBins for any jet phi
  lowBinBorderJet[1] = minPhi;   // low bin border for any jet phi
  highBinBorderJet[1] = maxPhi;  // high bin border for any jet phi
  
  // Axis 2 for the any jet histogram: jet eta
  nBinsJet[2] = nEtaBins;        // nBins for any jet eta
  lowBinBorderJet[2] = minEta;   // low bin border for any jet eta
  highBinBorderJet[2] = maxEta;  // high bin border for any jet eta
  
  // Axis 3 for the any jet histogram: centrality
  nBinsJet[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderJet[3] = minCentrality;  // low bin border for centrality
  highBinBorderJet[3] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the jet histogram: jet flavor (quark/gluon/undetermined)
  nBinsJet[4] = nClosureParticleTypeBins;        // nBins for jet flavor
  lowBinBorderJet[4] = minClosureParticleType;   // low bin border for jet flavor
  highBinBorderJet[4] = maxClosureParticleType;  // high bin border for jet flavor

  // Axis 5 for the jet histogram: matching jet flag (doesn't have/has matching jet)
  nBinsJet[5] = nMatchingJetFlags;            // nBins for matching jet flag
  lowBinBorderJet[5] = minMatchingJetFlag;    // low bin border for matching jet flag
  highBinBorderJet[5] = maxMatchingJetFlag;   // high bin border for matching jet flag
  
  // Create the histogram for all jets using the above binning information
  fhInclusiveJet = new THnSparseF("inclusiveJet", "inclusiveJet", nAxesJet, nBinsJet, lowBinBorderJet, highBinBorderJet); fhInclusiveJet->Sumw2();
  fhLeadingJet = new THnSparseF("leadingJet", "leadingJet", nAxesJet, nBinsJet, lowBinBorderJet, highBinBorderJet); fhLeadingJet->Sumw2();
  fhCalorimeterJet = new THnSparseF("calorimeterJet", "calorimeterJet", nAxesJet, nBinsJet, lowBinBorderJet, highBinBorderJet); fhCalorimeterJet->Sumw2();

  // Set custom centrality bins for histograms
  fhInclusiveJet->SetBinEdges(3,wideCentralityBins);
  fhLeadingJet->SetBinEdges(3,wideCentralityBins);
  fhCalorimeterJet->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparses for jet pT closures ========
  
  // Axis 0 for the jet pT closure histogram: generator level jet pT
  nBinsJetClosure[0] = nClosurePtBins;       // nBins for generator level pT bins in closure plots
  lowBinBorderJetClosure[0] = minClosurePt;  // low bin border generator level pT in closure plots
  highBinBorderJetClosure[0] = maxClosurePt; // high bin border generator level pT in closure plots
  
  // Axis 1 for the jet pT closure histogram: reconstructed jet pT
  nBinsJetClosure[1] = nClosurePtBins;       // nBins for reconstructed jet pT bins in closure plots
  lowBinBorderJetClosure[1] = minClosurePt;  // low bin border for reconstructed jet pT in closure plots
  highBinBorderJetClosure[1] = maxClosurePt; // high bin border for reconstructed jet pT in closure plots
  
  // Axis 2 for the jet pT closure histogram: generator level jet eta
  nBinsJetClosure[2] = nEtaBins;             // nBins for jet eta
  lowBinBorderJetClosure[2] = minEta;        // low bin border for jet eta
  highBinBorderJetClosure[2] = maxEta;       // high bin border for jet eta
  
  // Axis 3 for the jet pT closure histogram: centrality
  nBinsJetClosure[3] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderJetClosure[3] = minCentrality;    // low bin border for centrality
  highBinBorderJetClosure[3] = maxCentrality;   // high bin border for centrality
  
  // Axis 4 for the jet pT closure histogram: ref parton = quark/gluon/undetermined
  nBinsJetClosure[4] = nClosureParticleTypeBins;         // nBins for reference parton
  lowBinBorderJetClosure[4] = minClosureParticleType;    // low bin border for reference parton
  highBinBorderJetClosure[4] = maxClosureParticleType;   // high bin border for reference parton
  
  // Axis 5 for the jet pT closure histogram: reco/gen ratio for closure
  nBinsJetClosure[5] = nClosureRatioBins;        // nBins for closure ratio
  lowBinBorderJetClosure[5] = minClosureRatio;   // low bin border for closure ratio
  highBinBorderJetClosure[5] = maxClosureRatio;  // high bin border for closure ratio

  // Axis 6 for the jet pT closure histogram: generator level jet phi
  nBinsJetClosure[6] = nPhiBins;             // nBins for jet phi
  lowBinBorderJetClosure[6] = minPhi;        // low bin border for jet phi
  highBinBorderJetClosure[6] = maxPhi;       // high bin border for jet phi
  
  // Create histograms for jet pT closure
  fhJetPtClosure = new THnSparseF("jetPtClosure", "jetPtClosure", nAxesJetClosure, nBinsJetClosure, lowBinBorderJetClosure, highBinBorderJetClosure); fhJetPtClosure->Sumw2();
  
  // Set custom centrality bins for histograms
  fhJetPtClosure->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparses for jet-event plane correlation study ========
  
  // Axis 0 for the additional histogram: DeltaPhi between jet and event plane angle
  nBinsJetPtEventPlaneCorrelation[0] = nDeltaPhiBinsJetEventPlane;       // nBins for deltaPhi between jet and event plane
  lowBinBorderJetEventPlaneCorrelation[0] = minDeltaPhiJetEventPlane;    // low bin border for deltaPhi between jet and event plane
  highBinBorderJetEventPlaneCorrelation[0] = maxDeltaPhiJetEventPlane;   // high bin border for deltaPhi between jet and event plane
  
  // Axis 1 for the additional histogram: jet pT
  nBinsJetPtEventPlaneCorrelation[1] = nJetPtBinsEventPlane;       // nBins for wide multiplicity bins
  lowBinBorderJetEventPlaneCorrelation[1] = minJetPtEventPlane;    // low bin border for wide multiplicity
  highBinBorderJetEventPlaneCorrelation[1] = maxJetPtEventPlane;   // high bin border for wide multiplicity
  
  // Axis 2 for the additional histogram: centrality
  nBinsJetPtEventPlaneCorrelation[2] = nWideCentralityBins;           // nBins for centrality
  lowBinBorderJetEventPlaneCorrelation[2] = minCentrality;          // low bin border for centrality
  highBinBorderJetEventPlaneCorrelation[2] = maxCentrality;         // high bin border for centrality
  
  // Create histograms for event plane study
  for(int iEventPlane = 0; iEventPlane < knEventPlanes; iEventPlane++){
    fhInclusiveJetEventPlane[iEventPlane] = new THnSparseF(Form("inclusiveJetEventPlaneOrder%d", iEventPlane+2), Form("inclusiveJetEventPlaneOrder%d", iEventPlane+2), nAxesJetEventPlaneCorrelation, nBinsJetPtEventPlaneCorrelation, lowBinBorderJetEventPlaneCorrelation, highBinBorderJetEventPlaneCorrelation); fhInclusiveJetEventPlane[iEventPlane]->Sumw2();
    fhLeadingJetEventPlane[iEventPlane] = new THnSparseF(Form("leadingJetEventPlaneOrder%d", iEventPlane+2), Form("leadingJetEventPlaneOrder%d", iEventPlane+2), nAxesJetEventPlaneCorrelation, nBinsJetPtEventPlaneCorrelation, lowBinBorderJetEventPlaneCorrelation, highBinBorderJetEventPlaneCorrelation); fhLeadingJetEventPlane[iEventPlane]->Sumw2();
    fhCalorimeterJetEventPlane[iEventPlane] = new THnSparseF(Form("calorimeterJetEventPlaneOrder%d", iEventPlane+2), Form("calorimeterJetEventPlaneOrder%d", iEventPlane+2), nAxesJetEventPlaneCorrelation, nBinsJetPtEventPlaneCorrelation, lowBinBorderJetEventPlaneCorrelation, highBinBorderJetEventPlaneCorrelation); fhCalorimeterJetEventPlane[iEventPlane]->Sumw2();
    
    // Set custom centrality bins for histograms
    fhInclusiveJetEventPlane[iEventPlane]->SetBinEdges(1,jetPtBinsEventPlane);
    fhLeadingJetEventPlane[iEventPlane]->SetBinEdges(1,jetPtBinsEventPlane);
    fhCalorimeterJetEventPlane[iEventPlane]->SetBinEdges(1,jetPtBinsEventPlane);
    fhInclusiveJetEventPlane[iEventPlane]->SetBinEdges(2,wideCentralityBins);
    fhLeadingJetEventPlane[iEventPlane]->SetBinEdges(2,wideCentralityBins);
    fhCalorimeterJetEventPlane[iEventPlane]->SetBinEdges(2,wideCentralityBins);
  }

}

/*
 * Write the histograms to file
 */
void JetBackgroundHistograms::Write() const{
  
  // Write the histograms to file
  fhVertexZ->Write();
  fhVertexZWeighted->Write();
  fhEvents->Write();
  fhCentrality->Write();
  fhCentralityWeighted->Write();
  fhPtHat->Write();
  fhPtHatWeighted->Write();
  fhInclusiveJet->Write();
  fhLeadingJet->Write();
  fhCalorimeterJet->Write();
  fhJetPtClosure->Write();

  for(int iEventPlane = 0; iEventPlane < knEventPlanes; iEventPlane++){
    fhInclusiveJetEventPlane[iEventPlane]->Write();
    fhLeadingJetEventPlane[iEventPlane]->Write();
    fhCalorimeterJetEventPlane[iEventPlane]->Write();
  }
}

/*
 * Write the histograms to a given file
 */
void JetBackgroundHistograms::Write(TString outputFileName) const{
  
  // Define the output file
  TFile* outputFile = new TFile(outputFileName, "RECREATE");
  
  // Write the histograms to file
  Write();
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
}