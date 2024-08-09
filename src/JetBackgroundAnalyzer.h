// Class for the main analysis for jet background subtraction

#ifndef JETBACKGROUNDANALYZER_H
#define JETBACKGROUNDANALYZER_H

// C++ includes
#include <vector>
#include <bitset>
#include <assert.h>   // Standard c++ debugging tool. Terminates the program if expression given evaluates to 0.
#include <tuple>      // For returning several arguments in a transparent manner
#include <fstream>
#include <string>

// Root includes
#include <TString.h>
#include <TRandom3.h>
#include <TMath.h>

// Own includes
#include "ConfigurationCard.h"
#include "JetBackgroundHistograms.h"
#include "MonteCarloForestReader.h"
#include "JetCorrector.h"
#include "JetUncertainty.h"
#include "JetMetScalingFactorManager.h"

class JetBackgroundAnalyzer{
  
private:
  
  enum enumFilledHistograms{kFillEventInformation, kFillJets, kFillTracks, kFillJetConeHistograms, kFillEnergyEnergyCorrelators, kFillEnergyEnergyCorrelatorsSystematics, kFillJetPtClosure, kFillJetPtUnfoldingResponse, kFillTrackParticleMatchingHistograms, knFillTypes}; // Histograms to fill
  
public:
  
  // Constructors and destructor
  JetBackgroundAnalyzer(); // Default constructor
  JetBackgroundAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard* newCard); // Custom constructor
  JetBackgroundAnalyzer(const JetBackgroundAnalyzer& in); // Copy constructor
  virtual ~JetBackgroundAnalyzer(); // Destructor
  JetBackgroundAnalyzer& operator=(const JetBackgroundAnalyzer& obj); // Equal sign operator
  
  // Methods
  void RunAnalysis();                     // Run the dijet analysis
  JetBackgroundHistograms* GetHistograms() const;   // Getter for histograms

 private:
  
  // Private methods
  void ReadConfigurationFromCard(); // Read all the configuration from the input card
  
  Bool_t PassEventCuts(MonteCarloForestReader* eventReader, const Bool_t fillHistograms); // Check if the event passes the event cuts
  Double_t GetVzWeight(const Double_t vz) const;  // Get the proper vz weighting depending on analyzed system
  Double_t GetCentralityWeight(const Int_t hiBin) const; // Get the proper centrality weighting depending on analyzed system
  Double_t GetSmearingFactor(Double_t jetPt, Double_t jetEta, const Double_t centrality); // Getter for jet pT smearing factor
  Int_t GetCentralityBin(const Double_t centrality) const; // Getter for centrality bin
  Double_t GetDeltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2) const; // Get deltaR between two objects
  
  // Private data members
  MonteCarloForestReader* fEventReader;            // Reader for jets in the event
  std::vector<TString> fFileNames;               // Vector for all the files to loop over
  ConfigurationCard* fCard;                      // Configuration card for the analysis
  JetBackgroundHistograms* fHistograms;                    // Filled histograms
  TF1* fVzWeightFunction;                        // Weighting function for vz. Needed for MC.
  TF1* fCentralityWeightFunctionCentral;         // Weighting function for central centrality classes. Needed for MC.
  TF1* fCentralityWeightFunctionPeripheral;      // Weighting function for peripheral centrality classes. Needed for MC.
  TF1* fSmearingFunction;                        // Additional smearing for jets. Needed in systematic uncertainty study.
  JetCorrector* fJetCorrector2018;               // Class for making jet energy correction for 2018 data
  JetMetScalingFactorManager* fEnergyResolutionSmearingFinder; // Manager to find proper jet energy resolution scaling factors provided by the JetMet group
  TRandom3* fRng;                                // Random number generator
  
  // Analyzed data and forest types
  Int_t fJetType;                    // Type of jets used for analysis. 0 = Reconstructed jets, 1 = Generator level jets
  Int_t fJetSubtraction;             // Background subtraction algorithm. 0 = Calo jets with PU, 1 = PF jets with CS, 2 = PF jets with flow CS 
  Int_t fDebugLevel;                 // Amount of debug messages printed to console
  
  // Weights for filling the MC histograms
  Double_t fVzWeight;                // Weight for vz in MC
  Double_t fCentralityWeight;        // Weight for centrality in MC
  Double_t fPtHatWeight;             // Weight for pT hat in MC
  Double_t fTotalEventWeight;        // Combined weight factor for MC

  // Event plane calculation cuts
  Double_t fMaxParticleEtaEventPlane;  // Maximum eta value for particles used to determine the event plane
  Double_t fMaxParticlePtEventPlane;   // Maximum pT value for particles used to determine the event plane
  
  // Jet selection cuts
  Int_t fJetAxis;                      // Used jet axis type. 0 = Anti-kT jet axis, 1 = Axis from leading PF candidate
  Double_t fVzCut;                     // Cut for vertez z-position in an event
  Double_t fMinimumPtHat;              // Minimum accepted pT hat value
  Double_t fMaximumPtHat;              // Maximum accepted pT hat value
  Double_t fJetEtaCut;                 // Eta cut around midrapidity
  Double_t fJetMinimumPtCut;           // Minimum pT cut for jets
  Double_t fJetMaximumPtCut;           // Maximum pT accepted for jets (and tracks)
  Double_t fMinimumMaxTrackPtFraction; // Cut for jets consisting only from soft particles
  Double_t fMaximumMaxTrackPtFraction; // Cut for jets consisting only from one high pT
  Double_t fJetClosureMinimumPt;       // Minimum jet pT for jet pT closure plots
  
  // Jet pT closure histogram filling is optional
  Bool_t fFillJetPtClosure;            // Fill jet pT closure histograms

};

#endif
