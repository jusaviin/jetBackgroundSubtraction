// Reader for jet and particle trees from CMS Monte Carlo simulation
//
//===========================================================
// MonteCarloForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef MONTECARLOFORESTREADER_H
#define MONTECARLOFORESTREADER_H

// C++ includes
#include <iostream>
#include <assert.h>
#include <vector>

// Root includes
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>

// Own includes
#include "JetBackgroundHistograms.h"

using namespace std;

class MonteCarloForestReader{
  
private:
  static const Int_t fnMaxJet = 250;        // Maximum number of jets in an event
  
public:
  
  // Possible data types to be read with the reader class
  enum enumJetType {kReconstructedJet, kGeneratorLevelJet, knJetTypes};
  
  // Constructors and destructors
  MonteCarloForestReader();                                              // Default constructor
  MonteCarloForestReader(Int_t jetType, Int_t jetAxis);                  // Custom constructor
  MonteCarloForestReader(const MonteCarloForestReader& in);              // Copy constructor
  ~MonteCarloForestReader();                                             // Destructor
  MonteCarloForestReader& operator=(const MonteCarloForestReader& obj);  // Equal sign operator
  
  // Methods
  void GetEvent(Int_t iEvent);                 // Get the i:th event in tree
  Int_t GetNEvents() const;                    // Get the number of events
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  void ReadForestFromFileList(std::vector<TString> fileList);   // Read the forest from a file list
  void BurnForest();                           // Burn the forest
  
  // Getters for leaves in heavy ion tree
  Float_t GetVz() const;              // Getter for vertex z position
  Float_t GetCentrality() const;      // Getter for centrality
  Int_t GetHiBin() const;             // Getter for CMS hiBin
  Float_t GetPtHat() const;           // Getter for pT hat
  Float_t GetEventWeight() const;     // Getter for event weight in MC

  // Getters for leaves in skim tree
  Int_t GetPrimaryVertexFilterBit() const;           // Getter for primary vertex filter bit
  Int_t GetHfCoincidenceFilterBit() const;           // Getter for hadronic forward coincidence filter bit
  Int_t GetClusterCompatibilityFilterBit() const;    // Getter for cluster compatibility filter bit
  
  // Getters for leaves in the defined jet tree
  Int_t GetNJets(Int_t jetType) const;                     // Getter for number of jets
  Float_t GetJetPt(Int_t jetType, Int_t iJet) const;       // Getter for jet pT
  Float_t GetJetPhi(Int_t jetType, Int_t iJet) const;        // Getter for jet phi
  Float_t GetJetEta(Int_t jetType, Int_t iJet) const;        // Getter for jet eta
  Float_t GetJetRawPt(Int_t jetType, Int_t iJet) const;      // Getter for jet raw pT
  Float_t GetJetMaxTrackPt(Int_t jetType, Int_t iJet) const; // Getter for maximum track pT inside a jet

  // Matching between gen and reco
  Bool_t HasMatchingJet(Int_t jetType, Int_t iJet) const;  // Check if there is matching reconstructed/generator level jet
  Int_t GetMatchingIndex(Int_t jetType, Int_t iJet) const; // Get the index of the matched jet
  Float_t GetMatchedPt(Int_t jetType, Int_t iJet) const;   // Getter for matched jet pT
  Float_t GetMatchedEta(Int_t jetType, Int_t iJet) const;  // Getter for matched jet eta
  Float_t GetMatchedPhi(Int_t jetType, Int_t iJet) const;  // Getter for matched jet phi

  // Getters for leaves in reconstructed jet tree
  Int_t GetNJets() const;                     // Getter for number of jets
  Float_t GetJetPt(Int_t iJet) const;         // Getter for jet pT
  Float_t GetJetPhi(Int_t iJet) const;        // Getter for jet phi
  Float_t GetJetEta(Int_t iJet) const;        // Getter for jet eta
  Float_t GetJetRawPt(Int_t iJet) const;      // Getter for jet raw pT
  Float_t GetJetMaxTrackPt(Int_t iJet) const; // Getter for maximum track pT inside a jet
  
  // Reco to gen matching
  Bool_t HasMatchingGenJet(Int_t iJet) const; // Check if reconstructed jet has a matching generator level jet
  Int_t GetMatchingGenIndex(Int_t iJet) const; // Get the matching generator level jet index for the given reconstructed jet
  Float_t GetMatchedGenPt(Int_t iJet) const;  // Getter for matched generator level jet pT
  Float_t GetMatchedGenEta(Int_t iJet) const; // Getter for matched generator level jet eta
  Float_t GetMatchedGenPhi(Int_t iJet) const; // Getter for matched generator level jet phi

  // Getters for leaves in generator level jet tree
  Int_t GetNGeneratorJets() const;                   // Getter for number of generator level jets
  Float_t GetGeneratorJetPt(Int_t iJet) const;       // Getter for generator level jet pT
  Float_t GetGeneratorJetPhi(Int_t iJet) const;      // Getter for generator level jet phi
  Float_t GetGeneratorJetEta(Int_t iJet) const;      // Getter for generator level jet eta

  // Gen to reco matching
  Bool_t HasMatchingRecoJet(Int_t iJet) const;  // Check if generator level has a matching reconstructed jet
  Int_t GetMatchingRecoIndex(Int_t iJet) const; // Get the matching reconstructed jet index for the given generator level jet
  Float_t GetMatchedRecoPt(Int_t iJet) const;   // Getter for matched reconstructed jet pT
  Float_t GetMatchedRecoEta(Int_t iJet) const;  // Getter for matched reconstructed jet eta
  Float_t GetMatchedRecoPhi(Int_t iJet) const;  // Getter for matched reconstructed jet phi

  // Jet flavor
  Int_t GetJetFlavor(Int_t jetType, Int_t iJet) const; // Getter for the jet flavor for input jet type
  Int_t GetRecoJetFlavor(Int_t iJet) const;            // Getter for reconstructed jet flavor
  Int_t GetGenJetFlavor(Int_t iJet) const;             // Getter for generator level jet flavor
  
  // Getters for leaves in track tree
  Int_t GetNTracks() const;                                  // Getter for number of tracks
  Float_t GetTrackPt(Int_t iTrack) const;                    // Getter for track pT
  Float_t GetTrackPtError(Int_t iTrack) const;               // Getter for track pT error
  Float_t GetTrackPhi(Int_t iTrack) const;                   // Getter for track phi
  Float_t GetTrackEta(Int_t iTrack) const;                   // Getter for track eta
  Bool_t GetTrackHighPurity(Int_t iTrack) const;             // Getter for the high purity of the track
  Float_t GetTrackVertexDistanceZ(Int_t iTrack) const;       // Getter for track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceZError(Int_t iTrack) const;  // Getter for error of track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceXY(Int_t iTrack) const;      // Getter for track distance from primary vertex in xy-direction
  Float_t GetTrackVertexDistanceXYError(Int_t iTrack) const; // Getter for error of track distance from primary vertex in xy-direction
  Float_t GetTrackNormalizedChi2(Int_t iTrack) const;        // Getter for normalized track chi2 value from reconstruction fit
  Int_t GetNHitsTrackerLayer(Int_t iTrack) const;            // Getter for number of hits in tracker layers
  Int_t GetNHitsTrack(Int_t iTrack) const;                   // Getter for number of hits for the track
  Float_t GetTrackEnergyEcal(Int_t iTrack) const;            // Getter for track energy in ECal
  Float_t GetTrackEnergyHcal(Int_t iTrack) const;            // Getter for track energy in HCal
  Int_t GetTrackCharge(Int_t iTrack) const;                  // Getter for track charge
  
  // Getters for leaves in generator level particle tree
  Int_t GetNGenParticles() const;                            // Getter for number of generator level particles
  Float_t GetGenParticlePt(Int_t iTrack) const;              // Getter for generator level particle pT
  Float_t GetGenParticlePhi(Int_t iTrack) const;             // Getter for generator level particle phi
  Float_t GetGenParticleEta(Int_t iTrack) const;             // Getter for generator level particle eta
  Int_t GetGenParticleCharge(Int_t iTrack) const;            // Getter for generator level particle charge
  Int_t GetGenParticleSubevent(Int_t iTrack) const;          // Getter for generator level particle subevent index
  
private:
  
  // Methods
  void Initialize();      // Connect the branches to the tree
    
  Int_t fJetType;         // Choose the type of jets used for analysis. 0 = Calo PU jets, 1 = PF CS jets, 2 = Flow subtracted Pf CS jets
  Int_t fJetAxis;         // Jet axis used for the jets. 0 = Anti-kT, 1 = WTA
  
  // Trees in the forest
  TTree* fHeavyIonTree;    // Tree for heavy ion event information
  TTree* fSkimTree;        // Tree for event cuts
  TTree* fJetTree;         // Tree for jet information
  TTree* fTrackTree;       // Tree for reconstructed tracks
  TTree* fGenParticleTree; // Tree for generator level particles

  // Branches for heavy ion tree
  TBranch* fHiVzBranch;                   // Branch for vertex z-position
  TBranch* fHiBinBranch;                  // Branch for centrality
  TBranch* fPtHatBranch;                  // Branch for pT hat
  TBranch* fEventWeightBranch;            // Branch for event weight

  // Branches for skim tree
  TBranch* fPrimaryVertexBranch;          // Branch for primary vertex filter bit
  TBranch* fHfCoincidenceBranch;          // Branch for energy recorded in HF calorimeter towers
  TBranch* fClusterCompatibilityBranch;   // Branch for cluster compatibility
  
  // Branches for jet tree
  TBranch* fnJetsBranch;          // Branch for number of jets
  TBranch* fJetPtBranch;          // Branch for jet pT
  TBranch* fJetPhiBranch;         // Branch for jet phi
  TBranch* fJetWTAPhiBranch;      // Branch for jet phi
  TBranch* fJetEtaBranch;         // Branch for jet eta
  TBranch* fJetWTAEtaBranch;      // Branch for jet eta
  TBranch* fJetRawPtBranch;       // Branch for raw jet pT
  TBranch* fJetMaxTrackPtBranch;  // Maximum pT for a track inside a jet

  TBranch* fJetRefPtBranch;      // Branch for reference generator level pT for a reconstructed jet
  TBranch* fJetRefEtaBranch;     // Branch for reference generator level eta for a reconstructed jet
  TBranch* fJetRefPhiBranch;     // Branch for reference generator level phi for a reconstructed jet
  TBranch* fJetRefFlavorBranch;  // Branch for reference generator level jet flavor

  TBranch* fnGenJetsBranch;      // Branch for number of generator level jets
  TBranch* fGenJetPtBranch;      // Branch for generator level jet pT
  TBranch* fGenJetPhiBranch;     // Branch for generator level jet phi
  TBranch* fGenJetWTAPhiBranch;  // Branch for generator level jet WTA phi
  TBranch* fGenJetEtaBranch;     // Branch for generator level jet eta
  TBranch* fGenJetWTAEtaBranch;  // Branch for generator level jet eta

  // Branches for track tree
  TBranch* fnTracksBranch;                     // Branch for number of tracks
  TBranch* fTrackPtBranch;                     // Branch for track pT
  TBranch* fTrackPtErrorBranch;                // Branch for track pT error
  TBranch* fTrackPhiBranch;                    // Branch for track phi
  TBranch* fTrackEtaBranch;                    // Branch for track eta
  TBranch* fHighPurityTrackBranch;             // Branch for high purity of the track
  TBranch* fTrackVertexDistanceZBranch;        // Branch for track distance from primary vertex in z-direction
  TBranch* fTrackVertexDistanceZErrorBranch;   // Branch for error for track distance from primary vertex in z-direction
  TBranch* fTrackVertexDistanceXYBranch;       // Branch for track distance from primary vertex in xy-direction
  TBranch* fTrackVertexDistanceXYErrorBranch;  // Branch for error for track distance from primary vertex in xy-direction
  TBranch* fTrackChi2Branch;                   // Branch for track chi2 value from reconstruction fit
  TBranch* fnTrackDegreesOfFreedomBranch;      // Branch for number of degrees of freedom in reconstruction fit
  TBranch* fnHitsTrackerLayerBranch;           // Branch for number of hits in tracker layers
  TBranch* fnHitsTrackBranch;                  // Branch for number of hits for the track
  TBranch* fTrackEnergyEcalBranch;             // Branch for track energy in ECal
  TBranch* fTrackEnergyHcalBranch;             // Branch for track energy in HCal
  TBranch* fTrackChargeBranch;                 // Branch for track charge

  // Branches for genenerator level particle tree
  TBranch* fGenParticlePtBranch;        // Branch for generator level particle pT:s
  TBranch* fGenParticlePhiBranch;       // Branch for generator level particle phis
  TBranch* fGenParticleEtaBranch;       // Branch for generator level particle etas
  TBranch* fGenParticleChargeBranch;    // Branch for generator level particle charges
  TBranch* fGenParticleSubeventBranch;  // Branch for generator level particle subevent indices (0 = PYTHIA, (>0) = HYDJET)

  // Leaves for heavy ion tree
  Float_t fVertexZ;    // Vertex z-position
  Int_t fHiBin;        // HiBin = Centrality percentile * 2
  Float_t fPtHat;      // pT hat

  // Leaves for the skim tree
  Int_t fPrimaryVertexFilterBit;           // Filter bit for primary vertex
  Int_t fHfCoincidenceFilterBit;           // Filter bit for energy recorded in HF calorimeter towers
  Int_t fClusterCompatibilityFilterBit;    // Filter bit for cluster compatibility
  
  // Leaves for jet tree
  Int_t fnJets;          // number of jets in an event
  Int_t fnGenJets;       // Number of generator level jets in an event
  Float_t fEventWeight;  // jet weight in the MC tree
  
  Float_t fJetPtArray[fnMaxJet] = {0};         // pT:s of all the jets in an event
  Float_t fJetPhiArray[fnMaxJet] = {0};        // phis of all the jets in an event
  Float_t fJetWTAPhiArray[fnMaxJet] = {0};     // WTA phis of all the jets in an event
  Float_t fJetEtaArray[fnMaxJet] = {0};        // etas of all the jets in an event
  Float_t fJetWTAEtaArray[fnMaxJet] = {0};     // WTA etas of all the jets in an event
  Float_t fJetRawPtArray[fnMaxJet] = {0};      // raw jet pT for all the jets in an event
  Float_t fJetMaxTrackPtArray[fnMaxJet] = {0}; // maximum track pT inside a jet for all the jets in an event
   
  Float_t fJetRefPtArray[fnMaxJet] = {0};      // reference generator level pT for a reconstructed jet
  Float_t fJetRefEtaArray[fnMaxJet] = {0};     // reference generator level eta for a reconstructed jet
  Float_t fJetRefPhiArray[fnMaxJet] = {0};     // reference generator level phi for a reconstructed jet
  Int_t fJetRefFlavorArray[fnMaxJet] = {0};    // flavor for initiating parton for the reference gen jet

  Float_t fGenJetPtArray[fnMaxJet] = {0};       // pT:s of the generator level jets in an event
  Float_t fGenJetPhiArray[fnMaxJet] = {0};      // phis of the generator level jets in an event
  Float_t fGenJetWTAPhiArray[fnMaxJet] = {0};   // WTA phis of the generator level jets in an event
  Float_t fGenJetEtaArray[fnMaxJet] = {0};      // etas of the generator level jets in an event
  Float_t fGenJetWTAEtaArray[fnMaxJet] = {0};   // WTA etas of the generator level jets in an event
  
  // Leaves for the track tree regardless of forest type
  Int_t fnTracks;  // Number of tracks

  // Leaves for the track tree
  vector<float>* fTrackPtVector;                     // Vector for track pT:s
  vector<float>* fTrackPtErrorVector;                // Vector for track pT errors
  vector<float>* fTrackPhiVector;                    // Vector for track phis
  vector<float>* fTrackEtaVector;                    // Vector for track etas
  vector<bool>* fHighPurityTrackVector;              // Vector for the high purity of tracks
  vector<float>* fTrackVertexDistanceZVector;        // Vector for track distance from primary vertex in z-direction
  vector<float>* fTrackVertexDistanceZErrorVector;   // Vector for error for track distance from primary vertex in z-direction
  vector<float>* fTrackVertexDistanceXYVector;       // Vector for track distance from primary vertex in xy-direction
  vector<float>* fTrackVertexDistanceXYErrorVector;  // Vector for error for track distance from primary vertex in xy-direction
  vector<float>* fTrackNormalizedChi2Vector;         // Vector for normalized track chi2 value from reconstruction fit
  vector<char>* fnHitsTrackerLayerVector;            // Vector for number of hits in tracker layers
  vector<char>* fnHitsTrackVector;                   // Vector for number of hits for the track
  vector<float>* fTrackEnergyEcalVector;             // Vector for track energy in ECal
  vector<float>* fTrackEnergyHcalVector;             // Vector for track energy in HCal
  vector<char>* fTrackChargeVector;                  // Vector for track charges

  // Leaves for the generator level particle tree
  Int_t fnGenParticles;                    // Number of generator level particles
  vector<float>* fGenParticlePtArray;      // Array for generator level particle pT:s
  vector<float>* fGenParticlePhiArray;     // Array for generator level particle phis
  vector<float>* fGenParticleEtaArray;     // Array for generator level particle etas
  vector<int>* fGenParticleChargeArray;    // Array for generator level particle charges
  vector<int>* fGenParticleSubeventArray;  // Array for generator level particle subevent indices (0 = PYTHIA, (>0) = HYDJET)
};

#endif
