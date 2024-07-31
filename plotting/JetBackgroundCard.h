#ifndef JETBACKGROUNDCARD_H
#define JETBACKGROUNDCARD_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TObjString.h>
#include <TVectorT.h>

/*
 * JetBackgroundCard class
 *
 * This class reads the ConfigurationCard from the input root file and decodes the
 * necessary information for the analysis.
 */
class JetBackgroundCard {
  
public:
 
  // Indices for card entries
  enum enumCardEntries{
    kMaxParticleEtaEventPlane,    // Maximum eta for particles included in the event plane calculation
    kMaxParticlePtEventPlane,     // Maximum pT for particles included in the event plane calculation
    kMatchJets,                   // 0 = Reco and Gen jets are not matched, 1 = They are matched, 2 = Anti-matching
    kJetType,                     // 0 = Reconstructed jets, 1 = Generator level
    kJetSubtraction,              // 0 = Calo PU jets, 1 = csPF jets, 2 = flowPuCsPF jets
    kJetAxis,                     // 0 = E-scheme axis, 1 = WTA axis
    kJetEtaCut,                   // Eta cut for jets
    kMinPtCut,                    // Minimum allowed pT for the inclusive jets
    kMaxPtCut,                    // Maximum allowed pT for the inclusive jets
    kMinMaxTrackPtFraction,       // Minimum fraction of jet pT taken by the highest pT track in jet
    kMaxMaxTrackPtFraction,       // Maximum fraction of jet pT taken by the highest pT track in jet
    kMinJetPtClosure,             // Minimum generator level jet pT for closure histograms
    kZVertexCut,                  // Maximum accepted vz in the event
    kLowPtHatCut,                 // Minimum accepted pT hat
    kHighPtHatCut,                // Maximum accepted pT hat
    kCentralityBinEdges,          // Centrality bin edges
    kJetPtBinEdges,               // Jet pT bin edges
    kTrackPtBinEdges,             // Track pT bin edges
    kPtHatBinEdges,               // pT hat bin edges
    knEntries};                   // Number of entries in the card
  
private:
  
  // Names for each entry read from the configuration card
  const char* fCardEntryNames[knEntries] = {"MaxParticleEtaEventPlane", "MaxParticlePtEventPlane", "MatchJets", "JetType", "JetSubtraction", "JetAxis", "JetEtaCut", "MinJetPtCut", "MaxJetPtCut", "MinMaxTrackPtFraction", "MaxMaxTrackPtFraction", "MinJetPtClosure", "ZVertexCut", "LowPtHatCut", "HighPtHatCut", "CentralityBinEdges", "JetPtBinEdges", "TrackPtBinEdges", "PtHatBinEdges"};

  const char* fInputFileSaveName = "InputFile";
  
  TFile* fInputFile;         // Input file from which all the data is read
  TString fCardDirectory;    // Path to the ConfigurationCard directory
  
  void ReadVectors();        // Read the vectors from the file
  
  // Strings for git hash
  TObjString* fGitHash;
  TObjString* fProjectionGitHash;
  
  // Vectors for all the lines inside the card
  TVectorT<float>* fCardEntries[knEntries];   // Array of all the vectors in the card
  TObjString* fInputFileName;                 // Input file name from which the histogram are projected
  
  // Private methods
  int GetNBins(const int index) const;                                // Get the number of bins for internal index
  double GetLowBinBorder(const int index, const int iBin) const;      // Get the low border of i:th bin from internal index
  double GetHighBinBorder(const int index, const int iBin) const;     // Get the high border of i:th bin from internal index
  int GetBinIndex(const int index, const double value) const;         // Get the bin index in the i:th bin from internal index based on given value
  int FindBinIndex(const int index, const double lowBinBorder, const double highBinBorder) const; // Get the bin index in the i:th bin from internal index based on the provided bin borders
   
public:
  
  JetBackgroundCard(TFile* inFile); // Contructor with input file
  ~JetBackgroundCard() = default;   // Destructor
  
  void Write(TDirectory* file);            // Write the contents of the card to a file
  void Print() const;                      // Print the contents of the card to the console
  
  int GetNCentralityBins() const; // Get the number of centrality bins
  int GetNTrackPtBins() const;    // Get the number of track pT bins
  int GetNJetPtBins() const;   // Get the number of jet pT bins
  double GetLowBinBorderCentrality(const int iBin) const;  // Get the low border of i:th centrality bin
  double GetLowBinBorderTrackPt(const int iBin) const;     // Get the low border of i:th track pT bin
  double GetLowBinBorderJetPt(const int iBin) const;    // Get the low border of i:th jet pT bin
  double GetHighBinBorderCentrality(const int iBin) const; // Get the high border of i:th centrality bin
  double GetHighBinBorderTrackPt(const int iBin) const;    // Get the high border of i:th track pT bin
  double GetHighBinBorderJetPt(const int iBin) const;   // Get the high border of i:th jet pT bin
  std::pair<double,double> GetBinBordersCentrality(const int iBin) const; // Get the bin borders of the i:th centrality bin
  std::pair<double,double> GetBinBordersTrackPt(const int iBin) const; // Get the bin borders of the i:th track pT bin
  std::pair<double,double> GetBinBordersJetPt(const int iBin) const; // Get the bin borders of the i:th jet pT bin
  int GetBinIndexCentrality(const double value) const;     // Get the bin index for a given centrality value
  int GetBinIndexTrackPt(const double value) const;        // Get the bin index for a given track pT value
  int GetBinIndexJetPt(const double value) const;          // Get the bin index for a given jet pT value
  int FindBinIndexCentrality(const double lowBorder, const double highBorder) const; // Find if a centrality bin with given borders exists and return its index
  int FindBinIndexCentrality(const std::pair<double,double> binBorders) const; // Find if a centrality bin with given borders exists and return its index
  int FindBinIndexTrackPt(const double lowBorder, const double highBorder) const;    // Find if a track pT bin with given borders exists and return its index
  int FindBinIndexTrackPt(const std::pair<double,double> binBorders) const;    // Find if a track pT bin with given borders exists and return its index
  int FindBinIndexJetPt(const double lowBorder, const double highBorder) const;   // Find if a jet pT bin with given borders exists and return its index
  int FindBinIndexJetPt(const std::pair<double,double> binBorders) const;   // Find if a jet pT bin with given borders exists and return its index
  int GetJetType() const;          // Get the jet type index
  int GetJetSubtraction() const;   // Get the background subtraction method used for jets
  double GetJetPtCut() const;      // Get the minimum jet pT cut
  
  void AddOneDimensionalVector(int entryIndex, float entryContent); // Add one dimensional vector to the card
  void AddVector(int entryIndex, int dimension, double* contents); // Add a vector to the card
  void AddFileName(TString fileName); // Add a file name to the card
  void AddProjectionGitHash(const char* gitHash); // Add a git hash used to project the histograms to the file
  
};

#endif
