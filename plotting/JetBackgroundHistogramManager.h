#ifndef JETBACKGROUNDHISTOGRAMMANAGER_H
#define JETBACKGROUNDHISTOGRAMMANAGER_H

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>

// Own includes
#include "JetBackgroundCard.h"
#include "../src/JetBackgroundHistograms.h"
#include "AlgorithmLibrary.h"

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class JetBackgroundHistogramManager {

  // Possible data types to be read with the reader class
  enum enumJetType {kInclusiveJet, kLeadingJet, knJetTypes};
 
public:
  
  // Dimensions for histogram arrays
  static const int kMaxCentralityBins = 5;       // Maximum allowed number of centrality bins
  static const int kMaxJetPtBins = 20;           // Maximum allowed number of jet pT bins for closure histograms
  static const int knGenJetPtBins = 45;          // Number of generator level jet pT bins for jet pT closures
  static const int knJetEtaBins = 50;            // Number of jet eta bins for jet pT closures
  static const int knJetPhiBins = 64;            // Number of jet phi bins for jet pT closures
  
private:
  
  // Naming for jet histograms
  const char* fJetHistogramName[knJetTypes] = {"inclusiveJet", "leadingJet"};
  const char* fJetAxisName[knJetTypes] = {"Inclusive jet", "Leading jet"};
  
  // Naming for jet initiating particle
  const char* fInitialPartonName[JetBackgroundHistograms::knInitialPartonTypes+1] = {"_quark","_gluon","_undetermined",""};

  // Naming for matching
  const char* fMatchingName[JetBackgroundHistograms::knMatchingTypes+1] = {"_unmatched", "_matched", ""};
  
public:
  
  JetBackgroundHistogramManager();                                    // Default constructor
  JetBackgroundHistogramManager(TFile* inputFile);                    // Constructor with input file
  JetBackgroundHistogramManager(TFile* inputFile, JetBackgroundCard* card);     // Constructor with input file and card
  JetBackgroundHistogramManager(JetBackgroundCard* card);                       // Constructor with card
  JetBackgroundHistogramManager(const JetBackgroundHistogramManager& in);       // Copy constructor
  ~JetBackgroundHistogramManager();                                   // Destructor
  
  void LoadHistograms();          // Load the histograms from the inputfile
  void Write(const char* fileName, const char* fileOption);          // Write all the loaded histograms into a file
  void LoadProcessedHistograms(); // Load processed histograms from the inputfile
  
  // Setters for binning information
  void SetCentralityBins(const bool readBinsFromFile, const int nBins, const double* binBorders, bool setIndices = true); // Set up centrality bin indices according to provided bin borders
  void SetJetPtBins(const bool readBinsFromFile, const int nBins, const double* binBorders, bool setIndices = true);   // Set up jet pT bin indices for energy-energy correlator according to provided bin borders
  
  // Setters for event information and dijets
  void SetLoadEventInformation(const bool loadOrNot); // Setter for loading event information
  
  // Setters for jets
  void SetLoadJetHistograms(const bool loadOrNot);        // Setter for loading all jet histograms
  
  // Setter for loading additional histograms
  void SetLoad2DHistograms(const bool loadOrNot);           // Setter for loading two-dimensional histograms
  void SetLoadJetPtClosureHistograms(const bool loadOrNot); // Setter for loading jet pT closure histograms

  // Setter for jet-event plane correlation histograms
  void SetLoadJetEventPlaneHistograms(const bool loadOrNot);
  
  // Setters for ranges for different bins
  void SetCentralityBinRange(const int first, const int last);          // Setter for centrality bin range
  void SetJetPtBinRange(const int first, const int last);            // Setter for jet pT bin range in energy-energy correlator histograms

  // Getters for number of bins in histograms
  int GetNCentralityBins() const;          // Getter for the number of centrality bins
  int GetNJetPtBins() const;            // Getter for the number of jet pT bins in energy-energy correlator histograms
  double GetCentralityBinBorder(const int iCentrality) const;  // Getter for i:th centrality bin border
  double GetJetPtBinBorder(const int iJetPt) const;         // Getter for i:th jet pT bin border in energy-energy correlator histograms
  
  // Getters for histogram and axis naming
  const char* GetJetHistogramName(const int iJetType) const; // Getter for the jet histogram name
  const char* GetJetAxisName(const int iJetType) const;      // Getter for name suitable for x-axis in a jet histogram
  
  TString GetSystem() const;  // Getter for collision system
  
  // Getters for event information histograms
  TH1D* GetHistogramVertexZ() const;            // Getter for z-vertex histogram
  TH1D* GetHistogramVertexZWeighted() const;    // Getter for weighted z-vertex histogram
  TH1D* GetHistogramEvents() const;             // Getter for histogram for number of events surviving different event cuts
  TH1D* GetHistogramCentrality() const;         // Getter for centrality histogram in all events
  TH1D* GetHistogramCentralityWeighted() const; // Getter for weighted centrality histogram in all events
  
  // Getters for jet histograms
  TH1D* GetHistogramJetPt(int iCentrality, int iJetType, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;     // Jet pT histograms
  TH1D* GetHistogramInclusiveJetPt(int iCentrality, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;     // Inclusive jet pT histograms
  TH1D* GetHistogramLeadingJetPt(int iCentrality, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;     // Leading jet pT histograms
  TH1D* GetHistogramJetPhi(int iCentrality, int iJetType, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;    // Jet phi histograms
  TH1D* GetHistogramInclusiveJetPhi(int iCentrality, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;   // Inclusive jet phi histograms
  TH1D* GetHistogramLeadingJetPhi(int iCentrality, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;    // Leading jet phi histograms
  TH1D* GetHistogramJetEta(int iCentrality, int iJetType, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;    // Jet eta histograms
  TH1D* GetHistogramInclusiveJetEta(int iCentrality, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;   // Inclusive jet eta histograms
  TH1D* GetHistogramLeadingJetEta(int iCentrality, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;    // Leading jet eta histograms
  TH2D* GetHistogramJetEtaPhi(int iCentrality, int iJetType, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const; // 2D eta-phi histogram for jets
  TH2D* GetHistogramInclusiveJetEtaPhi(int iCentrality, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;  // 2D eta-phi histogram for inclusive jets
  TH2D* GetHistogramLeadingJetEtaPhi(int iCentrality, int iParton = JetBackgroundHistograms::knInitialPartonTypes, int iMatch = JetBackgroundHistograms::knMatchingTypes) const;  // 2D eta-phi histogram for leading jets 

  // Getter for jet pT closure histograms
  TH1D* GetHistogramJetPtClosure(const int iGenPtBin, const int iEtaBin, const int iPhiBin, const int iCentrality, const int iClosureParticle) const; // Jet pT closure

  // Getter for jet-event plane histograms
  TH1D* GetHistogramJetEventPlane(int iOrder, int iJetType, int iCentrality, int iJetPt = -1);
  TH1D* GetHistogramInclusiveJetEventPlane(int iOrder, int iCentrality, int iJetPt = -1);
  TH1D* GetHistogramLeadingJetEventPlane(int iOrder, int iCentrality, int iJetPt = -1);
  
  // Getters for the loaded centrality and track pT bins
  int GetFirstCentralityBin() const;          // Get the first loaded centrality bin
  int GetLastCentralityBin() const;           // Get the last loaded centrality bin
  int GetFirstJetPtBin() const;            // Get the first loaded energy-energy correlator jet pT bin
  int GetLastJetPtBin() const;             // Get the last loaded energy-energy correlator jet pT bin
  
  // Getters for normalization information
  int GetNEvents() const;                      // Getter for the number of events passing the cuts
  
  // Getter for the card
  JetBackgroundCard* GetCard() const;  // Getter for the JCard
  
private:
  
  // Data members
  TFile* fInputFile;                  // File from which the histograms are read
  JetBackgroundCard* fCard;                     // Card inside the data file for binning, cut collision system etc. information
  
  // ==============================================
  // ======== Flags for histograms to load ========
  // ==============================================
  
  bool fLoadEventInformation;                   // Load the event information histograms
  bool fLoadJets;                               // Load the jet histograms
  bool fLoad2DHistograms;                       // Load also two-dimensional (eta,phi) histograms
  bool fLoadJetPtClosureHistograms;             // Load jet pT closure histograms
  bool fLoadJetEventPlaneCorrelationHistograms; // Load jet-event plane correlation histograms
  
  // ==============================================
  // ======== Ranges of histograms to load ========
  // ==============================================
  
  int fFirstLoadedCentralityBin;  // First centrality bin that is loaded
  int fLastLoadedCentralityBin;   // Last centrality bin that is loaded
  int fFirstLoadedJetPtBin;    // First loaded jet pT bin for the energy-energy correlator histograms
  int fLastLoadedJetPtBin;     // Last loaded jet pT bin for the energy-energy correlator histograms
  
  // =============================================
  // ============ Binning information ============
  // =============================================
  int fCentralityBinIndices[kMaxCentralityBins+1];           // Indices for centrality bins in centrality binned histograms
  double fCentralityBinBorders[kMaxCentralityBins+1];        // Centrality bin borders, from which bin indices are obtained
  int fJetPtIndices[kMaxJetPtBins+1];                  // Indices for jet pT bins in energy-energy correlator histograms
  double fJetPtBinBorders[kMaxJetPtBins+1];            // Jet pT bin borders in energy-energy correlator histograms
  int fnCentralityBins;                                      // Number of centrality bins in the JCard of the data file
  int fnJetPtBins;                                        // Number of jet pT bins for the energy-energy correlator histograms

  // =============================================
  // ===== Histograms for the dijet analysis =====
  // =============================================

  // Event information histograms
  TH1D* fhVertexZ;            // Vertex z position
  TH1D* fhVertexZWeighted;    // Weighted vertex z-position (only meaningfull for MC)
  TH1D* fhEvents;             // Number of events surviving different event cuts
  TH1D* fhCentrality;         // Centrality of all events
  TH1D* fhCentralityWeighted; // Weighted centrality distribution in all events (only meaningful for MC)
  TH1D* fhPtHat;              // pT hat for MC events (only meaningful for MC)
  TH1D* fhPtHatWeighted;      // Weighted pT hat distribution (only meaningful for MC)

  // Histograms for jets
  TH1D* fhJetPt[knJetTypes][kMaxCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1][JetBackgroundHistograms::knMatchingTypes+1];      // Jet pT histograms
  TH1D* fhJetPhi[knJetTypes][kMaxCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1][JetBackgroundHistograms::knMatchingTypes+1];     // Jet phi histograms
  TH1D* fhJetEta[knJetTypes][kMaxCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1][JetBackgroundHistograms::knMatchingTypes+1];     // Jet eta histograms
  TH2D* fhJetEtaPhi[knJetTypes][kMaxCentralityBins][JetBackgroundHistograms::knInitialPartonTypes+1][JetBackgroundHistograms::knMatchingTypes+1];  // 2D eta-phi histogram for jets

  // Histograms for jet pT closure
  TH1D* fhJetPtClosure[knGenJetPtBins+1][knJetEtaBins+1][knJetPhiBins+1][kMaxCentralityBins][JetBackgroundHistograms::knInitialPartonTypes];  // Jet pT closure

  // Histograms for jet-event plane correlation
  TH1D* fhJetEventPlane[knJetTypes][JetBackgroundHistograms::knEventPlanes][kMaxCentralityBins][kMaxJetPtBins];

  // Private methods
  void InitializeFromCard(); // Initialize several member variables from JetBackgroundCard
  
  // Binning related methods
  void SetBinIndices(const char* histogramName, const int nBins, int* binIndices, const double* binBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetBinBordersAndIndices(const char* histogramName, const int nBins, double *copyBinBorders, int* binIndices, const double* binBorders, const int iAxis, const bool setIndices); // Read the bin indices for given bin borders
  
  // Finders for histograms with different amount of restrictions
  TH2D* FindHistogram2D(THnSparseD* histogramArray, int xAxis, int yAxis, int nAxes, int* axisNumber, int* lowBinIndex, int* highBinIndex, const bool normalizeToBinWidth = true); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH2D* FindHistogram2D(THnSparseD* histogramArray, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0, const bool normalizeToBinWidth = true); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(THnSparseD* histogramArray, int xAxis, int nAxes, int* axisNumber, int* lowBinIndex, int* highBinIndex, const bool normalizeToBinWidth = true); // Extract a histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(THnSparseD* histogramArray, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0, const bool normalizeToBinWidth = true); // Extract a histogram using given axis restrictions from THnSparseD
  
  // Loaders for different groups of histograms
  void LoadJetHistograms(); // Loader for jet histograms
  void LoadJetPtClosureHistograms(); // Loader for jet pT closure histograms
  void LoadJetEventPlaneHistograms(); // Loader for jet-event plane correlation histograms
  
  // Generic setter for bin indice and borders
  void SetGenericBins(const bool readBinsFromFile, const char* histogramName, const int iAxis, int nSetBins, double* setBinBorders, int* setBinIndices, const int nBins, const double* binBorders, const char* errorMessage, const int maxBins, const bool setIndices); // Generic bin setter
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int& first, int& last); // Sanity check for given binning
  int BinIndexCheck(const int nBins, const int binIndex) const; // Check that given index is in defined range
  
  // Methods for histogram writing
  void WriteJetHistograms();            // Write the jet histograms to the file that is currently open
  void WriteClosureHistograms();        // Write the closure histograms to the file that is currently open
  void WriteJetEventPlaneHistograms();  // Write the jet-event plane correlation histograms to the file that is currently open
  
};

#endif
