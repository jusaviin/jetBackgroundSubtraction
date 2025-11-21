#ifndef JETPHIFLATTENER_H
#define JETPHIFLATTENER_H

// Root includes
#include <TH1.h>
#include <TFile.h>

/*
 * JetPhiFlattener class
 *
 * For testing purposes, I want to see if event plane correlation becomes more even if jet phi distribution is flat
 */
class JetPhiFlattener {
  
public:
 
  // Maximum dimensions for the arrays
  static const int kMaxCentralityBins = 5;      // Maximum allowed number of centrality bins
  
  JetPhiFlattener();                        // Contructor
  JetPhiFlattener(TString inputFileName);   // Custom constructor
  ~JetPhiFlattener() = default;        // Destructor
  
  double GetJetPhiFlattener(const double centrality, const double phi) const;

  void ReadCorrectionTables(); // Initialize the correction tables from file
  
private:
  
  // File from which the histograms are read
  TFile* fInputFile;
  
  // Binning information to be read from the input file
  double fCentralityBinBorders[5] = {4,14,34,54,94};
  int fnCentralityBins = 4;
  
  // Histograms from which the corrections are read
  TH1D* fCorrectionTable[kMaxCentralityBins];
  
  // Methods
  int FindCentralityBin(const double centrality) const; // Find centrality bin index
  
};

#endif
