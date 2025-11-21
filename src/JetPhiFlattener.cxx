/*
 * Implementation of the JetPhiFlattener class
 */

// Own includes
#include "JetPhiFlattener.h"

/*
 * Contructor
 */
JetPhiFlattener::JetPhiFlattener() :
  fInputFile(nullptr)
{
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    fCorrectionTable[iCentrality] = nullptr;
  } // Centrality loop
}

/*
 * Custom constructor
 */
JetPhiFlattener::JetPhiFlattener(TString inputFileName)
{
  fInputFile = TFile::Open(inputFileName);
  ReadCorrectionTables();
}

/*
 * Read the track pair efficiency correction tables from the file
 */
void JetPhiFlattener::ReadCorrectionTables(){

  // Read all the tables from the file
  for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
    fCorrectionTable[iCentrality] = (TH1D*) fInputFile->Get(Form("phiCorrection_C%d", iCentrality));
  }
  
}

/*
 * Find centrality bin index
 */
int JetPhiFlattener::FindCentralityBin(const double centrality) const{

  // If the value is below the lowest index, just use the lowest index
  if(centrality < fCentralityBinBorders[0]) return 0;
  
  // Find the value from the arrey and return the bin index
  for(int iBin = 0; iBin < fnCentralityBins; iBin++){
    if(centrality < fCentralityBinBorders[iBin+1]) return iBin;
  }
  
  // If the value is higher than the upper bin border, just use the highest index
  return fnCentralityBins-1;
}


/*
 * Get the weight to flatten the jet phi distribution
 *
 *  const double centrality = Centrality of the event
 *  const double phi = Jet phi value
 *
 *  return: Flattening factor for this centrality and phi value
 */
double JetPhiFlattener::GetJetPhiFlattener(const double centrality, const double phi) const{
  

  // Determine the correct bin from the centrality and track pT values
  int iCentrality = FindCentralityBin(centrality);
  return fCorrectionTable[iCentrality]->GetBinContent(fCorrectionTable[iCentrality]->FindBin(phi));

}