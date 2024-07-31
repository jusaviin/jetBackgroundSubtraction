/*
 * Implementation of the JetBackgroundCard class
 */

// Own includes
#include "JetBackgroundCard.h"

/*
 * Contructor with input file
 *
 *  TFile *inFile = Input file
 */
JetBackgroundCard::JetBackgroundCard(TFile *inFile):
  fInputFile(inFile),
  fCardDirectory("JCard"),
  fGitHash(0),
  fProjectionGitHash(0)
{

  // Read the vectors from the input file
  fInputFile->cd(fCardDirectory.Data());
  ReadVectors();
}

/*
 * Reader for all the vectors from the input file
 */
void JetBackgroundCard::ReadVectors(){
  
  // Read the git hash
  fGitHash = (TObjString*) gDirectory->Get("GitHash");
  fProjectionGitHash = (TObjString*) gDirectory->Get("ProjectionGitHash");

  // Read the TVectorT<float>:s
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    fCardEntries[iEntry] = (TVectorT<float>*) gDirectory->Get(fCardEntryNames[iEntry]);
  }
  
  // Read the input file name
  fInputFileName = (TObjString*) gDirectory->Get(fInputFileSaveName);

}

/*
 *  Getter for jet type
 */
int JetBackgroundCard::GetJetType() const{
  return (*fCardEntries[kJetType])[1];
}

/*
 *  Getter for jet background subtration algorithm
 */
int JetBackgroundCard::GetJetSubtraction() const{
  return (*fCardEntries[kJetSubtraction])[1];
}

// Getter for the minimum jet pT cut
double JetBackgroundCard::GetJetPtCut() const{
  return (*fCardEntries[kMinPtCut])[1];
}

/*
 * Get the number of bins for internal index
 * If no vector is found in the index, return 0.
 */
int JetBackgroundCard::GetNBins(const int index) const{
  if(fCardEntries[index]) return fCardEntries[index]->GetNoElements()-1;
  return 0;
}

// Get the number of centrality bins
int JetBackgroundCard::GetNCentralityBins() const{
  return GetNBins(kCentralityBinEdges);
}

// Get the number of track pT bins
int JetBackgroundCard::GetNTrackPtBins() const{
  return GetNBins(kTrackPtBinEdges);
}

// Get the number of jet pT bins
int JetBackgroundCard::GetNJetPtBins() const{
  return GetNBins(kJetPtBinEdges);
}

/*
 * Get a bin index based on a given value.
 * If value is out of bounds, return -1
 */
int JetBackgroundCard::GetBinIndex(const int index, const double value) const{
  
  // Find the number of bins in the array
  int nBins = GetNBins(index);
  
  // If the number given is smaller than the first value, it is out of bounds. Return -1 in this case.
  if(value < (*fCardEntries[index])[1]) return -1;
  
  // Find the bin in which the given value is and return it
  for(int iBin = 2; iBin <= nBins+1; iBin++){
    if(value < (*fCardEntries[index])[iBin]) return iBin-2;
  }
  
  // If a value was not found, it must be larger than anything in the array. Return -1 in this case.
  return -1;
}

// Get the bin index for a given centrality value
int JetBackgroundCard::GetBinIndexCentrality(const double value) const{
  return GetBinIndex(kCentralityBinEdges,value);
}

// Get the bin index for a given track pT value
int JetBackgroundCard::GetBinIndexTrackPt(const double value) const{
  return GetBinIndex(kTrackPtBinEdges,value);
}

// Get the bin index for a given jet pT
int JetBackgroundCard::GetBinIndexJetPt(const double value) const{
  return GetBinIndex(kJetPtBinEdges,value);
}

/*
 * Find if a bin with given borders exist in an array defined by internal index
 *
 *  const int index = Internal index defining the array to search the bin from
 *  const double lowBorder = Low bin border of the searched bin
 *  const double highBorder = High bin border of the searched bin
 *
 *  return: Return the index of the corresponding bin. If the bin is not found, return -1
 */
int JetBackgroundCard::FindBinIndex(const int index, const double lowBorder, const double highBorder) const{
  
  // Find the number of bins in the array
  int nBins = GetNBins(index);

  // Define a small number
  double epsilon = 0.001;
 
  // Go through all the bins and check if there exist a bin that matches with in given bins
  for(int iBin = 0; iBin < nBins; iBin++){
    if(TMath::Abs(GetLowBinBorder(index, iBin) - lowBorder) < epsilon){
      if(TMath::Abs(GetHighBinBorder(index, iBin) - highBorder) < epsilon){
        return iBin;
      }
    }
  }
  
  // If a match for bin borders was not found, the given bin is not in the array. Return -1 to show that.
  return -1;
}

// Find if a centrality bin with given borders exists and return its index
int JetBackgroundCard::FindBinIndexCentrality(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kCentralityBinEdges,lowBorder,highBorder);
}

// Find if a centrality bin with given borders exists and return its index
int JetBackgroundCard::FindBinIndexCentrality(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kCentralityBinEdges,binBorders.first,binBorders.second);
}

// Find if a track pT bin with given borders exists and return its index
int JetBackgroundCard::FindBinIndexTrackPt(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kTrackPtBinEdges,lowBorder,highBorder);
}

// Find if a track pT bin with given borders exists and return its index
int JetBackgroundCard::FindBinIndexTrackPt(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kTrackPtBinEdges,binBorders.first,binBorders.second);
}

// Find if a jet pT bin with given borders exists and return its index
int JetBackgroundCard::FindBinIndexJetPt(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kJetPtBinEdges,lowBorder,highBorder);
}

// Find if a jet pT bin with given borders exists and return its index
int JetBackgroundCard::FindBinIndexJetPt(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kJetPtBinEdges,binBorders.first,binBorders.second);
}

// Get the low border of i:th bin from internal index
double JetBackgroundCard::GetLowBinBorder(const int index, const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin > GetNBins(index)) return -1;
  
  // Return the asked bin index
  if(fCardEntries[index]) return (*fCardEntries[index])[iBin+1];
  return -1;
}

// Get the low border of i:th centrality bin
double JetBackgroundCard::GetLowBinBorderCentrality(const int iBin) const{
  return GetLowBinBorder(kCentralityBinEdges,iBin);
}

// Get the low border of i:th track pT bin
double JetBackgroundCard::GetLowBinBorderTrackPt(const int iBin) const{
  return GetLowBinBorder(kTrackPtBinEdges,iBin);
}

// Get the low border of i:th jet pT bin
double JetBackgroundCard::GetLowBinBorderJetPt(const int iBin) const{
  return GetLowBinBorder(kJetPtBinEdges,iBin);
}

// Get the high border of i:th bin from internal index
double JetBackgroundCard::GetHighBinBorder(const int index, const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin >= GetNBins(index)) return -1;
  
  // Return the asked bin index
  if(fCardEntries[index]) return (*fCardEntries[index])[iBin+2];
  return -1;
}

// Get the high border of i:th centrality bin
double JetBackgroundCard::GetHighBinBorderCentrality(const int iBin) const{
  return GetHighBinBorder(kCentralityBinEdges,iBin);
}

// Get the high border of i:th track pT bin
double JetBackgroundCard::GetHighBinBorderTrackPt(const int iBin) const{
  return GetHighBinBorder(kTrackPtBinEdges,iBin);
}

// Get the high border of i:th jet pT bin
double JetBackgroundCard::GetHighBinBorderJetPt(const int iBin) const{
  return GetHighBinBorder(kJetPtBinEdges,iBin);
}

// Get the bin borders of the i:th centrality bin
std::pair<double,double> JetBackgroundCard::GetBinBordersCentrality(const int iBin) const{
  return std::make_pair(GetLowBinBorderCentrality(iBin), GetHighBinBorderCentrality(iBin)); 
}

// Get the bin borders of the i:th track pT bin
std::pair<double,double> JetBackgroundCard::GetBinBordersTrackPt(const int iBin) const{
  return std::make_pair(GetLowBinBorderTrackPt(iBin), GetHighBinBorderTrackPt(iBin)); 
}

// Get the bin borders of the i:th jet pT bin
std::pair<double,double> JetBackgroundCard::GetBinBordersJetPt(const int iBin) const{
  return std::make_pair(GetLowBinBorderJetPt(iBin), GetHighBinBorderJetPt(iBin)); 
}

/*
 * Add one-dimensional vector to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added vector in entry array
 *  float entryContent = Content to be put into the vector
 */
void JetBackgroundCard::AddOneDimensionalVector(int entryIndex, float entryContent){
  
  // Make a new one dimensional vector to the desired index with given content
  float contents[1] = {entryContent};
  fCardEntries[entryIndex] = new TVectorT<float>(1,1,contents);

}

/*
 * Add one-dimensional vector to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added vector in entry array
 *  int dimension = Number of entries in the given array
 *  float *contents = Content to be put into the vector
 */
void JetBackgroundCard::AddVector(int entryIndex, int dimension, double *contents){
  
  // Convert double pointer to float pointer
  float* convertedContents = new float[dimension];
  for(int i = 0; i < dimension; i++){
    convertedContents[i] = contents[i];
  }
  
  // Make a new vector to the desired index with given content
  fCardEntries[entryIndex] = new TVectorT<float>(1,dimension,convertedContents);
  
  // Delete the converted contents array
  delete[] convertedContents;
  
}

/*
 * Add file name to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added file name in file name array
 *  TString fileName = Added file name
 */
void JetBackgroundCard::AddFileName(TString fileName){
  
  // Make convert the string to TObjString and add it to the file name object
  fInputFileName = new TObjString(fileName.Data());
  
}

/*
 * Add git hash for the projecting code used to project the histograms to the card
 *
 * Arguments:
 *  const char* gitHash = Git hash to be added for projection
 */
void JetBackgroundCard::AddProjectionGitHash(const char* gitHash){
  
  // Convert the const char* to TObjString and assign it to projection git hash
  fProjectionGitHash = new TObjString(gitHash);
  
}

/*
 * Print the contents of the card to the console
 */
void JetBackgroundCard::Print() const{

  const char* gitHash = "Unknown";
  if(fGitHash != NULL) gitHash = fGitHash->String().Data();
  
  std::cout<<std::endl<<"========================= JetBackgroundCard =========================="<<std::endl;
  std::cout << Form("%25s","GitHash"); //print keyword
  std::cout << ": " << gitHash << std::endl;
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    if(fCardEntries[iEntry]){
      std::cout << Form("%25s",fCardEntryNames[iEntry]); //print keyword
      std::cout << " (dim = "<<fCardEntries[iEntry]->GetNoElements() << ") ";//print size of TVector
      for(int iElement = 1; iElement <= fCardEntries[iEntry]->GetNoElements(); iElement++){
        std::cout << (*fCardEntries[iEntry])[iElement] << " ";//TVector components
      }
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;
  
  if(fProjectionGitHash != NULL) std::cout << "Git hash for projections: " << fProjectionGitHash->String().Data() << std::endl;
  
  if(fInputFileName) std::cout << "Used input file: " << fInputFileName->String().Data() << std::endl;
    
}

/*
 * Write the contents of the JetBackgroundCard to a file
 */
void JetBackgroundCard::Write(TDirectory* file){
  
  // Create a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write the git hash to the file
  if(fGitHash) fGitHash->Write("GitHash");
  
  // Write all the vectors to the file. Not all of these exist in older versions of cards, thus check if exists before writing.
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    if(fCardEntries[iEntry])  fCardEntries[iEntry]->Write(fCardEntryNames[iEntry]);
  }
   
  // Write the git hash used for projections to the file
  if(fProjectionGitHash) fProjectionGitHash->Write("ProjectionGitHash");
  
  // Write all the data names to the file.
  if(fInputFileName) fInputFileName->Write(fInputFileSaveName);
  
  // Return back to the main directory
  file->cd("../");
}
