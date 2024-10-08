/*
 * Library of useful algorithms, easily accessible from a class
 */

// Own includes
#include "AlgorithmLibrary.h"

/*
 *  Constructor
 */
AlgorithmLibrary::AlgorithmLibrary()
{
  // Constructor
}

/*
 *  Copy constructor
 */
AlgorithmLibrary::AlgorithmLibrary(const AlgorithmLibrary& in)
{
  // Copy constructor
}

/*
 * Assingment operator
 */
AlgorithmLibrary& AlgorithmLibrary::operator=(const AlgorithmLibrary& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  return *this;
}

/*
 *  Destructor
 */
AlgorithmLibrary::~AlgorithmLibrary()
{
  // Destructor
}

/*
 * Find the minimum and maximum values from a histogram. They must be more extreme than the current values
 *
 *  Arguments: TH1D* histogram = Histogram from which the minimum and maximum values are searched
 *             std::pair<double,double> currentMinMax = The found values need to be more extreme than these to be accepted
 */
std::pair<double,double> AlgorithmLibrary::FindHistogramMinMax(TH1D* histogram, std::pair<double,double> currentMinMax){
  
  // Use the whole histogram as a search range
  std::pair<int,int> searchRange = std::make_pair(1, histogram->GetNbinsX());
  return FindHistogramMinMax(histogram, currentMinMax, searchRange);
  
}

/*
 * Find the minimum and maximum values from a histogram. They must be more extreme than the current values
 *
 *  Arguments: TH1D* histogram = Histogram from which the minimum and maximum values are searched
 *             std::pair<double,double> currentMinMax = The found values need to be more extreme than these to be accepted
 *             std::pair<double,double> searchRange = bin value range from which the minimum and maximum values are searched
 */
std::pair<double,double> AlgorithmLibrary::FindHistogramMinMax(TH1D* histogram, std::pair<double,double> currentMinMax, std::pair<double,double> searchRange){
  
  // Find the bin range from which the minimum and maximum values are searched
  double epsilon = 0.00001;
  int firstBin = histogram->GetXaxis()->FindBin(searchRange.first+epsilon);
  int lastBin = histogram->GetXaxis()->FindBin(searchRange.second-epsilon);

  std::pair<int,int> searchBins = std::make_pair(firstBin, lastBin);

  // Return the minimum and maximum values from the histogram
  return FindHistogramMinMax(histogram, currentMinMax, searchBins);
  
}

/*
 * Find the minimum and maximum values from a histogram. They must be more extreme than the current values
 *
 *  Arguments: TH1D* histogram = Histogram from which the minimum and maximum values are searched
 *             std::pair<double,double> currentMinMax = The found values need to be more extreme than these to be accepted
 *             std::pair<int,int> searchRange = Bin index range from which the minimum and maximum values are searched
 */
std::pair<double,double> AlgorithmLibrary::FindHistogramMinMax(TH1D* histogram, std::pair<double,double> currentMinMax, std::pair<int,int> searchRange){
  
  // As initial guess, take the given minimum and maximum values
  std::pair<double,double> newMinMax = std::make_pair(currentMinMax.first, currentMinMax.second);
  
  // Loop through all the bins in the histogram and update the minimum and maximum values
  double currentValue, currentError;
  for(int iBin = searchRange.first; iBin <= searchRange.second; iBin++){
    currentValue = histogram->GetBinContent(iBin);
    currentError = histogram->GetBinError(iBin);
    if((currentValue-currentError) < newMinMax.first) newMinMax.first = currentValue-currentError;
    if((currentValue+currentError) > newMinMax.second) newMinMax.second = currentValue+currentError;
  }
  
  // Return the minimum and maximum values from the histogram
  return newMinMax;
  
}

/*
 * Rebin one dimensional histogram with asymmetric bin edges
 *
 * Arguments:
 *  TH1D* histogramInNeedOfRebinning = Histogram to be rebinned
 *  const int nBins = Number of bins for the rebinned histogram
 *  const double* binEdges = Bin edges for the rebinned histogram
 */
TH1D* AlgorithmLibrary::RebinAsymmetric(TH1D* histogramInNeedOfRebinning, const int nBins, const double* binEdges){
  
  // First, check that the new bin boundaries are also bin boundaries in the original histogram
  bool binsGood = CheckBinBoundaries(nBins,binEdges,histogramInNeedOfRebinning->GetXaxis());
  if(!binsGood){
    std::cout << "Cannot rebin histogram " << histogramInNeedOfRebinning->GetName() << " because given bin borders do not match with the bin borders of the original histogram!" << std::endl;
    return histogramInNeedOfRebinning;
  }
  
  // Clone the original histogram
  TH1D *clonedHistogram = (TH1D*) histogramInNeedOfRebinning->Clone("_rebinned");
  
  // Set the new binning for histogram (destroys content in each bin)
  clonedHistogram->SetBins(nBins,binEdges);
  
  // Make sure that each bin is set to zero
  for(int iBin = 1; iBin <= clonedHistogram->GetNbinsX(); iBin++){
    clonedHistogram->SetBinContent(iBin,0);
    clonedHistogram->SetBinError(iBin,0);
  }
  
  // Add the contents back to the histogram that was rebinned
  double binContent, binError, binCenter, oldContent, oldError;
  int newBin;
  double binWidth;
  for(int iBin = 1; iBin <= histogramInNeedOfRebinning->GetNbinsX(); iBin++){
    
    // Read the contents from the non-rebinned histogram
    binWidth = histogramInNeedOfRebinning->GetBinWidth(iBin);
    binContent = histogramInNeedOfRebinning->GetBinContent(iBin)*binWidth;  // Remove previous bin width normalization
    binError = histogramInNeedOfRebinning->GetBinError(iBin)*binWidth;      // Remove previous bin width normalization
    binCenter = histogramInNeedOfRebinning->GetBinCenter(iBin);
    
    // Add the contents to the rebinned histgram
    newBin = clonedHistogram->FindBin(binCenter);
    oldContent = clonedHistogram->GetBinContent(newBin);
    oldError = clonedHistogram->GetBinError(newBin);
    clonedHistogram->SetBinContent(newBin,binContent+oldContent);
    clonedHistogram->SetBinError(newBin,TMath::Sqrt(binError*binError+oldError*oldError));
  }
  
  // Normalize the bin contents to bin width
  for(int iBin = 1; iBin <= clonedHistogram->GetNbinsX(); iBin++){
    binWidth = clonedHistogram->GetBinWidth(iBin);
    binContent = clonedHistogram->GetBinContent(iBin);
    binError = clonedHistogram->GetBinError(iBin);
    clonedHistogram->SetBinContent(iBin,binContent/binWidth);
    clonedHistogram->SetBinError(iBin,binError/binWidth);
  }
  
  // Return the rebinned histogram
  return clonedHistogram;
}

/*
 * Rebin a two-dimensional histogram
 *
 *  Arguments:
 *   TH2D* histogramInNeedOfRebinning = The two dimensional histogram that is going to be rebinned
 *   const int nBinsX = Number of bins after rebinning in the x-axis
 *   const double* binBordersX = Bin borders for rebinning for the x-axis
 *   const int nBinsY = Number of bins after rebinning in the y-axis
 *   const double* binBordersY = Bin borders for rebinning for the y-axis
 *   const bool undoBinArea = Undo bin area normalization from the input histogram
 *   const bool normalizeBinArea = Normalize the output histogram contents with bin area
 *
 *   return: Rebinned histogram
 */
TH2D* AlgorithmLibrary::RebinHistogram(TH2D* histogramInNeedOfRebinning, const int nBinsX, const double* binBordersX, const int nBinsY, const double* binBordersY, const bool undoBinArea, const bool normalizeBinArea){

  // First, check that the new bin boundaries are also bin boundaries in the original histogram
  bool binsGood = CheckBinBoundaries(nBinsX, binBordersX, histogramInNeedOfRebinning->GetXaxis());
  if(!binsGood){
    std::cout << "Cannot rebin histogram " << histogramInNeedOfRebinning->GetName() << " because of a bin edge problem in x-axis!" << std::endl;
    return histogramInNeedOfRebinning;
  }
  
  binsGood = CheckBinBoundaries(nBinsY, binBordersY, histogramInNeedOfRebinning->GetYaxis());
  if(!binsGood){
    std::cout << "Cannot rebin histogram " << histogramInNeedOfRebinning->GetName() << " because of a bin edge problem in y-axis!" << std::endl;
    return histogramInNeedOfRebinning;
  }
  
  // Root does not offer a method to directly rebin a two-dimensional histogram, so I have implemented my own
  // Helper variables for rebinning
  double currentBinContent;
  double currentBinError;
  double nonRebinnedContent;
  double nonRebinnedError;
  int newHistogramIndex;
  double xValue;
  double yValue;
  double binWidthX;
  double binWidthY;
  double binArea = 1;
  
  // Create the histogram with new binning
  TString newName = Form("%sRebinned",histogramInNeedOfRebinning->GetName());
  TH2D *rebinnedHistogram = new TH2D(newName,newName,nBinsX,binBordersX,nBinsY,binBordersY);
  
  // Loop over all the bins in the old histogram and insert the content to the new histogram
  for(int iBinX = 1; iBinX <= histogramInNeedOfRebinning->GetNbinsX(); iBinX++){
    xValue = histogramInNeedOfRebinning->GetXaxis()->GetBinCenter(iBinX);
    binWidthX = histogramInNeedOfRebinning->GetXaxis()->GetBinWidth(iBinX);
    
    for(int iBinY = 1; iBinY <= histogramInNeedOfRebinning->GetNbinsY(); iBinY++){
      yValue = histogramInNeedOfRebinning->GetYaxis()->GetBinCenter(iBinY);
      binWidthY = histogramInNeedOfRebinning->GetYaxis()->GetBinWidth(iBinY);
      
      // If requested, undo bin area normalization from the input histogram
      if(undoBinArea) binArea = binWidthX*binWidthY;
      
      // Find the global bin index from the new histogram correcponding to the bin in the old histogram
      newHistogramIndex = rebinnedHistogram->FindBin(xValue,yValue);
      
      // Add the bin content from the old histogram to the new histogram, adding errors in quadrature
      nonRebinnedContent = histogramInNeedOfRebinning->GetBinContent(iBinX,iBinY) * binArea;
      nonRebinnedError = histogramInNeedOfRebinning->GetBinError(iBinX,iBinY) * binArea;
      currentBinContent = rebinnedHistogram->GetBinContent(newHistogramIndex);
      currentBinError = rebinnedHistogram->GetBinError(newHistogramIndex);
      rebinnedHistogram->SetBinContent(newHistogramIndex,currentBinContent+nonRebinnedContent);
      rebinnedHistogram->SetBinError(newHistogramIndex,TMath::Sqrt(TMath::Power(currentBinError,2)+TMath::Power(nonRebinnedError,2)));
    } // DeltaEta loop
  } // DeltaPhi loop
  
  // After rebinning, do bin area normalization for each bin if requested
  if(normalizeBinArea){
    for(int iBinX = 1; iBinX <= rebinnedHistogram->GetNbinsX(); iBinX++){
      binWidthX = rebinnedHistogram->GetXaxis()->GetBinWidth(iBinX);
      for(int iBinY = 1; iBinY <= rebinnedHistogram->GetNbinsY(); iBinY++){
        binWidthY = rebinnedHistogram->GetYaxis()->GetBinWidth(iBinY);
        currentBinContent = rebinnedHistogram->GetBinContent(iBinX,iBinY);
        currentBinError = rebinnedHistogram->GetBinError(iBinX,iBinY);
        rebinnedHistogram->SetBinContent(iBinX,iBinY,currentBinContent/(binWidthY*binWidthX));
        rebinnedHistogram->SetBinError(iBinX,iBinY,currentBinError/(binWidthY*binWidthX));
      } // y-axis loop
    } // x-axis loop
  } // Normalizing by bin area
  
  return rebinnedHistogram;
  
}

/*
 * Checker that new bin boundaries are compatible with the old ones
 *
 *  Arguments:
 *   const int nCheckedBins = Number of new bins to be checked
 *   const double* checkedBins = Bin boundaries to be checked
 *   const TAxis* originalAxis = Original axis against with the new bins are checked
 *
 *   return: True, if all the new bin boundaries can be found from the original axis. False if not.
 */
bool AlgorithmLibrary::CheckBinBoundaries(const int nCheckedBins, const double* checkedBins, const TAxis* originalAxis){
  
  // Flag, if the current bin is a bin boundary in the histogram to be rebinned
  bool binOK = false;
  
  // First, check that the bin boundaries for the rebinned histogram match with the old bin boundaries
  for(int iCheckedBin = 0; iCheckedBin < nCheckedBins + 1; iCheckedBin++){
    binOK = false;
    for(int iOldBin = 1; iOldBin <= originalAxis->GetNbins()+1; iOldBin++){
      
      // We the bin edge is close enough to one original bin, accept the bin
      if(TMath::Abs(originalAxis->GetBinLowEdge(iOldBin)-checkedBins[iCheckedBin]) < 1e-4){
        binOK = true;
        break;
      }
      
    } // Loop over bins in the original axis
    if(!binOK){ // If the bin is not in original histogram, print error message and return false
      std::cout << "The bin boundary " << checkedBins[iCheckedBin] << " is not a bin boundary in the original histogram!" << std::endl;
      return false;
    }
  } // Loop over bins to be checked
  
  // If all is good, return true
  return true;
}

// Normalize all the columns of a 2-D histogram to a given value
void AlgorithmLibrary::NormalizeMatrix(TH2D* histogramInNeedOfNormalization, const double value, const int direction){
  
  if(direction == 1){
    NormalizeColumns(histogramInNeedOfNormalization, value);
  } else {
    NormalizeRows(histogramInNeedOfNormalization, value);
  }
  
}

// Normalize all the columns of a 2-D histogram to a given value
void AlgorithmLibrary::NormalizeColumns(TH2D* histogramInNeedOfNormalization, const double value){
  
  // Helper variables
  double binContent;
  double binError;
  double binSum;
  double normalizationFactor;
  
  // Loop over all the y-bins for a given x-bin and normalize the content
  for(int iBinX = 1; iBinX <= histogramInNeedOfNormalization->GetNbinsX(); iBinX++){
    
    binSum = 0;
    
    // First calculate the sum of contents in the y-bins
    for(int iBinY = 1; iBinY <= histogramInNeedOfNormalization->GetNbinsY(); iBinY++){
      binSum += histogramInNeedOfNormalization->GetBinContent(iBinX, iBinY);
    }
    
    // If there is no content, no need to normalize
    if(binSum == 0) continue;
    
    // Get the normalization factor
    normalizationFactor = value / binSum;
    
    // Normalize the bin content in each bin
    for(int iBinY = 1; iBinY <= histogramInNeedOfNormalization->GetNbinsY(); iBinY++){
      binContent = histogramInNeedOfNormalization->GetBinContent(iBinX, iBinY);
      binError = histogramInNeedOfNormalization->GetBinError(iBinX, iBinY);
      binContent *= normalizationFactor;
      binError *= normalizationFactor;
      histogramInNeedOfNormalization->SetBinContent(iBinX, iBinY, binContent);
      histogramInNeedOfNormalization->SetBinError(iBinX, iBinY, binError);
    }
  }
  
}

// Normalize all the rows of a 2-D histogram to a given value
void AlgorithmLibrary::NormalizeRows(TH2D* histogramInNeedOfNormalization, const double value){
  
  // Helper variables
  double binContent;
  double binError;
  double binSum;
  double normalizationFactor;
  
  // Loop over all the x-bins for a given y-bin and normalize the content
  for(int iBinY = 1; iBinY <= histogramInNeedOfNormalization->GetNbinsY(); iBinY++){
    
    binSum = 0;
    
    // First calculate the sum of contents in the y-bins
    for(int iBinX = 1; iBinX <= histogramInNeedOfNormalization->GetNbinsX(); iBinX++){
      binSum += histogramInNeedOfNormalization->GetBinContent(iBinX, iBinY);
    }
    
    // If there is no content, no need to normalize
    if(binSum == 0) continue;
    
    // Get the normalization factor
    normalizationFactor = value / binSum;
    
    // Normalize the bin content in each bin
    for(int iBinX = 1; iBinX <= histogramInNeedOfNormalization->GetNbinsX(); iBinX++){
      binContent = histogramInNeedOfNormalization->GetBinContent(iBinX, iBinY);
      binError = histogramInNeedOfNormalization->GetBinError(iBinX, iBinY);
      binContent *= normalizationFactor;
      binError *= normalizationFactor;
      histogramInNeedOfNormalization->SetBinContent(iBinX, iBinY, binContent);
      histogramInNeedOfNormalization->SetBinError(iBinX, iBinY, binError);
    }
  }
  
}

/*
 * Method for rotating a two dimensional histogram. Assumes constant bin size.
 *
 *  Arguments:
 *   TH2D* originalHistogram = Histogram to be rotated
 *
 *   return: Histogram where x and y axes have been switched
 */
TH2D* AlgorithmLibrary::RotateHistogram(TH2D* originalHistogram){
  
  // Find the binning information
  double nBinsX = originalHistogram->GetNbinsX();
  double nBinsY = originalHistogram->GetNbinsY();
  double lowBinX = originalHistogram->GetXaxis()->GetBinLowEdge(1);
  double highBinX = originalHistogram->GetXaxis()->GetBinUpEdge(nBinsX);
  double lowBinY = originalHistogram->GetYaxis()->GetBinLowEdge(1);
  double highBinY = originalHistogram->GetYaxis()->GetBinUpEdge(nBinsY);
  
  // Create a new histogram with inverted axes
  TString newName = Form("%sRotated", originalHistogram->GetName());
  TString newTitle = Form("%sRotated", originalHistogram->GetTitle());
  TH2D *rotatedHistogram = new TH2D(newName, newTitle, nBinsY, lowBinY, highBinY, nBinsX, lowBinX, highBinX);
  
  // Fill the new histogram
  double binValue, binError;
  for(int iX = 1; iX <= nBinsX; iX++){
    for(int iY = 1; iY <= nBinsY; iY++){
      binValue = originalHistogram->GetBinContent(iX, iY);
      binError = originalHistogram->GetBinError(iX, iY);
      rotatedHistogram->SetBinContent(iY, iX, binValue);
      rotatedHistogram->SetBinError(iY, iX, binError);
    }
  }
  
  // Return the new histogram
  return rotatedHistogram;
  
}

/*
 * Transform histogram to another one that shows the relative uncertainty of the original histogram in each bin
 *
 *  TH1D* transformedHistogram = Histogram in need of transformation
 *  const bool centerAtOne = True: Put the central value at one and relative systematic uncertainty around it.
 *                           False: Show relative systematic uncertainty as the histogram value. Disable errors.
 */
void AlgorithmLibrary::TransformToRelativeUncertainty(TH1D* transformedHistogram, const bool centerAtOne){
  
  double binContent;
  double binError;
  double relativeError;

  for(int iBin = 1; iBin <= transformedHistogram->GetNbinsX(); iBin++){
    binContent = transformedHistogram->GetBinContent(iBin);
    binError = transformedHistogram->GetBinError(iBin);
    relativeError = binError / binContent;

    if(centerAtOne){
      transformedHistogram->SetBinContent(iBin, 1);
      transformedHistogram->SetBinError(iBin, relativeError);
    } else {
      transformedHistogram->SetBinContent(iBin, relativeError);
      transformedHistogram->SetBinError(iBin, 0);
    }
  }

}

/*
 * Transform histogram describing relative uncertainties to one describing absolute uncertainties
 *
 *  TH1D* transformedHistogram = Histogram in need of transformation
 *  TH1D* absoluteScaleHistogram = Histogram giving the absolute scale for transformation
 *  const bool centerAtOne = True: The central value  in the relative uncertainty histogram is at one and relative systematic uncertainty around it.
 *                           False: The relative uncertainty is desribed in the histogram value in the relative uncertainty histogram
 */
void AlgorithmLibrary::TransformToAbsoluteUncertainty(TH1D* transformedHistogram, TH1D* absoluteScaleHistogram, const bool centerAtOne){

  double binContent;
  double binError;
  double relativeError;

  for(int iBin = 1; iBin <= transformedHistogram->GetNbinsX(); iBin++){

    if(centerAtOne){
      relativeError = transformedHistogram->GetBinError(iBin);
    } else {
      relativeError = transformedHistogram->GetBinContent(iBin);
    }

    binContent = absoluteScaleHistogram->GetBinContent(iBin);
    binError = binContent*relativeError;

    transformedHistogram->SetBinContent(iBin, binContent);
    transformedHistogram->SetBinError(iBin, binError);
  }

}

/*
 * Suppress single bin fluctuations in the fluctuating histogram
 *
 *  TH1D* fluctuatingHistogram = Histogram suffering from fluctuations
 *  const double lowRange = Lowest value in x-axis taken into account in fluctuation suppression
 *  const double highRange = Highest value on x-axis taken into account in fluctuation suppression
 *  const double threshold = Definition of a fluctuation that is too big for a single bin
 *  const double suppressionLevel = Definition how much fluctuating bins are suppressed
 */
void AlgorithmLibrary::SuppressSingleBinFluctuations(TH1D* fluctuatingHistogram, const double lowRange, const double highRange, const double threshold, const double suppressionLevel){

  double epsilon = 0.00001;
  int firstBin = fluctuatingHistogram->GetXaxis()->FindBin(lowRange+epsilon);
  int lastBin = fluctuatingHistogram->GetXaxis()->FindBin(highRange-epsilon);

  // Fluctuations are not suppressed in the first and last bin of the range, but they are in all the bins in between
  double previousValue, currentValue, nextValue, biggestNeighbor;
  for(int iBin = firstBin+1; iBin < lastBin; iBin++){
    previousValue = fluctuatingHistogram->GetBinContent(iBin-1);
    currentValue = fluctuatingHistogram->GetBinContent(iBin);
    nextValue = fluctuatingHistogram->GetBinContent(iBin+1);

    if(currentValue > threshold*previousValue && currentValue > threshold*nextValue){
      biggestNeighbor = previousValue;
      if(nextValue > previousValue) biggestNeighbor = nextValue;
      fluctuatingHistogram->SetBinContent(iBin, biggestNeighbor+biggestNeighbor*suppressionLevel);
    }
  }

}

/*
 * Get a const char* with today's date
 */
TString AlgorithmLibrary::GetToday(){

  // Get the current time_point
  std::chrono::system_clock::time_point now = std::chrono::system_clock::now();

  // Convert the time_point to time_t
  std::time_t now_t = std::chrono::system_clock::to_time_t(now);

  // Convert the time_t to tm in local timezone
  std::tm* now_tm = std::localtime(&now_t);

  // Create an output string stream
  std::ostringstream oss;

  // Write the date in year-month-day format to the stream
  oss << std::put_time(now_tm, "%Y-%m-%d");

  // Get the string from the stream
  std::string date_str = oss.str();

  // Return the string as TString
  return date_str;

}

/*
 * Do a Fourier fit for a deltaPhi distribution
 *
 *  TH1D* hDeltaPhi = Fitted deltaPhi histogram
 *  const int maxVn = Largest vn included in the Fourier fit
 *  const bool onlyNearSideFit = Do the fit only on the near side of the distribution
 *  const char* fitOption = Option given for the fitter
 *
 *  return = The Fourier fit function
 */
TF1* AlgorithmLibrary::FourierFit(TH1D* hDeltaPhi, const int maxVn, const bool onlyNearSideFit, const char* fitOption){
  
  // Define the fourier fit. Use fit parameters up to maxVn
  TString fourierFormula = "[0]*(1";
  for(int vn = 1; vn <= maxVn; vn++){
    fourierFormula = Form("%s+2.0*[%d]*TMath::Cos(%d.0*x)", fourierFormula.Data(), vn, vn);
  }
  fourierFormula = Form("%s)", fourierFormula.Data());
  
  TF1 *fourier = new TF1("fourier", fourierFormula.Data(), -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);

  // Set names for the parameters
  fourier->SetParName(0, "DistLevel");
  for(int vn = 1; vn <= maxVn; vn++){
    fourierFormula = Form("v%d",vn);
    fourier->SetParName(vn,fourierFormula.Data());
  }

  // Set initial values for the parameters
  fourier->SetParameter(0, hDeltaPhi->GetBinContent(hDeltaPhi->FindBin(TMath::Pi()/2)));
  fourier->SetParameter(1, 0);
  fourier->SetParameter(2, 0.03);
  for(int vn = 3; vn <= maxVn; vn++){
    fourier->SetParameter(vn, 0);
  }
  
  // Set limits such that the parameters must remain sensible
  fourier->SetParLimits(1, -1.0, 1.0);
  for(int vn = 2; vn <= maxVn; vn++){
    fourier->SetParLimits(vn, -1.0, 1.0);
  }
  
  // Do the fit!
  double lowerFitBound = -TMath::Pi()/2.0;
  double upperFitBound = 3.0*TMath::Pi()/2.0;
  if(onlyNearSideFit) upperFitBound = TMath::Pi()/2.0;
  hDeltaPhi->Fit("fourier", fitOption, "", lowerFitBound, upperFitBound);
  
  // Return the fitted distribution
  return fourier;
}

/*
 * Read a list of specifically formatted files
 *  Arguments: TString inputFileList = File containing a file name of the specifically formatted filei n each line
 *
 *  return: Tuple with vectors of files, legend comments and save comment
 */
std::tuple<std::vector<TFile*>, std::vector<TString>, std::vector<TString>, TString> AlgorithmLibrary::ReadFileList(TString inputFileList){

  // Define vectors for input files and legend string corresponding to said files
  std::vector<TFile*> inputFile;
  std::vector<TString> jetLegendString;
  std::vector<TString> saveNameString;
  TString saveComment;

  // Set up the input file list for reading
  std::ifstream file_stream(inputFileList);
  std::string line;
  TObjArray* lineContents;
  TObjString* lineItem;

  bool firstLine = true;

  // Open the input file list for reading
  if( file_stream.is_open() ) {
    
    // Loop over the lines in the input file list
    while( !file_stream.eof() ) {
      getline(file_stream, line);
      TString lineString(line);
      
      // If the line is non-empty, extract the file name and legend comment from it
      if( lineString.CompareTo("", TString::kExact) != 0 ) {

        if(firstLine){
          saveComment = lineString;
          firstLine = false;
        } else {
          // Other lines define the files that are compared, and a comment given to them in legend
          lineContents = lineString.Tokenize("&"); // Tokenize the string from '&' character

          // It is assumed that the line content before first '&' character gives the file name
          lineItem = (TObjString*)lineContents->At(0);
          inputFile.push_back(TFile::Open(lineItem->String()));

          // It is assumed that the line content between first and second '&' character gives the legend comment
          lineItem = (TObjString*)lineContents->At(1);
          jetLegendString.push_back(lineItem->String().Strip(TString::kBoth, ' '));

          // It is assumed that the line content after the second '&' character gives the save comment is only this figure is saved
          lineItem = (TObjString*)lineContents->At(2);
          saveNameString.push_back(lineItem->String().Strip(TString::kBoth, ' '));

          // Note: Strip command removes empty space from the beginning and ned of the string
        }
      } // Empty line if
      
    } // Loop over lines in the file
    
  // If cannot read the file, give error and end program
  } else {
    std::cout << "Error, could not open " << inputFileList.Data() << " for reading" << std::endl;
    std::cout << "Please check the file name!" << std::endl;
    inputFile.push_back(NULL);
    jetLegendString.push_back("ERROR");
    saveNameString.push_back("ERROR");
    saveComment = "ERROR";
  }

  // Collect the information to a tuple and return it
  return std::make_tuple(inputFile, jetLegendString, saveNameString, saveComment);

}
