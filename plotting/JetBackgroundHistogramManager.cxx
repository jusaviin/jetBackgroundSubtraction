/*
 * Implementation of JetBackgroundHistogramManager
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>

// Own includes
#include "JetBackgroundHistogramManager.h"

/*
 * Default constructor
 */
JetBackgroundHistogramManager::JetBackgroundHistogramManager() :
  fInputFile(NULL),
  fCard(NULL),
  fLoadEventInformation(false),
  fLoadJets(false),
  fLoad2DHistograms(false),
  fLoadJetPtClosureHistograms(false),
  fLoadJetPtResponseMatrix(false),
  fLoadJetEventPlaneCorrelationHistograms(false),
  fFirstLoadedCentralityBin(0),
  fLastLoadedCentralityBin(1),
  fFirstLoadedJetPtBin(0),
  fLastLoadedJetPtBin(1),
  fnCentralityBins(kMaxCentralityBins),
  fnJetPtBins(kMaxJetPtBins)
{
  
  // Default binning for centrality
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = iCentrality+1;
    fCentralityBinBorders[iCentrality] = 0;
  }
  
  // Default binning for jet pT in jet-event plane correlation histograms
  for(int iJetPt = 0; iJetPt < kMaxJetPtBins + 1; iJetPt++){
    fJetPtIndices[iJetPt] = iJetPt+1;
    fJetPtBinBorders[iJetPt] = 0;
  }
  
  // Initialize all the other histograms to null
  fhVertexZ = NULL;            // Vertex z position
  fhVertexZWeighted = NULL;    // Weighted vertex z-position (only meaningfull for MC)
  fhEvents = NULL;             // Number of events surviving different event cuts
  fhCentrality = NULL;         // Centrality of all events
  fhCentralityWeighted = NULL; // Weighted centrality distribution in all events (only meaningful for MC)
  fhPtHat = NULL;              // pT hat for MC events (only meaningful for MC)
  fhPtHatWeighted = NULL;      // Weighted pT hat distribution (only meaningful for MC)
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    
    // Jet histograms
    for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
      for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes+1; iParton++){
        for(int iMatch = 0; iMatch < JetBackgroundHistograms::knMatchingTypes+1; iMatch++){
          fhJetPt[iJetType][iCentrality][iParton][iMatch] = NULL;         // Jet pT histograms
          fhJetPhi[iJetType][iCentrality][iParton][iMatch] = NULL;       // Jet phi histograms
          fhJetEta[iJetType][iCentrality][iParton][iMatch] = NULL;       // Jet eta histograms
          fhJetEtaPhi[iJetType][iCentrality][iParton][iMatch] = NULL; // 2D eta-phi histogram for jets
        } // Matching jet loop
      } // Initiating parton type 
    } // Jet type loop
    
    // Jet pT closure histograms
    for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
      for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){
        for(int iJetPhi = 0; iJetPhi <= knJetPhiBins; iJetPhi++){
          for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes+1; iParton++){
            fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iParton] = NULL;
          } // Closure particle loop
        } // Jet phi bin loop 
      } // Jet eta bin loop
    } // Gen jet pT loop

    // Jet pT response matrix
    fhJetPtResponseMatrix[iCentrality] = NULL;

    // Jet-event plane correlation histograms
    for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
      for(int iJetPt = 0; iJetPt < kMaxJetPtBins; iJetPt++){
        for(int iOrder = 0; iOrder < JetBackgroundHistograms::knEventPlanes; iOrder++){
          fhJetEventPlane[iJetType][iOrder][iCentrality][iJetPt] = NULL;
        } // Event plane order loop
      } // Jet pT loop
    } // Jet type loop

  } // Centrality loop
}

/*
 * Constructor with input file
 */
JetBackgroundHistogramManager::JetBackgroundHistogramManager(TFile* inputFile) :
  JetBackgroundHistogramManager()
{
  fInputFile = inputFile;
  
  // Read card from inputfile
  fCard = new JetBackgroundCard(inputFile);
  
  // Initialize values using the information in card
  InitializeFromCard();
  
}

/*
 * Constructor with input file and input card
 */
JetBackgroundHistogramManager::JetBackgroundHistogramManager(TFile* inputFile, JetBackgroundCard* card) :
  JetBackgroundHistogramManager()
{
  fInputFile = inputFile;
  
  // Initialize values using the information in card
  fCard = card;
  InitializeFromCard();
  
}

/*
 * Constructor with input card
 */
JetBackgroundHistogramManager::JetBackgroundHistogramManager(JetBackgroundCard* card) :
  JetBackgroundHistogramManager()
{

  // Initialize values using the information in card
  fCard = card;
  InitializeFromCard();
  
}

/*
 * Initialize several member variables from JetBackgroundCard
 */
void JetBackgroundHistogramManager::InitializeFromCard(){
  
  // Read bins for centrality and jet pT from the card
  fnCentralityBins = fCard->GetNCentralityBins();
  fnJetPtBins = fCard->GetNJetPtBins();
  
  // Centrality binning
  for(int iCentrality = 0; iCentrality <= fnCentralityBins; iCentrality++){
    fCentralityBinBorders[iCentrality] = fCard->GetLowBinBorderCentrality(iCentrality);
  }
  fLastLoadedCentralityBin = fnCentralityBins-1;
  
  // Jet pT binning for energy-energy correlators
  for(int iJetPt = 0; iJetPt <= fnJetPtBins; iJetPt++){
    fJetPtBinBorders[iJetPt] = fCard->GetLowBinBorderJetPt(iJetPt);
  }
  fLastLoadedJetPtBin = fnJetPtBins-1;
  
}

/*
 * Copy constructor
 */
JetBackgroundHistogramManager::JetBackgroundHistogramManager(const JetBackgroundHistogramManager& in) :
  fInputFile(in.fInputFile),
  fCard(in.fCard),
  fLoadEventInformation(in.fLoadEventInformation),
  fLoadJets(in.fLoadJets),
  fLoad2DHistograms(in.fLoad2DHistograms),
  fLoadJetPtClosureHistograms(in.fLoadJetPtClosureHistograms),
  fLoadJetPtResponseMatrix(in.fLoadJetPtResponseMatrix),
  fLoadJetEventPlaneCorrelationHistograms(in.fLoadJetEventPlaneCorrelationHistograms),
  fFirstLoadedCentralityBin(in.fFirstLoadedCentralityBin),
  fLastLoadedCentralityBin(in.fLastLoadedCentralityBin),
  fFirstLoadedJetPtBin(in.fFirstLoadedJetPtBin),
  fLastLoadedJetPtBin(in.fLastLoadedJetPtBin),
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhEvents(in.fhEvents),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted)
{
  // Copy constructor
  
  // Copy binning for centrality
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = in.fCentralityBinIndices[iCentrality];
    fCentralityBinBorders[iCentrality] = in.fCentralityBinBorders[iCentrality];
  }

  // Copy binning for jet pT in jet-event plane histograms
  for(int iJetPt = 0; iJetPt < kMaxJetPtBins+1; iJetPt++){
    fJetPtIndices[iJetPt] = in.fJetPtIndices[iJetPt];
    fJetPtBinBorders[iJetPt] = in.fJetPtBinBorders[iJetPt];
  }
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    
    // Jet histograms
    for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
      for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes+1; iParton++){
        for(int iMatch = 0; iMatch < JetBackgroundHistograms::knMatchingTypes+1; iMatch++){
          fhJetPt[iJetType][iCentrality][iParton][iMatch] = in.fhJetPt[iJetType][iCentrality][iParton][iMatch];         // Jet pT histograms
          fhJetPhi[iJetType][iCentrality][iParton][iMatch] = in.fhJetPhi[iJetType][iCentrality][iParton][iMatch];       // Jet phi histograms
          fhJetEta[iJetType][iCentrality][iParton][iMatch] = in.fhJetEta[iJetType][iCentrality][iParton][iMatch];       // Jet eta histograms
          fhJetEtaPhi[iJetType][iCentrality][iParton][iMatch] = in.fhJetEtaPhi[iJetType][iCentrality][iParton][iMatch]; // 2D eta-phi histogram for jets
        } // Matching jet loop
      } // Initiating parton type 
    } // Jet type loop
    
    // Jet pT closure histograms
    for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
      for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){
        for(int iJetPhi = 0; iJetPhi <= knJetPhiBins; iJetPhi++){
          for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes+1; iParton++){
            fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iParton] = in.fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iParton];
          } // Closure particle loop
        } // Jet phi bin loop
      } // Jet eta bin loop
    } // Gen jet pT loop

    // Jet pT response matrix
    fhJetPtResponseMatrix[iCentrality] = in.fhJetPtResponseMatrix[iCentrality];

    // Jet-event plane correlation histograms
    for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
      for(int iJetPt = 0; iJetPt < kMaxJetPtBins; iJetPt++){
        for(int iOrder = 0; iOrder < JetBackgroundHistograms::knEventPlanes; iOrder++){
          fhJetEventPlane[iJetType][iOrder][iCentrality][iJetPt] = in.fhJetEventPlane[iJetType][iOrder][iCentrality][iJetPt];
        } // Event plane order loop
      } // Jet pT loop
    } // Jet type loop

  } // Centrality loop
}

/*
 * Destructor
 */
JetBackgroundHistogramManager::~JetBackgroundHistogramManager(){
  delete fCard;
}

/*
 * Load all the selected histograms from the inputfile
 */
void JetBackgroundHistogramManager::LoadHistograms(){
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                           // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                        // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");        // MC weighted vertex z position
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                  // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");  // MC weighted centrality in all events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                            // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");            // Weighted pT hat for MC events    
  }
  
  // Load jet histograms
  LoadJetHistograms();

  // Load jet pT closure histograms
  LoadJetPtClosureHistograms();

  // Load jet pT response matrices
  LoadJetPtResponseMatrix();

  // Load the jet-event plane correlation histograms
  LoadJetEventPlaneHistograms();
  
}

/*
 * Loader for jet histograms
 *
 * THnSparse for jets:
 *
 *   Histogram name: inclusiveJet/leadingJet
 *
 *     Axis index           Content of axis
 * --------------------------------------------------------
 *       Axis 0                 Jet pT
 *       Axis 1                 Jet phi
 *       Axis 2                 Jet eta
 *       Axis 3               Centrality
 *       Axis 4   Jet flavor: 1 = Quark, 2 = Gluon, 3 = Undetermined
 *       Axis 5       Matching jet exists: 1 = No, 2 = Yes
 */
void JetBackgroundHistogramManager::LoadJetHistograms(){

  // Only load the jet histograms is selected to do so
  if(!fLoadJets) return;
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  THnSparseD* histogramArray;
  
  // Define arrays to help find the histograms
  int axisIndices[3] = {0};
  int lowLimits[3] = {0};
  int highLimits[3] = {0};
  
  int nAxes = 1;           // Number of constraining axes for this iteration
  
  // Open the multidimensional histogram from which the histograms are projected
  for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
    histogramArray = (THnSparseD*) fInputFile->Get(fJetHistogramName[iJetType]);
  
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Reset the variable for axis constraints
      nAxes = 1;
    
      // Select the bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentrality];
      higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemoverCentrality;
    
      axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;  // Centrality

      // First, load the histograms without initiating parton or jet matching information
      fhJetPt[iJetType][iCentrality][JetBackgroundHistograms::knInitialPartonTypes][JetBackgroundHistograms::knMatchingTypes] = FindHistogram(histogramArray,0,nAxes,axisIndices,lowLimits,highLimits);
      fhJetPhi[iJetType][iCentrality][JetBackgroundHistograms::knInitialPartonTypes][JetBackgroundHistograms::knMatchingTypes] = FindHistogram(histogramArray,1,nAxes,axisIndices,lowLimits,highLimits);
      fhJetEta[iJetType][iCentrality][JetBackgroundHistograms::knInitialPartonTypes][JetBackgroundHistograms::knMatchingTypes] = FindHistogram(histogramArray,2,nAxes,axisIndices,lowLimits,highLimits);
      if(fLoad2DHistograms) fhJetEtaPhi[iJetType][iCentrality][JetBackgroundHistograms::knInitialPartonTypes][JetBackgroundHistograms::knMatchingTypes] = FindHistogram2D(histogramArray,1,2,nAxes,axisIndices,lowLimits,highLimits);

      // Next, load the histograms for different initiating partons, but without requirement for matching
      nAxes++;
      for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes; iParton++){

        axisIndices[1] = 4; lowLimits[1] = iParton+1; highLimits[1] = iParton+1;  // Initial parton type

        fhJetPt[iJetType][iCentrality][iParton][JetBackgroundHistograms::knMatchingTypes] = FindHistogram(histogramArray,0,nAxes,axisIndices,lowLimits,highLimits);
        fhJetPhi[iJetType][iCentrality][iParton][JetBackgroundHistograms::knMatchingTypes] = FindHistogram(histogramArray,1,nAxes,axisIndices,lowLimits,highLimits);
        fhJetEta[iJetType][iCentrality][iParton][JetBackgroundHistograms::knMatchingTypes] = FindHistogram(histogramArray,2,nAxes,axisIndices,lowLimits,highLimits);
        if(fLoad2DHistograms) fhJetEtaPhi[iJetType][iCentrality][iParton][JetBackgroundHistograms::knMatchingTypes] = FindHistogram2D(histogramArray,1,2,nAxes,axisIndices,lowLimits,highLimits);

        // Also load the histograms with matching/not matching jets
        nAxes++;
        for(int iMatch = 0; iMatch < JetBackgroundHistograms::knMatchingTypes; iMatch++){

          axisIndices[2] = 5; lowLimits[2] = iMatch+1; highLimits[2] = iMatch+1;  // Has matching jet

          fhJetPt[iJetType][iCentrality][iParton][iMatch] = FindHistogram(histogramArray,0,nAxes,axisIndices,lowLimits,highLimits);
          fhJetPhi[iJetType][iCentrality][iParton][iMatch] = FindHistogram(histogramArray,1,nAxes,axisIndices,lowLimits,highLimits);
          fhJetEta[iJetType][iCentrality][iParton][iMatch] = FindHistogram(histogramArray,2,nAxes,axisIndices,lowLimits,highLimits);
          if(fLoad2DHistograms) fhJetEtaPhi[iJetType][iCentrality][iParton][iMatch] = FindHistogram2D(histogramArray,1,2,nAxes,axisIndices,lowLimits,highLimits);

        } // Matching type loop
        histogramArray->GetAxis(5)->SetRange(0,0);
        nAxes--;

      } // Initiating parton loop

      // Finally, load histograms with or without jet matches, but disregarding initial parton type
      histogramArray->GetAxis(4)->SetRange(0,0);
    
      for(int iMatch = 0; iMatch < JetBackgroundHistograms::knMatchingTypes; iMatch++){

        axisIndices[1] = 5; lowLimits[1] = iMatch+1; highLimits[1] = iMatch+1;  // Has matching jet

        fhJetPt[iJetType][iCentrality][JetBackgroundHistograms::knInitialPartonTypes][iMatch] = FindHistogram(histogramArray,0,nAxes,axisIndices,lowLimits,highLimits);
        fhJetPhi[iJetType][iCentrality][JetBackgroundHistograms::knInitialPartonTypes][iMatch] = FindHistogram(histogramArray,1,nAxes,axisIndices,lowLimits,highLimits);
        fhJetEta[iJetType][iCentrality][JetBackgroundHistograms::knInitialPartonTypes][iMatch] = FindHistogram(histogramArray,2,nAxes,axisIndices,lowLimits,highLimits);
        if(fLoad2DHistograms) fhJetEtaPhi[iJetType][iCentrality][JetBackgroundHistograms::knInitialPartonTypes][iMatch] = FindHistogram2D(histogramArray,1,2,nAxes,axisIndices,lowLimits,highLimits);

      } // Matching type loop
    
    } // Centrality loop
  } // Jet type loop
}

/*
 * Loader for jet pT closure histograms
 *
 * THnSparse for closure histograms:
 *
 *   Histogram name: jetPtClosure
 *
 *     Axis index                  Content of axis
 * -----------------------------------------------------------
 *       Axis 0              Matched generator level jet pT
 *       Axis 1               Matched reconstructed jet pT
 *       Axis 2                         Jet eta
 *       Axis 3                       Centrality
 *       Axis 4                      Quark / gluon
 *       Axis 5             Matched reco to gen jet pT ratio
 *       Axis 6                         Jet phi
 */
void JetBackgroundHistogramManager::LoadJetPtClosureHistograms(){
  
  if(!fLoadJetPtClosureHistograms) return; // Do not load the histograms if they are not selected for loading
  
  // Define arrays to help find the histograms
  int axisIndices[4] = {0};
  int lowLimits[4] = {0};
  int highLimits[4] = {0};
  int nRestrictionAxes = 3;
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  THnSparseD* histogramArray = (THnSparseD*)fInputFile->Get("jetPtClosure");
  
  // Load all the histograms from the file
  for(int iGenJetPt = 0; iGenJetPt < knGenJetPtBins; iGenJetPt++){
    for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes+1; iParton++){
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Reset the ranges for all the axes in the histogram array
        for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
          histogramArray->GetAxis(iAxis)->SetRange(0,0);
        }

        // Select the bin indices
        lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
        higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
        
        // Setup the axes with restrictions
        nRestrictionAxes = 3;
        axisIndices[0] = 0; lowLimits[0] = iGenJetPt+1;    highLimits[0] = iGenJetPt+1;             // Gen jet pT
        axisIndices[1] = 3; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin; // Centrality
        axisIndices[2] = 4; lowLimits[2] = iParton+1; highLimits[2] = iParton+1;  // Quark/gluon
        
        // For the last closure particle bin no restrictions for quark/gluon jets
        if(iParton == JetBackgroundHistograms::knInitialPartonTypes){
          
          // Remove the closure particle requirement from the restriction axes
          nRestrictionAxes--;
        }
        
        fhJetPtClosure[iGenJetPt][knJetEtaBins][knJetPhiBins][iCentralityBin][iParton] = FindHistogram(histogramArray,5,nRestrictionAxes,axisIndices,lowLimits,highLimits);

        // Reset the range for the generator level jet pT axis
        histogramArray->GetAxis(0)->SetRange(0,0);
        
        // Eta binning for the closure histogram
        for(int iJetEta = 0; iJetEta < knJetEtaBins; iJetEta++){
          
          // For the last closure particle bin no restrictions for quark/gluon jets
          /*if(iParton == JetBackgroundHistograms::knInitialPartonTypes){
           nRestrictionAxes = 3;
           axisIndices[2] = 2; lowLimits[2] = iJetEta+1; highLimits[2] = iJetEta+1; // Jet eta
           } else {
           nRestrictionAxes = 4;
           axisIndices[3] = 2; lowLimits[3] = iJetEta+1; highLimits[3] = iJetEta+1; // Jet eta
           }
           
           fhJetPtClosure[iClosureType][iGenJetPt][iJetEta][knJetPhiBins][iCentralityBin][iParton] = FindHistogram(fInputFile,"jetPtClosure",5,nRestrictionAxes,axisIndices,lowLimits,highLimits);*/
          
          // Fill the pT integrated eta slices only once
          if(iGenJetPt == 0){
            
            // Setup the axes with restrictions
            nRestrictionAxes = 3;
            axisIndices[0] = 2; lowLimits[0] = iJetEta+1;    highLimits[0] = iJetEta+1;                 // Jet eta
            axisIndices[1] = 3; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin; // Centrality
            axisIndices[2] = 4; lowLimits[2] = iParton+1; highLimits[2] = iParton+1;  // Quark/gluon
            
            // For the last closure particle bin no restrictions for quark/gluon jets
            if(iParton == JetBackgroundHistograms::knInitialPartonTypes){
              
              // Remove the last set array bin
              nRestrictionAxes--;
            }
            
            fhJetPtClosure[knGenJetPtBins][iJetEta][knJetPhiBins][iCentralityBin][iParton] = FindHistogram(histogramArray,5,nRestrictionAxes,axisIndices,lowLimits,highLimits);
          }
          
        } // Jet eta bin loop

        // Reset the range for the jet eta axis
        histogramArray->GetAxis(2)->SetRange(0,0);

        // Phi binning for the closure histogram
        for(int iJetPhi = 0; iJetPhi < knJetPhiBins; iJetPhi++){
          
          // Fill the pT integrated phi slices only once
          if(iGenJetPt == 0){
            
            // Setup the axes with restrictions
            nRestrictionAxes = 3;
            axisIndices[0] = 6; lowLimits[0] = iJetPhi+1;    highLimits[0] = iJetPhi+1;                 // Jet phi
            axisIndices[1] = 3; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin; // Centrality
            axisIndices[2] = 4; lowLimits[2] = iParton+1; highLimits[2] = iParton+1;  // Quark/gluon
            
            // For the last closure particle bin no restrictions for quark/gluon jets
            if(iParton == JetBackgroundHistograms::knInitialPartonTypes){
              
              // Remove the last set array bin
              nRestrictionAxes--;
            }
            
            fhJetPtClosure[knGenJetPtBins][knJetEtaBins][iJetPhi][iCentralityBin][iParton] = FindHistogram(histogramArray,5,nRestrictionAxes,axisIndices,lowLimits,highLimits);
          }
          
        } // Jet eta bin loop
      } // Centrality loop
    } // Closure particle loop
  } // Gen jet pT loop
}

/*
 * Loader for jet pT response matrix
 *
 * THnSparse for closure histograms is used for this:
 *
 *   Histogram name: jetPtClosure
 *
 *     Axis index                  Content of axis
 * -----------------------------------------------------------
 *       Axis 0              Matched generator level jet pT
 *       Axis 1               Matched reconstructed jet pT
 *       Axis 2                         Jet eta
 *       Axis 3                       Centrality
 *       Axis 4                      Quark / gluon
 *       Axis 5             Matched reco to gen jet pT ratio
 *       Axis 6                         Jet phi
 */
void JetBackgroundHistogramManager::LoadJetPtResponseMatrix(){
  
  if(!fLoadJetPtResponseMatrix) return; // Do not load the histograms if they are not selected for loading
  
  // Define arrays to help find the histograms
  int axisIndices[1] = {0};
  int lowLimits[1] = {0};
  int highLimits[1] = {0};
  int nRestrictionAxes = 1;
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;

  // Find the histogram array from which the projections are made
  THnSparseD* histogramArray = (THnSparseD*)fInputFile->Get("jetPtClosure");

  // Load all the histograms from the file
  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
    // Select the bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentrality];
    higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemoverCentrality;
        
    // Setup centrality axis restrictions
    axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin; // Centrality

    // Project the response matrix from the closure histogram
    fhJetPtResponseMatrix[iCentrality] = FindHistogram2D(histogramArray,1,0,nRestrictionAxes,axisIndices,lowLimits,highLimits,false);
        
  } // Centrality loop
}

/*
 * Loader for jet-event plane correlation histograms
 *
 * THnSparse for jet-event plane correlation:
 *
 *   Histogram name: inclusiveJetEventPlaneOrder2/inclusiveJetEventPlaneOrder3/inclusiveJetEventPlaneOrder4
 *                    leadingJetEventPlaneOrder2 / leadingJetEventPlaneOrder3 / leadingJetEventPlaneOrder4
 *
 *     Axis index          Content of axis
 * --------------------------------------------------------
 *       Axis 0        Jet-event plane correlation
 *       Axis 1                 Jet pT
 *       Axis 2               Centrality
 */
void JetBackgroundHistogramManager::LoadJetEventPlaneHistograms(){

  // Only load the jet histograms is selected to do so
  if(!fLoadJetEventPlaneCorrelationHistograms) return;
  
  // Define helper variables
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerJetPtBin = 0;
  int higherJetPtBin = 0;
  THnSparseD* histogramArray;
  
  // Define arrays to help find the histograms
  int axisIndices[2] = {0};
  int lowLimits[2] = {0};
  int highLimits[2] = {0};
  
  int nAxes = 1;           // Number of constraining axes for this iteration
  
  // Open the multidimensional histogram from which the histograms are projected
  for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
    for(int iOrder = 0; iOrder < JetBackgroundHistograms::knEventPlanes; iOrder++){
      histogramArray = (THnSparseD*) fInputFile->Get(Form("%sEventPlaneOrder%d", fJetHistogramName[iJetType], iOrder+2));
  
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

        // Reset the variable for axis constraints
        nAxes = 1;
    
        // Select the bin indices
        lowerCentralityBin = fCentralityBinIndices[iCentrality];
        higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemover;
    
        axisIndices[0] = 2; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;  // Centrality

        // First, load the histograms without jet pT requirements
        fhJetEventPlane[iJetType][iOrder][iCentrality][fnJetPtBins] = FindHistogram(histogramArray,0,nAxes,axisIndices,lowLimits,highLimits);

        // Next, load the histograms for different jet pT bins
        nAxes++;
        for(int iJetPt = fFirstLoadedJetPtBin; iJetPt <= fLastLoadedJetPtBin; iJetPt++){

          // Select the bin indices
          lowerJetPtBin = fJetPtIndices[iJetPt];
          higherJetPtBin = fJetPtIndices[iJetPt+1]+duplicateRemover;

          axisIndices[1] = 1; lowLimits[1] = lowerJetPtBin; highLimits[1] = higherJetPtBin;  // Jet pT

          // Load the histograms with jet pT binning
          fhJetEventPlane[iJetType][iOrder][iCentrality][iJetPt] = FindHistogram(histogramArray,0,nAxes,axisIndices,lowLimits,highLimits);

        } // Jet pT loop
        histogramArray->GetAxis(1)->SetRange(0,0);
      } // Centrality loop
    } // Event plane order loop
  } // Jet type loop
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   THnSparseD* histogramArray = Inputfile containing the THnSparse to be read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int nAxes = Number of axes that are restained for the projection
 *   int* axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int* lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int* highBinIndex = Indices of the highest considered bins in the restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH2D* JetBackgroundHistogramManager::FindHistogram2D(THnSparseD* histogramArray, int xAxis, int yAxis, int nAxes, int* axisNumber, int* lowBinIndex, int* highBinIndex, const bool normalizeToBinWidth){
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  TString newName = histogramArray->GetName();
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    newName.Append(Form("_%d=%d-%d",axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]));
  }
  
  // Project out the histogram and give it the created unique name
  TH2D* projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName.Data());
  
  // Apply bin width normalization to the projected histogram
  if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   THnSparseD* histogramArray = Inputfile containing the THnSparse to be read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH2D* JetBackgroundHistogramManager::FindHistogram2D(THnSparseD* histogramArray, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2, const bool normalizeToBinWidth){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram2D(histogramArray,xAxis,yAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices,normalizeToBinWidth);
}

/*
 * Extract a histogram with given restrictions on other axes in THnSparse
 *
 *  Arguments:
 *   THnSparseD* histogramArray = Histogram array from which the desired histograms are projected
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   int nAxes = Number of axes that are restained for the projection
 *   int* axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int* lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int* highBinIndex = Indices of the highest considered bins in the restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH1D* JetBackgroundHistogramManager::FindHistogram(THnSparseD* histogramArray, int xAxis, int nAxes, int* axisNumber, int* lowBinIndex, int* highBinIndex, const bool normalizeToBinWidth){
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for each histogram that is read from the file
  TString newName = histogramArray->GetName();
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    newName.Append(Form("_%d=%d-%d",axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]));
  }
  
  // Project out the histogram and give it the created unique name
  TH1D* projectedHistogram = NULL;
  
  // Check that we are not trying to project a non-existing axis
  if(xAxis < histogramArray->GetNdimensions()){
    projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName.Data());
  
    // Apply bin width normalization to the projected histogram
    if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  }
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   THnSparseD* histogramArray = Histogram array from which the desired histograms are projected
 *   int xAxis = Index for the axis in THnSparse that is projected to TH1D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH1D* JetBackgroundHistogramManager::FindHistogram(THnSparseD* histogramArray, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2, const bool normalizeToBinWidth){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram(histogramArray,xAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices,normalizeToBinWidth);
}

/*
 * Write all the loaded histograms into a file
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void JetBackgroundHistogramManager::Write(const char* fileName, const char* fileOption){
  
  // Create the output file
  TFile* outputFile = new TFile(fileName,fileOption);
  
  // Helper variable for renaming the saved histograms
  TString histogramNamer;
  
  // Write the event information histograms to the output file
  if(fLoadEventInformation){
    fhEvents->Write("",TObject::kOverwrite);             // Number of events surviving different event cuts
    fhVertexZ->Write("",TObject::kOverwrite);            // Vertex z position
    fhVertexZWeighted->Write("",TObject::kOverwrite);    // MC weighted vertex z position
    fhCentrality->Write("",TObject::kOverwrite);         // Centrality in all events
    fhCentralityWeighted->Write("",TObject::kOverwrite); // MC weighted centrality in all events
    fhPtHat->Write("",TObject::kOverwrite);              // pT hat for MC events (only meaningful for MC)
    fhPtHatWeighted->Write("",TObject::kOverwrite);      // Weighted pT hat distribution (only meaningful for MC)
  }
 
  // Write the jet histograms to the output file
  WriteJetHistograms();

  // Write the jet pT closure histograms to the output file
  WriteClosureHistograms();

  // Write the jet pT response matrices to the output file
  WriteJetPtResponseMatrix(); 

  // Write the jet-event plane correlation histograms to the output file
  WriteJetEventPlaneHistograms();
  
  // Write the card to the output file if it is not already written
  if(!gDirectory->GetDirectory("JCard")) fCard->Write(outputFile);
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
}

/*
 * Write the jet histograms to the file that is currently open
 */
void JetBackgroundHistogramManager::WriteJetHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Write the jet histograms to the output file
  if(!fLoadJets) return;  // Only write the jet histograms if they are loaded

  for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
  
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fJetHistogramName[iJetType])) gDirectory->mkdir(fJetHistogramName[iJetType]);
    gDirectory->cd(fJetHistogramName[iJetType]);
  
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes+1; iParton++){
        for(int iMatch = 0; iMatch < JetBackgroundHistograms::knMatchingTypes+1; iMatch++){
      
          // Jet pT
          histogramNamer = Form("%sPt%s%s_C%d", fJetHistogramName[iJetType], fInitialPartonName[iParton], fMatchingName[iMatch], iCentrality);
          if(fhJetPt[iJetType][iCentrality][iParton][iMatch]) fhJetPt[iJetType][iCentrality][iParton][iMatch]->Write(histogramNamer.Data(), TObject::kOverwrite);
    
          // Jet phi
          histogramNamer = Form("%sPhi%s%s_C%d", fJetHistogramName[iJetType], fInitialPartonName[iParton], fMatchingName[iMatch], iCentrality);
          if(fhJetPhi[iJetType][iCentrality][iParton][iMatch]) fhJetPhi[iJetType][iCentrality][iParton][iMatch]->Write(histogramNamer.Data(), TObject::kOverwrite);
    
          // Jet eta
          histogramNamer = Form("%sEta%s%s_C%d", fJetHistogramName[iJetType], fInitialPartonName[iParton], fMatchingName[iMatch], iCentrality);
          if(fhJetEta[iJetType][iCentrality][iParton][iMatch]) fhJetEta[iJetType][iCentrality][iParton][iMatch]->Write(histogramNamer.Data(), TObject::kOverwrite);
    
          // Jet eta-phi
          histogramNamer = Form("%sEtaPhi%s%s_C%d", fJetHistogramName[iJetType], fInitialPartonName[iParton], fMatchingName[iMatch], iCentrality);
          if(fLoad2DHistograms && fhJetEtaPhi[iJetType][iCentrality][iParton][iMatch]) fhJetEtaPhi[iJetType][iCentrality][iParton][iMatch]->Write(histogramNamer.Data(), TObject::kOverwrite);
        } // Matching type loop
      } // Initiating parton type loop  
    } // Centrality loop
  
    // Return back to main directory
    gDirectory->cd("../");

  } // Jet type loop
  
}

/*
 * Write the closure histograms to the file that is currently open
 */
void JetBackgroundHistogramManager::WriteClosureHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Only write the jet pT closure histograms if they are previously loaded
  if(fLoadJetPtClosureHistograms){
    
    // Create a directory for the histograms if it does not already exist
    histogramNamer = Form("jetPtClosure_%s",fJetHistogramName[0]);
    if(!gDirectory->GetDirectory(histogramNamer.Data())) gDirectory->mkdir(histogramNamer.Data());
    gDirectory->cd(histogramNamer.Data());
    
    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Loop over closure particles (quark/gluon/no selection)
      for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes+1; iParton++){
        
        // Loop over generator level jet pT bins
        for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
          
          // Loop over jet eta bins
          for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){

            // Loop over jet phi bins
            for(int iJetPhi = 0; iJetPhi <= knJetPhiBins; iJetPhi++){
            
              // Only write histogram that are non-NULL
              if(fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iParton]){
                histogramNamer = Form("jetPtClosure_%s%s_C%d", fJetHistogramName[0], fInitialPartonName[iParton], iCentrality);
                if(iGenJetPt < knGenJetPtBins) histogramNamer.Append(Form("T%d",iGenJetPt));
                if(iJetEta < knJetEtaBins) histogramNamer.Append(Form("E%d",iJetEta));
                if(iJetPhi < knJetPhiBins) histogramNamer.Append(Form("P%d",iJetPhi));
                fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iParton]->Write(histogramNamer.Data(), TObject::kOverwrite);
              }
            
            } // JEt phi bin loop
          } // Jet eta bin loop
        } // Generator level jet pT loop
      } // Closure particle type (quark/gluon) loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Writing jet pT closure histograms
  
}

/*
 * Write the jet pT response matrices to the file that is currently open
 */
void JetBackgroundHistogramManager::WriteJetPtResponseMatrix(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Only write the jet pT closure histograms if they are previously loaded
  if(fLoadJetPtResponseMatrix){
    
    // Create a directory for the histograms if it does not already exist
    histogramNamer = "jetPtResponseMatrix";
    if(!gDirectory->GetDirectory(histogramNamer.Data())) gDirectory->mkdir(histogramNamer.Data());
    gDirectory->cd(histogramNamer.Data());

    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Only write histograms that are non-NULL
      if(fhJetPtResponseMatrix[iCentrality]) {
        histogramNamer = Form("jetPtResponseMatrix_C%d", iCentrality);
        fhJetPtResponseMatrix[iCentrality]->Write(histogramNamer.Data(), TObject::kOverwrite);
      }

    }  // Centrality loop

    // Return back to main directory
    gDirectory->cd("../");
    
  } // Writing jet pT response matrices
  
}

/*
 * Write the jet-event plane correlation histograms to the file that is currently open
 */
void JetBackgroundHistogramManager::WriteJetEventPlaneHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Write the jet histograms to the output file
  if(!fLoadJetEventPlaneCorrelationHistograms) return;  // Only write the histograms if they are loaded

  for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
    for(int iOrder = 0; iOrder < JetBackgroundHistograms::knEventPlanes; iOrder++){
  
      // Create a directory for the histograms if it does not already exist
      histogramNamer = Form("%sEventPlaneOrder%d", fJetHistogramName[iJetType], iOrder+2);
      if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
      gDirectory->cd(histogramNamer);
  
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

        // First write the histograms without jet pT selection
        histogramNamer = Form("%sEventPlaneOrder%d_C%d", fJetHistogramName[iJetType], iOrder+2, iCentrality);
        if(fhJetEventPlane[iJetType][iOrder][iCentrality][fnJetPtBins]) fhJetEventPlane[iJetType][iOrder][iCentrality][fnJetPtBins]->Write(histogramNamer.Data(), TObject::kOverwrite);

        for(int iJetPt = fFirstLoadedJetPtBin; iJetPt <= fLastLoadedJetPtBin; iJetPt++){

          // Then write the histograms with jet pT selection
          histogramNamer = Form("%sEventPlaneOrder%d_C%dT%d", fJetHistogramName[iJetType], iOrder+2, iCentrality, iJetPt);
          if(fhJetEventPlane[iJetType][iOrder][iCentrality][iJetPt]) fhJetEventPlane[iJetType][iOrder][iCentrality][iJetPt]->Write(histogramNamer.Data(), TObject::kOverwrite);

        } // Jet pT loop
      } // Centrality loop
  
      // Return back to main directory
      gDirectory->cd("../");
    } // Event plane order loop
  } // Jet type loop
  
}

/*
 * Load the selected histograms from a file containing readily processed histograms
 */
void JetBackgroundHistogramManager::LoadProcessedHistograms(){
  
  // Helper variable for finding names of loaded histograms
  TString histogramNamer;
  TString folderNamer;
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                           // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                        // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");        // MC weighted vertex z position
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                  // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");  // MC weighted centrality in all events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                            // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");            // Weighted pT hat for MC events
  }
  
  // Load the jet histograms from the processed file
  if(fLoadJets){
    for(int iJetType = 0; iJetType < knJetTypes; iJetType++){  
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes+1; iParton++){
          for(int iMatch = 0; iMatch < JetBackgroundHistograms::knMatchingTypes+1; iMatch++){

            // Load jet pT histograms
            histogramNamer = Form("%s/%sPt%s%s_C%d", fJetHistogramName[iJetType], fJetHistogramName[iJetType], fInitialPartonName[iParton], fMatchingName[iMatch], iCentrality);
            fhJetPt[iJetType][iCentrality][iParton][iMatch] = (TH1D*) fInputFile->Get(histogramNamer.Data());
    
            // Load jet phi histograms
            histogramNamer = Form("%s/%sPhi%s%s_C%d", fJetHistogramName[iJetType], fJetHistogramName[iJetType], fInitialPartonName[iParton], fMatchingName[iMatch], iCentrality);
            fhJetPhi[iJetType][iCentrality][iParton][iMatch] = (TH1D*) fInputFile->Get(histogramNamer.Data());
    
            // Load jet eta histograms
            histogramNamer = Form("%s/%sEta%s%s_C%d", fJetHistogramName[iJetType], fJetHistogramName[iJetType], fInitialPartonName[iParton], fMatchingName[iMatch], iCentrality);
            fhJetEta[iJetType][iCentrality][iParton][iMatch] = (TH1D*) fInputFile->Get(histogramNamer.Data());
    
            // Load jet eta-phi histograms
            histogramNamer = Form("%s/%sEtaPhi%s%s_C%d", fJetHistogramName[iJetType], fJetHistogramName[iJetType], fInitialPartonName[iParton], fMatchingName[iMatch], iCentrality);
            if(fLoad2DHistograms){
             fhJetEtaPhi[iJetType][iCentrality][iParton][iMatch] = (TH2D*) fInputFile->Get(histogramNamer.Data());
            } // Loading 2D histograms
          } // Matching status loop
        } // Initiating parton type loop
      } // Centrality loop
    } // Jet type loop
  } // Loading jet histograms

  // Load the jet pT closure histograms from a processed file
  if(fLoadJetPtClosureHistograms){
    
    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Loop over closure particles (quark/gluon/no selection)
      for(int iParton = 0; iParton < JetBackgroundHistograms::knInitialPartonTypes+1; iParton++){
        
        // Loop over generator level jet pT bins
        for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
          
          // Loop over jet eta bins
          for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){

            // Loop over jet phi bins
            for(int iJetPhi = 0; iJetPhi <= knJetPhiBins; iJetPhi++){
            
              histogramNamer = Form("jetPtClosure_%s/jetPtClosure_%s%s_C%d", fJetHistogramName[0], fJetHistogramName[0], fInitialPartonName[iParton], iCentrality);
              if(iGenJetPt < knGenJetPtBins) histogramNamer.Append(Form("T%d",iGenJetPt));
              if(iJetEta < knJetEtaBins) histogramNamer.Append(Form("E%d",iJetEta));
              if(iJetPhi < knJetPhiBins) histogramNamer.Append(Form("P%d",iJetPhi));
              fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iParton] = (TH1D*) fInputFile->Get(histogramNamer.Data());
            
            } // Jet phi bin loop
          } // Jet eta bin loop
        } // Generator level jet pT loop
      } // Closure particle type (quark/gluon) loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Opening jet pT closure histograms

  // Load the jet pT response matrices from a processed file
  if(fLoadJetPtResponseMatrix){

    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      histogramNamer = Form("jetPtResponseMatrix/jetPtResponseMatrix_C%d", iCentrality);
      fhJetPtResponseMatrix[iCentrality] = (TH2D*)fInputFile->Get(histogramNamer.Data());

    }  // Centrality loop
  } // Loading jet pT response matrices

  // Load the jet-event plane correlation histograms from the processed file
  if(fLoadJetEventPlaneCorrelationHistograms){
    for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
      for(int iOrder = 0; iOrder < JetBackgroundHistograms::knEventPlanes; iOrder++){
  
        // There are different folders for each jet type and each event plane order
        folderNamer = Form("%sEventPlaneOrder%d", fJetHistogramName[iJetType], iOrder+2);

        for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

          // Load the histograms without jet pT restrictions
          histogramNamer = Form("%s/%s_C%d", folderNamer.Data(), folderNamer.Data(), iCentrality);
          fhJetEventPlane[iJetType][iOrder][iCentrality][fnJetPtBins] = (TH1D*) fInputFile->Get(histogramNamer.Data());

          for(int iJetPt = fFirstLoadedJetPtBin; iJetPt <= fLastLoadedJetPtBin; iJetPt++){

            // Load the histograms without jet pT restrictions
            histogramNamer = Form("%s/%s_C%dT%d", folderNamer.Data(), folderNamer.Data(), iCentrality, iJetPt);
            fhJetEventPlane[iJetType][iOrder][iCentrality][iJetPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());

          } // Jet pT loop
        } // Centrality loop
      } // Event plane order loop
    } // Jet type loop
  } // Loading jet-event plane correlation histograms

}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const char* histogramName = Name of the histogram from which the bin indices are searched
 *   const int nBins = Number of bins for the indices
 *   int* binIndices = Array of integers to be filled with bin index information read from the file
 *   const double* binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 */
void JetBackgroundHistogramManager::SetBinIndices(const char* histogramName, const int nBins, int* binIndices, const double* binBorders, const int iAxis){
  THnSparseD* histogramArray = (THnSparseD*) fInputFile->Get(histogramName);
  TH1D* hBinner = FindHistogram(histogramArray,iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const char* histogramName = Name of the histogram from which the bin indices are searched
 *   const int nBins = Number of bins for the indices
 *   double* copyBinBorders = Array to which a copy of bin borders is made
 *   int* binIndices = Array of integers to be filled with bin index information read from the file
 *   const double* binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 *   const bool setIndices = true: Set both bin indices and bin borders, false: Set only bin borders
 */
void JetBackgroundHistogramManager::SetBinBordersAndIndices(const char* histogramName, const int nBins, double* copyBinBorders, int* binIndices, const double* binBorders, const int iAxis, const bool setIndices){
  TH1D* hBinner;
  THnSparseD* histogramArray;
  if(setIndices) {
    histogramArray = (THnSparseD*) fInputFile->Get(histogramName);
    hBinner = FindHistogram(histogramArray,iAxis,0,0,0);
  }
  for(int iBin = 0; iBin < nBins+1; iBin++){
    copyBinBorders[iBin] = binBorders[iBin];
    if(setIndices) binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Set up generic bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const char* histogramName = Name of the histogram from which the bin indices are searched
 *  const int iAxis = Axis from which the set indices can be found
 *  int nSetBins = Number of bins that is set
 *  double* setBinBorders = Bin borders that are set
 *  int* setBinIndices = Bin indices that are set
 *  const int nBins = New number of bins that is given
 *  const double* binBorders = New bin borders that are given
 *  const char* errorMessage = Type of the set bins to be printed in possible error message
 *  const int maxBins = Maximum number of allowed bins of this type
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void JetBackgroundHistogramManager::SetGenericBins(const bool readBinsFromFile, const char* histogramName, const int iAxis, int nSetBins, double* setBinBorders, int* setBinIndices, const int nBins, const double* binBorders, const char* errorMessage, const int maxBins, const bool setIndices){
  
  // If bins are read from file, do not use the given bin borders
  if(readBinsFromFile){
    if(setIndices) SetBinIndices(histogramName, nSetBins, setBinIndices, setBinBorders, iAxis);
  } else { // If the given bin borders are use, update the number of bins and bin borders according to input
    if(nBins <= maxBins){
      nSetBins = nBins;
      SetBinBordersAndIndices(histogramName, nSetBins, setBinBorders, setBinIndices, binBorders, iAxis, setIndices);
    } else {
      cout << "Error! Too many " << errorMessage << " bins given. Maximum number is " << maxBins << ". Will not set bins." << endl;
    }
  }
}

/*
 * Set up centrality bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given centrality bins
 *  const double* binBorders = New bin borders for centrality
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void JetBackgroundHistogramManager::SetCentralityBins(const bool readBinsFromFile, const int nBins, const double* binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, "inclusiveJetEventPlaneOrder2", 2, fnCentralityBins, fCentralityBinBorders, fCentralityBinIndices, nBins, binBorders, "centrality", kMaxCentralityBins, setIndices);
  
}


/*
 * Set up jet pT bin indices for energy-energy correlator according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double* binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void JetBackgroundHistogramManager::SetJetPtBins(const bool readBinsFromFile, const int nBins, const double* binBorders, bool setIndices){
  
  SetGenericBins(readBinsFromFile, "inclusiveJetEventPlaneOrder2", 1, fnJetPtBins, fJetPtBinBorders, fJetPtIndices, nBins, binBorders, "jet pT", kMaxJetPtBins, setIndices);
  
}

// Setter for loading event information
void JetBackgroundHistogramManager::SetLoadEventInformation(const bool loadOrNot){
  fLoadEventInformation = loadOrNot;
}

// Setter for loading jet histograms
void JetBackgroundHistogramManager::SetLoadJetHistograms(const bool loadOrNot){
  fLoadJets = loadOrNot;
}

// Setter for loading two-dimensional histograms
void JetBackgroundHistogramManager::SetLoad2DHistograms(const bool loadOrNot){
  fLoad2DHistograms = loadOrNot;
}

// Setter for loading jet pT closure histograms
void JetBackgroundHistogramManager::SetLoadJetPtClosureHistograms(const bool loadOrNot){
  fLoadJetPtClosureHistograms = loadOrNot;
}

// Setter for loading jet pT response matrix
void JetBackgroundHistogramManager::SetLoadJetPtResponseMatrix(const bool loadOrNot){
  fLoadJetPtResponseMatrix = loadOrNot;
}

// Setter for loading jet-event plane correlation histograms
void JetBackgroundHistogramManager::SetLoadJetEventPlaneHistograms(const bool loadOrNot){
  fLoadJetEventPlaneCorrelationHistograms = loadOrNot;
}


// Setter for loaded centrality bins
void JetBackgroundHistogramManager::SetCentralityBinRange(const int first, const int last){
  fFirstLoadedCentralityBin = first;
  fLastLoadedCentralityBin = last;
  
  // Sanity check for centrality bins
  BinSanityCheck(fnCentralityBins,fFirstLoadedCentralityBin,fLastLoadedCentralityBin);
}

// Setter for jet pT bin range in jet-event plane correlation histograms
void JetBackgroundHistogramManager::SetJetPtBinRange(const int first, const int last){
  fFirstLoadedJetPtBin = first;
  fLastLoadedJetPtBin = last;
  
  // Sanity check for jet pT bins in energy-energy correlator histograms
  BinSanityCheck(fnJetPtBins,fFirstLoadedJetPtBin,fLastLoadedJetPtBin);
}

// Sanity check for set bins
void JetBackgroundHistogramManager::BinSanityCheck(const int nBins, int& first, int& last){
  if(first < 0) first = 0;
  if(last < first) last = first;
  if(last > nBins-1) last = nBins-1;
}

// Sanity check for input bin index
int JetBackgroundHistogramManager::BinIndexCheck(const int nBins, const int binIndex) const{
  if(binIndex < 0) return 0;
  if(binIndex > nBins-1) return nBins-1;
  return binIndex;
}

// Getter for the number of centrality bins
int JetBackgroundHistogramManager::GetNCentralityBins() const{
  return fnCentralityBins;
}

// Getter for the number of jet pT bins in jet-event plane histograms
int JetBackgroundHistogramManager::GetNJetPtBins() const{
  return fnJetPtBins;
}

// Getter for the jet histogram name
const char* JetBackgroundHistogramManager::GetJetHistogramName(const int iJetType) const{
  return fJetHistogramName[iJetType];
}

// Getter for name suitable for x-axis in a given jet histogram
const char* JetBackgroundHistogramManager::GetJetAxisName(const int iJetType) const{
  return fJetAxisName[iJetType];
}

// Getter for i:th centrality bin border
double JetBackgroundHistogramManager::GetCentralityBinBorder(const int iCentrality) const{
  return fCentralityBinBorders[iCentrality];
}

// Getter for i:th jet pT bin border in energy-energy correlator histograms
double JetBackgroundHistogramManager::GetJetPtBinBorder(const int iJetPt) const{
  return fJetPtBinBorders[iJetPt];
}

// Getters for event information histograms

// Getter for z-vertex histogram
TH1D* JetBackgroundHistogramManager::GetHistogramVertexZ() const{
  return fhVertexZ;
}

// Getter for z-vertex histogram
TH1D* JetBackgroundHistogramManager::GetHistogramVertexZWeighted() const{
  return fhVertexZWeighted;
}

// Getter for histogram for number of events surviving different event cuts
TH1D* JetBackgroundHistogramManager::GetHistogramEvents() const{
  return fhEvents;
}

// Getter for centrality histogram in all events
TH1D* JetBackgroundHistogramManager::GetHistogramCentrality() const{
  return fhCentrality;
}

// Getter for weighted centrality histogram in all events
TH1D* JetBackgroundHistogramManager::GetHistogramCentralityWeighted() const{
  return fhCentralityWeighted;
}

// Getters for jet histograms

// Getter for jet pT histograms
TH1D* JetBackgroundHistogramManager::GetHistogramJetPt(int iCentrality, int iJetType, int iParton, int iMatch) const{
  return fhJetPt[iJetType][iCentrality][iParton][iMatch];
}

// Getter for inclusive jet pT histograms
TH1D* JetBackgroundHistogramManager::GetHistogramInclusiveJetPt(int iCentrality, int iParton, int iMatch) const{
  return GetHistogramJetPt(iCentrality, kInclusiveJet, iParton, iMatch);
}

// Getter for leading jet pT histograms
TH1D* JetBackgroundHistogramManager::GetHistogramLeadingJetPt(int iCentrality, int iParton, int iMatch) const{
  return GetHistogramJetPt(iCentrality, kLeadingJet, iParton, iMatch);
}

// Getter for jet phi histograms
TH1D* JetBackgroundHistogramManager::GetHistogramJetPhi(int iCentrality, int iJetType, int iParton, int iMatch) const{
  return fhJetPhi[iJetType][iCentrality][iParton][iMatch];
}

// Getter for inclusive jet phi histograms
TH1D* JetBackgroundHistogramManager::GetHistogramInclusiveJetPhi(int iCentrality, int iParton, int iMatch) const{
  return GetHistogramJetPhi(iCentrality, kInclusiveJet, iParton, iMatch);
}

// Getter for leading jet phi histograms
TH1D* JetBackgroundHistogramManager::GetHistogramLeadingJetPhi(int iCentrality, int iParton, int iMatch) const{
  return GetHistogramJetPhi(iCentrality, kLeadingJet, iParton, iMatch);
}

// Getter for jet eta histograms
TH1D* JetBackgroundHistogramManager::GetHistogramJetEta(int iCentrality, int iJetType, int iParton, int iMatch) const{
  return fhJetEta[iJetType][iCentrality][iParton][iMatch];
}

// Getter for inclusive jet eta histograms
TH1D* JetBackgroundHistogramManager::GetHistogramInclusiveJetEta(int iCentrality, int iParton, int iMatch) const{
  return GetHistogramJetEta(iCentrality, kInclusiveJet, iParton, iMatch);
}

// Getter for leading jet eta histograms
TH1D* JetBackgroundHistogramManager::GetHistogramLeadingJetEta(int iCentrality, int iParton, int iMatch) const{
  return GetHistogramJetEta(iCentrality, kLeadingJet, iParton, iMatch);
}

// Getter for 2D eta-phi histogram for jets
TH2D* JetBackgroundHistogramManager::GetHistogramJetEtaPhi(int iCentrality, int iJetType, int iParton, int iMatch) const{
  return fhJetEtaPhi[iJetType][iCentrality][iParton][iMatch];
}

// Getter for inclusive jet phi histograms
TH2D* JetBackgroundHistogramManager::GetHistogramInclusiveJetEtaPhi(int iCentrality, int iParton, int iMatch) const{
  return GetHistogramJetEtaPhi(iCentrality, kInclusiveJet, iParton, iMatch);
}

// Getter for leading jet phi histograms
TH2D* JetBackgroundHistogramManager::GetHistogramLeadingJetEtaPhi(int iCentrality, int iParton, int iMatch) const{
  return GetHistogramJetEtaPhi(iCentrality, kLeadingJet, iParton, iMatch);
}

// Getter for jet pT closure histograms
TH1D* JetBackgroundHistogramManager::GetHistogramJetPtClosure(const int iGenPtBin, const int iEtaBin, const int iPhiBin, const int iCentrality, const int iParton) const{
  return fhJetPtClosure[iGenPtBin][iEtaBin][iPhiBin][iCentrality][iParton];
}

// Getter for jet pT response matrix
TH2D* JetBackgroundHistogramManager::GetHistogramJetPtResponseMatrix(const int iCentrality) const{
  return fhJetPtResponseMatrix[iCentrality];
}

// Getter for jet-event plane histograms
TH1D* JetBackgroundHistogramManager::GetHistogramJetEventPlane(int iOrder, int iJetType, int iCentrality, int iJetPt){
  if(iJetPt < 0) iJetPt = fnJetPtBins;
  return fhJetEventPlane[iJetType][iOrder-2][iCentrality][iJetPt];
}

// Getter for inclusive jet-event plane histograms
TH1D* JetBackgroundHistogramManager::GetHistogramInclusiveJetEventPlane(int iOrder, int iCentrality, int iJetPt){
  return GetHistogramJetEventPlane(iOrder, kInclusiveJet, iCentrality, iJetPt);
}


// Getter for leading jet-event plane histograms
TH1D* JetBackgroundHistogramManager::GetHistogramLeadingJetEventPlane(int iOrder, int iCentrality, int iJetPt){
  return GetHistogramJetEventPlane(iOrder, kLeadingJet, iCentrality, iJetPt);
}


// Get the first loaded centrality bin
int JetBackgroundHistogramManager::GetFirstCentralityBin() const{
  return fFirstLoadedCentralityBin;
}

// Get the last loaded centrality bin
int JetBackgroundHistogramManager::GetLastCentralityBin() const{
  return fLastLoadedCentralityBin;
}

// Get the first loaded jet pT bin for jet-event plane histograms
int JetBackgroundHistogramManager::GetFirstJetPtBin() const{
  return fFirstLoadedJetPtBin;
}

// Get the last loaded jet pT bin for jet-event plane histograms
int JetBackgroundHistogramManager::GetLastJetPtBin() const{
  return fLastLoadedJetPtBin;
}

// Getter for the number of events passing the cuts
int JetBackgroundHistogramManager::GetNEvents() const{
  return fhEvents->GetBinContent(fhEvents->FindBin(JetBackgroundHistograms::kVzCut));
}

// Getter for the JCard
JetBackgroundCard* JetBackgroundHistogramManager::GetCard() const{
  return fCard;
}
