// Class for the main analysis for jet background subtraction

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "JetBackgroundAnalyzer.h"

using namespace std;

/*
 * Default constructor
 */
JetBackgroundAnalyzer::JetBackgroundAnalyzer() :
  fFileNames(0),
  fCard(0),
  fHistograms(0),
  fVzWeightFunction(0),
  fCentralityWeightFunctionCentral(0),
  fCentralityWeightFunctionPeripheral(0),
  fSmearingFunction(0),
  fJetCorrector2018(),
  fCaloJetCorrector2018(),
  fRng(0),
  fJetType(0),
  fJetSubtraction(2),
  fDebugLevel(0),
  fSmearResolution(false),
  fDoCalorimeterJets(false),
  fVzWeight(1),
  fCentralityWeight(1),
  fPtHatWeight(1),
  fTotalEventWeight(1),
  fMaxParticleEtaEventPlane(2),
  fMaxParticlePtEventPlane(5),
  fJetAxis(0),
  fVzCut(0),
  fMinimumPtHat(0),
  fMaximumPtHat(0),
  fJetEtaCut(0),
  fJetMinimumPtCut(0),
  fJetMaximumPtCut(0),
  fMinimumMaxTrackPtFraction(0),
  fMaximumMaxTrackPtFraction(0),
  fJetClosureMinimumPt(0),
  fFillJetPtClosure(false)
{
  // Default constructor
  fHistograms = new JetBackgroundHistograms();
  fHistograms->CreateHistograms();
  
  // Initialize readers to null
  fEventReader = NULL;

  // Create a manager for jet energy resolution smearing in MC
  fEnergyResolutionSmearingFinder = new JetMetScalingFactorManager();
}

/*
 * Custom constructor
 */
JetBackgroundAnalyzer::JetBackgroundAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard) :
  fFileNames(fileNameVector),
  fCard(newCard),
  fHistograms(0),
  fJetCorrector2018(),
  fCaloJetCorrector2018(),
  fVzWeight(1),
  fCentralityWeight(1),
  fPtHatWeight(1),
  fTotalEventWeight(1)
{
  // Custom constructor
  fHistograms = new JetBackgroundHistograms(fCard);
  fHistograms->CreateHistograms();
  
  // Initialize readers to null
  fEventReader = NULL;
  
  // Configurure the analyzer from input card
  ReadConfigurationFromCard();
  
  // Function for smearing the jet pT for systemtic uncertainties
  fSmearingFunction = new TF1("fSmearingFunction","pol4",0,500);
    
  // The vz weight function is rederived from the miniAOD dataset.
  // Macro used for derivation: deriveMonteCarloWeights.C, Git hash: d4eab1cd188da72f5a81b8902cb6cc55ea1baf23
  // Input files: eecAnalysis_akFlowJet_onlyJets_weightEventInfo_combinedTriggers_processed_2023-03-06.root
  //              PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_wtaAxis_onlyJets_noTrigger_ptHatWeight_processed_2023-03-06.root
  fVzWeightFunction = new TF1("fvz","pol6",-15,15);
  fVzWeightFunction->SetParameters(1.00591, -0.0193751, 0.000961142, -2.44303e-05, -8.24443e-06, 1.66679e-07, 1.11028e-08);
    
  // The centrality weight function is rederived for the miniAOD dataset.
  // Macro used for derivation: deriveMonteCarloWeights.C, Git hash: d4eab1cd188da72f5a81b8902cb6cc55ea1baf23
  // Input files: eecAnalysis_akFlowJet_onlyJets_weightEventInfo_combinedTriggers_processed_2023-03-06.root
  //              PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_wtaAxis_onlyJets_noTrigger_ptHatWeight_processed_2023-03-06.root
  fCentralityWeightFunctionCentral = new TF1("fCentralWeight","pol6",0,30);
  fCentralityWeightFunctionCentral->SetParameters(4.73421, -0.0477343, -0.0332804, 0.00355699, -0.00017427, 4.18398e-06, -3.94746e-08);
  fCentralityWeightFunctionPeripheral = new TF1("fPeripheralWeight","pol6",30,90);
  fCentralityWeightFunctionPeripheral->SetParameters(3.38091, -0.0609601, -0.00228529, 9.43076e-05, -1.39593e-06, 9.85435e-09, -2.77153e-11);

  // Jet energy resolution smearing scale factor manager
  fEnergyResolutionSmearingFinder = new JetMetScalingFactorManager(true, JetMetScalingFactorManager::kNominal);

  // Initialize the random number generator with a random seed
  fRng = new TRandom3();
  fRng->SetSeed(0);
  
}

/*
 * Copy constructor
 */
JetBackgroundAnalyzer::JetBackgroundAnalyzer(const JetBackgroundAnalyzer& in) :
  fEventReader(in.fEventReader),
  fFileNames(in.fFileNames),
  fCard(in.fCard),
  fHistograms(in.fHistograms),
  fVzWeightFunction(in.fVzWeightFunction),
  fCentralityWeightFunctionCentral(in.fCentralityWeightFunctionCentral),
  fCentralityWeightFunctionPeripheral(in.fCentralityWeightFunctionPeripheral),
  fSmearingFunction(in.fSmearingFunction),
  fRng(in.fRng),
  fJetType(in.fJetType),
  fJetSubtraction(in.fJetSubtraction),
  fDebugLevel(in.fDebugLevel),
  fSmearResolution(in.fSmearResolution),
  fDoCalorimeterJets(in.fDoCalorimeterJets),
  fVzWeight(in.fVzWeight),
  fCentralityWeight(in.fCentralityWeight),
  fPtHatWeight(in.fPtHatWeight),
  fTotalEventWeight(in.fTotalEventWeight),
  fMaxParticleEtaEventPlane(in.fMaxParticleEtaEventPlane),
  fMaxParticlePtEventPlane(in.fMaxParticlePtEventPlane),
  fJetAxis(in.fJetAxis),
  fVzCut(in.fVzCut),
  fMinimumPtHat(in.fMinimumPtHat),
  fMaximumPtHat(in.fMaximumPtHat),
  fJetEtaCut(in.fJetEtaCut),
  fJetMinimumPtCut(in.fJetMinimumPtCut),
  fJetMaximumPtCut(in.fJetMaximumPtCut),
  fMinimumMaxTrackPtFraction(in.fMinimumMaxTrackPtFraction),
  fMaximumMaxTrackPtFraction(in.fMaximumMaxTrackPtFraction),
  fJetClosureMinimumPt(in.fJetClosureMinimumPt),
  fFillJetPtClosure(in.fFillJetPtClosure)
{
  // Copy constructor
}

/*
 * Assingment operator
 */
JetBackgroundAnalyzer& JetBackgroundAnalyzer::operator=(const JetBackgroundAnalyzer& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fEventReader = in.fEventReader;
  fFileNames = in.fFileNames;
  fCard = in.fCard;
  fHistograms = in.fHistograms;
  fVzWeightFunction = in.fVzWeightFunction;
  fCentralityWeightFunctionCentral = in.fCentralityWeightFunctionCentral;
  fCentralityWeightFunctionPeripheral = in.fCentralityWeightFunctionPeripheral;
  fSmearingFunction = in.fSmearingFunction;
  fRng = in.fRng;
  fJetType = in.fJetType;
  fJetSubtraction = in.fJetSubtraction;
  fDebugLevel = in.fDebugLevel;
  fSmearResolution = in.fSmearResolution;
  fDoCalorimeterJets = in.fDoCalorimeterJets;
  fVzWeight = in.fVzWeight;
  fCentralityWeight = in.fCentralityWeight;
  fPtHatWeight = in.fPtHatWeight;
  fTotalEventWeight = in.fTotalEventWeight;
  fMaxParticleEtaEventPlane = in.fMaxParticleEtaEventPlane;
  fMaxParticlePtEventPlane = in.fMaxParticlePtEventPlane;
  fJetAxis = in.fJetAxis;
  fVzCut = in.fVzCut;
  fMinimumPtHat = in.fMinimumPtHat;
  fMaximumPtHat = in.fMaximumPtHat;
  fJetEtaCut = in.fJetEtaCut;
  fJetMinimumPtCut = in.fJetMinimumPtCut;
  fJetMaximumPtCut = in.fJetMaximumPtCut;
  fMinimumMaxTrackPtFraction = in.fMinimumMaxTrackPtFraction;
  fMaximumMaxTrackPtFraction = in.fMaximumMaxTrackPtFraction;
  fJetClosureMinimumPt = in.fJetClosureMinimumPt;
  fFillJetPtClosure = in.fFillJetPtClosure;
  
  return *this;
}

/*
 * Destructor
 */
JetBackgroundAnalyzer::~JetBackgroundAnalyzer(){
  // destructor
  delete fHistograms;
  if(fVzWeightFunction) delete fVzWeightFunction;
  if(fJetCorrector2018) delete fJetCorrector2018;
  if(fCaloJetCorrector2018) delete fCaloJetCorrector2018;
  if(fEnergyResolutionSmearingFinder) delete fEnergyResolutionSmearingFinder;
  if(fCentralityWeightFunctionCentral) delete fCentralityWeightFunctionCentral;
  if(fCentralityWeightFunctionPeripheral) delete fCentralityWeightFunctionPeripheral;
  if(fSmearingFunction) delete fSmearingFunction;
  if(fRng) delete fRng;
  if(fEventReader) delete fEventReader;
}

/*
 * Read all the configuration from the input card
 */
void JetBackgroundAnalyzer::ReadConfigurationFromCard(){
  
  //****************************************
  //         Event selection cuts
  //****************************************
  
  fVzCut = fCard->Get("ZVertexCut");          // Event cut vor the z-position of the primary vertex
  fMinimumPtHat = fCard->Get("LowPtHatCut");  // Minimum accepted pT hat value
  fMaximumPtHat = fCard->Get("HighPtHatCut"); // Maximum accepted pT hat value

  //****************************************
  //      Event plane calculation cuts
  //****************************************

  fMaxParticleEtaEventPlane = fCard->Get("MaxParticleEtaEventPlane"); // Maximum eta value for particles used to determine the event plane
  fMaxParticlePtEventPlane = fCard->Get("MaxParticlePtEventPlane");   // Maximum pT value for particles used to determine the event plane
  
  //****************************************
  //          Jet selection cuts
  //****************************************
  
  fJetEtaCut = fCard->Get("JetEtaCut");           // Eta cut around midrapidity
  fJetMinimumPtCut = fCard->Get("MinJetPtCut");   // Minimum pT cut for jets
  fJetMaximumPtCut = fCard->Get("MaxJetPtCut");   // Maximum pT accepted for jets (and tracks)
  fMinimumMaxTrackPtFraction = fCard->Get("MinMaxTrackPtFraction");  // Cut for jets consisting only from soft particles
  fMaximumMaxTrackPtFraction = fCard->Get("MaxMaxTrackPtFraction");  // Cut for jets consisting only from one high pT particle
  fJetClosureMinimumPt = fCard->Get("MinJetPtClosure"); // Minimum jet pT for closure histograms

  //****************************************
  //            Jet selection
  //****************************************
  fJetType = fCard->Get("JetType");               // Select the type of analyzed jets (Reconstructed / Generator level)
  fJetSubtraction = fCard->Get("JetSubtraction"); // Select the background subtracted algorithm (Calo PU, PFCS, flow PFCS)
  fJetAxis = fCard->Get("JetAxis");               // Select between E-escheme and WTA axes
  fSmearResolution = (fCard->Get("SmearResolution") == 1); // Flag for smearing the jet resolution in MC
  fDoCalorimeterJets = (fCard->Get("DoCaloJets") == 1);    // Flag for filling calorimeter jet histograms

  //***************************************
  //            Jet pT closure
  //***************************************
  fFillJetPtClosure = (fCard->Get("FillJetPtClosure") == 1); // Flag to fill jet pT closure histograms
  
  //************************************************
  //              Debug messages
  //************************************************
  fDebugLevel = fCard->Get("DebugLevel");
}

/*
 * Main analysis loop
 */
void JetBackgroundAnalyzer::RunAnalysis(){
  
  //************************************************
  //  Define variables needed in the analysis loop
  //************************************************
  
  // Input files and forest readers for analysis
  TFile* inputFile;
  
  // Event variables
  Int_t nEvents = 0;                // Number of events
  Double_t vz = 0;                  // Vertex z-position
  Double_t centrality = 0;          // Event centrality
  Int_t hiBin = 0;                  // CMS hiBin (centrality * 2)
  Double_t ptHat = 0;               // pT hat for MC events
  
  // Variables for jets
  Int_t nJets = 0;                  // Number of jets in an event
  Double_t jetPt = 0;               // pT of the i:th jet in the event
  Double_t jetPhi = 0;              // phi of the i:th jet in the event
  Double_t jetEta = 0;              // eta of the i:th jet in the event
  Int_t jetFlavor = 0;              // Flavor of the jet. 0 = Quark jet. 1 = Gluon jet.
  Int_t matchingJetExists = 0;      // Flag for having a matching jet. 0 = No match found. 1 = Match exist

  // Variables for leading jet
  Double_t leadingJetPt = 0;        // pT of the leading jet
  Double_t leadingJetPhi = 0;       // phi of the leading jet
  Double_t leadingJetEta = 0;       // eta of the leading jet
  Int_t leadingJetFlavor = 0;       // Flavor of the leading jet. 0 = Quark jet. 1 = Gluon jet, 2 = Flavor undetermined
  Int_t leadingJetMatch = 0;        // 0 = No match found. 1 = Match exist

  // Variables for matched reconstructed jet
  Double_t reconstructedJetPt = 0;   // pT of the reconstructed jet
  Double_t reconstructedJetPhi = 0;  // phi of the reconstructed jet
  Double_t reconstructedJetEta = 0;  // eta of the reconstructed jet

  // Variables for particles
  Int_t nParticles = 0;             // Number of generator level particles
  Double_t particlePt = 0;          // pT of a generator level particle
  Double_t particleEta = 0;         // eta of a generator level particle
  Double_t particlePhi = 0;         // phi of a generator level particle
  
  // Variables for smearing study
  Double_t smearingFactor = 0;       // Larger of the JEC uncertainties
  
  // Variables for jet matching and closure
  Int_t partonFlavor = -999;        // Code for parton flavor in Monte Carlo

  // Event plane study related variables
  const Int_t nFlowComponentsEP = 3;                  // Number of flow component to which the event plane is determined
  Double_t eventPlaneQ[nFlowComponentsEP] = {0};      // Magnitude of the event plane Q-vector
  Double_t eventPlaneMultiplicity = 0;                // Particle multiplicity in the event plane
  Double_t eventPlaneQx[nFlowComponentsEP] = {0};     // x-component of the event plane vector
  Double_t eventPlaneQy[nFlowComponentsEP] = {0};     // y-component of the event plane vector
  Double_t jetEventPlaneDeltaPhi = 0;                 // DeltaPhi between jet and event plane angle
  Double_t eventPlaneAngle[nFlowComponentsEP] = {0};  // Manually calculated event plane angle
  
  // File name helper variables
  TString currentFile;
  
  // Fillers for THnSparses
  const Int_t nFillJet = 6;         // Inclusive and leading jets
  const Int_t nFillEventPlane = 3;  // Correlation between inclusive and leading jets with event plane
  const Int_t nAxesClosure = 7;     // Jet pT closure
  Double_t fillerJet[nFillJet];
  Double_t fillerEventPlane[nFillEventPlane];
  Double_t fillerClosure[nAxesClosure];
  
  // For 2018 PbPb and 2017 pp data, we need to correct jet pT
  std::string correctionFileRelative = "jetEnergyCorrections/Autumn18_HI_V8_MC_L2Relative_AK4PF.txt";
  std::string correctionFileCalo = "jetEnergyCorrections/Autumn18_HI_V8_MC_L2Relative_AK4Calo.txt";
  
  vector<string> correctionFiles;
  correctionFiles.push_back(correctionFileRelative);
  fJetCorrector2018 = new JetCorrector(correctionFiles);

  vector<string> correctionFilesCalo;
  correctionFilesCalo.push_back(correctionFileCalo);
  fCaloJetCorrector2018 = new JetCorrector(correctionFilesCalo);
  
  //************************************************
  //      Find forest readers for data files
  //************************************************

  fEventReader = new MonteCarloForestReader(fJetSubtraction, fJetAxis);
  
  //************************************************
  //       Main analysis loop over all files
  //************************************************
  
  // Loop over files
  Int_t nFiles = fFileNames.size();
  for(Int_t iFile = 0; iFile < nFiles; iFile++) {
    
    //************************************************
    //              Find and open files
    //************************************************
    
    // Find the filename and open the input file
    currentFile = fFileNames.at(iFile);
    inputFile = TFile::Open(currentFile);
    
    // Check that the file exists
    if(!inputFile){
      cout << "Error! Could not find the file: " << currentFile.Data() << endl;
      assert(0);
    }

    // Check that the file is open
    if(!inputFile->IsOpen()){
      cout << "Error! Could not open the file: " << currentFile.Data() << endl;
      assert(0);
    }
    
    // Check that the file is not zombie
    if(inputFile->IsZombie()){
      cout << "Error! The following file is a zombie: " << currentFile.Data() << endl;
      assert(0);
    }
    
    // Print the used files
    if(fDebugLevel > 0) cout << "Reading from file: " << currentFile.Data() << endl;
    
    //************************************************
    //            Read forest from file
    //************************************************
    
    // If file is good, read the forest from the file
    fEventReader->ReadForestFromFile(inputFile);  // There might be a memory leak in handling the forest...

    nEvents = fEventReader->GetNEvents();

    //************************************************
    //         Main event loop for each file
    //************************************************
    
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){ // nEvents

      // For each event, chack that the file stays open:
      // This is to try to combat file read errors occasionally happening during CRAB running.
      // Will need to monitor the situation and see if this really works.
      if(!inputFile->IsOpen() || inputFile->IsZombie()){
        cout << "Error! Lost access to the file: " << currentFile.Data() << endl;
        assert(0);
      }
      
      //************************************************
      //         Read basic event information
      //************************************************
      
      // Print to console how the analysis is progressing
      if(fDebugLevel > 1 && iEvent % 1000 == 0) cout << "Analyzing event " << iEvent << endl;
      
      // Read the event to memory
      fEventReader->GetEvent(iEvent);

      // Get vz, centrality and pT hat information
      vz = fEventReader->GetVz();
      centrality = fEventReader->GetCentrality();
      hiBin = fEventReader->GetHiBin();
      ptHat = fEventReader->GetPtHat();
      
      // We need to apply pT hat cuts before getting pT hat weight. There might be rare events above the upper
      // limit from which the weights are calculated, which could cause the code to crash.
      if(ptHat < fMinimumPtHat || ptHat >= fMaximumPtHat) continue;
      
      // Get the weighting for the event
      fVzWeight = GetVzWeight(vz);  // vz weight
      fCentralityWeight = GetCentralityWeight(hiBin); // centrality weight

      // Event weight for 2018 MC
      fPtHatWeight = fEventReader->GetEventWeight(); // 2018 MC
      fTotalEventWeight = fVzWeight*fCentralityWeight*fPtHatWeight;
      
      // Fill event counter histogram
      fHistograms->fhEvents->Fill(JetBackgroundHistograms::kAll);          // All the events looped over
      
      //  ============================================
      //  ===== Apply all the event quality cuts =====
      //  ============================================

      // Check event cuts
      if(!PassEventCuts(fEventReader, true)) continue;
      
      // Fill the event information histograms for the events that pass the event cuts
      fHistograms->fhVertexZ->Fill(vz,fPtHatWeight);                         // z vertex distribution from all events
      fHistograms->fhVertexZWeighted->Fill(vz,fTotalEventWeight);            // z-vertex distribution weighted with the weight function
      fHistograms->fhCentrality->Fill(centrality, fPtHatWeight);             // Centrality filled from all events
      fHistograms->fhCentralityWeighted->Fill(centrality,fTotalEventWeight); // Centrality weighted with the centrality weighting function
      fHistograms->fhPtHat->Fill(ptHat);                                     // pT hat histogram
      fHistograms->fhPtHatWeighted->Fill(ptHat,fTotalEventWeight);           // pT het histogram weighted with corresponding cross section and event number
      
      // ======================================
      // ===== Event quality cuts applied =====
      // ======================================

      //******************************************************************
      //    Determine the event plane from generator level information
      //******************************************************************

      // Before calculating the event plane, reset all the variables from the previous event
      eventPlaneMultiplicity = 0;
      for(Int_t iFlow = 0; iFlow < nFlowComponentsEP; iFlow++){
        eventPlaneQ[iFlow] = 0;
        eventPlaneQx[iFlow] = 0;
        eventPlaneQy[iFlow] = 0;
      }

      // Loop over all generator level particles in the event
      nParticles = fEventReader->GetNGenParticles();
      for(Int_t iParticle = 0; iParticle < nParticles; iParticle++){

        // Get the particle information
        particlePt = fEventReader->GetGenParticlePt(iParticle);
        particleEta = fEventReader->GetGenParticleEta(iParticle);
        particlePhi = fEventReader->GetGenParticlePhi(iParticle);

        // Cuts for particles used in event plane calculation
        if(TMath::Abs(particleEta) > fMaxParticleEtaEventPlane) continue;  // Only consider particles from mid-rapidity
        if(fEventReader->GetGenParticleSubevent(iParticle) == 0) continue; // Only use Hydjet-particles for event plane calculation
        if(particlePt > fMaxParticlePtEventPlane) continue;  // Ignore high-pT particles for event plane calculation

        // Determine the event planes from order 2 to order 2+nFlowComponentsEP-1
        for(int iFlow = 0; iFlow < nFlowComponentsEP; iFlow++){
          eventPlaneQx[iFlow] += TMath::Cos((iFlow+2.0)*(particlePhi));
          eventPlaneQy[iFlow] += TMath::Sin((iFlow+2.0)*(particlePhi));
        }
        eventPlaneMultiplicity += 1;

      }

      // Do not allow zero multiplicity to avoid dividing by zero problems
      if(eventPlaneMultiplicity == 0) eventPlaneMultiplicity += 1;
          
      // Calculate the Q-vector magnitudes and event plane angles for orders 2 tp 2+nFlowComponentsEP-1
      for(int iFlow = 0; iFlow < nFlowComponentsEP; iFlow++){
        eventPlaneQ[iFlow] = TMath::Sqrt(eventPlaneQx[iFlow]*eventPlaneQx[iFlow] + eventPlaneQy[iFlow]*eventPlaneQy[iFlow]);
        eventPlaneAngle[iFlow] = (1.0/(iFlow+2.0)) * TMath::ATan2(eventPlaneQy[iFlow], eventPlaneQx[iFlow]);
      }

      // Normalize the Q-vector with multiplicity
      for(int iFlow = 0; iFlow < nFlowComponentsEP; iFlow++){
        eventPlaneQ[iFlow] /= TMath::Sqrt(eventPlaneMultiplicity);
      }

      //***********************************************************
      //       First jet loop for event plane correlations
      //***********************************************************

      // Reser the leading jet variables for this event
      leadingJetPt = 0;
      leadingJetEta = -999;
      leadingJetPhi = -999;
      leadingJetFlavor = -999;
      leadingJetMatch = -999;

      // Jet loop
      nJets = fEventReader->GetNJets(fJetType);
      for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++){
        
        jetPt = fEventReader->GetJetRawPt(fJetType, jetIndex);  // Get the raw pT and do manual correction later
        jetPhi = fEventReader->GetJetPhi(fJetType, jetIndex);
        jetEta = fEventReader->GetJetEta(fJetType, jetIndex);
        jetFlavor = -999;
          
        //  ========================================
        //  ======== Apply jet quality cuts ========
        //  ========================================
          
        if(TMath::Abs(jetEta) >= fJetEtaCut) continue; // Cut for jet eta
          
        // No jet quality cuts for generator level jets
        if(!(fJetType == MonteCarloForestReader::kGeneratorLevelJet)){              
          if(fMinimumMaxTrackPtFraction >= fEventReader->GetJetMaxTrackPt(jetIndex)/fEventReader->GetJetRawPt(jetIndex)){
            continue; // Cut for jets with only very low pT particles
          }
          if(fMaximumMaxTrackPtFraction <= fEventReader->GetJetMaxTrackPt(jetIndex)/fEventReader->GetJetRawPt(jetIndex)){
            continue; // Cut for jets where all the pT is taken by one track
          }
        }
          
        //  ========================================
        //  ======= Jet quality cuts applied =======
        //  ========================================

        // No jet pT correction or smearing for generator level jets
        if(!(fJetType == MonteCarloForestReader::kGeneratorLevelJet)){
    
          // For reconstructed jets do a correction for the jet pT
          fJetCorrector2018->SetJetPT(jetPt);
          fJetCorrector2018->SetJetEta(jetEta);
          fJetCorrector2018->SetJetPhi(jetPhi);
          
          jetPt = fJetCorrector2018->GetCorrectedPT();
          
          // Apply gaussian smearing to take into account overly optimistic jet energy resolution
          if(fSmearResolution){
            smearingFactor = GetSmearingFactor(jetPt, jetEta, centrality);
            jetPt = jetPt * fRng->Gaus(1,smearingFactor);
          }
            
        } // Jet pT correction
          
        // After the jet pT can been corrected, apply analysis jet pT cuts
        if(jetPt < fJetMinimumPtCut) continue;
        if(jetPt > fJetMaximumPtCut) continue;

        // Check if the current jet has a matching jet
        matchingJetExists = 0;
        if(fEventReader->HasMatchingJet(fJetType, jetIndex)){
          // Require that one pT is not less than half of the other pT
          if(jetPt*0.5 < fEventReader->GetMatchedPt(fJetType, jetIndex) && fEventReader->GetMatchedPt(fJetType, jetIndex) * 0.5 < jetPt){
            matchingJetExists = 1;
          }
        }

        // Find the jet flavor and translate it into a quark [-6,-1] U [1,6] or gluon (21)
        // In the jet flavor is not any of these values, it remains undeterined
        jetFlavor = JetBackgroundHistograms::kUndetermined;
        partonFlavor = fEventReader->GetJetFlavor(fJetType, jetIndex);
        if(TMath::Abs(partonFlavor) == 21) jetFlavor = JetBackgroundHistograms::kGluon;
        if(TMath::Abs(partonFlavor) < 7){
          if(partonFlavor != 0) jetFlavor = JetBackgroundHistograms::kQuark;
        }

        // After the event selection, update the leading jet variables
        if(jetPt > leadingJetPt){
          leadingJetPt = jetPt;
          leadingJetPhi = jetPhi;
          leadingJetEta = jetEta;
          leadingJetFlavor = jetFlavor;
          leadingJetMatch = matchingJetExists;
        }
        
        //************************************************
        //         Fill histograms for all jets
        //************************************************
          
        // Fill the axes in correct order
        fillerJet[0] = jetPt;             // Axis 0 = any jet pT
        fillerJet[1] = jetPhi;            // Axis 1 = any jet phi
        fillerJet[2] = jetEta;            // Axis 2 = any jet eta
        fillerJet[3] = centrality;        // Axis 3 = centrality
        fillerJet[4] = jetFlavor;         // Axis 4 = flavor of the jet
        fillerJet[5] = matchingJetExists; // Axis 5 = flag is matching jet exists
          
        fHistograms->fhInclusiveJet->Fill(fillerJet,fTotalEventWeight); // Fill the data point to histogram

        //**********************************************************************
        //      Fill histograms for inclusive jet - event plane correlation
        //**********************************************************************

        for(int iFlow = 0; iFlow < nFlowComponentsEP; iFlow++){
          
          // Determine the deltaPhi between jet axis and the event plane in the interval [-pi/2,3pi/2]
          jetEventPlaneDeltaPhi = jetPhi - eventPlaneAngle[iFlow];
          while(jetEventPlaneDeltaPhi > (1.5*TMath::Pi())){jetEventPlaneDeltaPhi += -2*TMath::Pi();}
          while(jetEventPlaneDeltaPhi < (-0.5*TMath::Pi())){jetEventPlaneDeltaPhi += 2*TMath::Pi();}

          // Require matching generator level jet with at least 80 GeV
          if(fEventReader->HasMatchingGenJet(jetIndex)){

            // Fill the jet - event plane correlation histograms
            fillerEventPlane[0] = jetEventPlaneDeltaPhi;  // Axis 0: DeltaPhi between jet and event plane
            fillerEventPlane[1] = jetPt;                  // Axis 1: Jet pT
            fillerEventPlane[2] = centrality;             // Axis 2: centrality

            fHistograms->fhInclusiveJetEventPlane[iFlow]->Fill(fillerEventPlane, fTotalEventWeight);
          }

        }
        
      } // End of jet loop

      // If a leading jet exists, fill the leading jet histograms
      if(leadingJetPt > 0){

        //***************************************************
        //         Fill histograms for leading jets
        //***************************************************

        // Fill the axes in correct order
        fillerJet[0] = leadingJetPt;       // Axis 0 = leading jet pT
        fillerJet[1] = leadingJetPhi;      // Axis 1 = leading jet phi
        fillerJet[2] = leadingJetEta;      // Axis 2 = leading jet eta
        fillerJet[3] = centrality;         // Axis 3 = centrality
        fillerJet[4] = leadingJetFlavor;   // Axis 4 = flavor of the leading jet
        fillerJet[5] = leadingJetMatch;    // Axis 5 = flag if matching jet exists
          
        fHistograms->fhLeadingJet->Fill(fillerJet,fTotalEventWeight); // Fill the data point to histogram

        //**********************************************************************
        //      Fill histograms for leading jet - event plane correlation
        //**********************************************************************

        for(int iFlow = 0; iFlow < nFlowComponentsEP; iFlow++){
          
          // Determine the deltaPhi between jet axis and the event plane in the interval [-pi/2,3pi/2]
          jetEventPlaneDeltaPhi = leadingJetPhi - eventPlaneAngle[iFlow];
          while(jetEventPlaneDeltaPhi > (1.5*TMath::Pi())){jetEventPlaneDeltaPhi += -2*TMath::Pi();}
          while(jetEventPlaneDeltaPhi < (-0.5*TMath::Pi())){jetEventPlaneDeltaPhi += 2*TMath::Pi();}

          // Fill the jet - event plane correlation histograms
          fillerEventPlane[0] = jetEventPlaneDeltaPhi;  // Axis 0: DeltaPhi between jet and event plane
          fillerEventPlane[1] = leadingJetPt;           // Axis 1: Leading jet pT
          fillerEventPlane[2] = centrality;             // Axis 2: centrality

          fHistograms->fhLeadingJetEventPlane[iFlow]->Fill(fillerEventPlane, fTotalEventWeight);

        }
      } // Filling leading jet histograms

      //*******************************************************************
      //     If selected, fill the histograms also for calorimeter jets
      //*******************************************************************
      if(fDoCalorimeterJets){

        nJets = fEventReader->GetNJets(fJetType);
        for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++){

          // Find the calorimeter jet kinematics
          jetPt = fEventReader->GetCalorimeterJetPt(jetIndex);
          jetPhi = fEventReader->GetCalorimeterJetPhi(jetIndex);
          jetEta = fEventReader->GetCalorimeterJetEta(jetIndex);

          // Select the jets from a defined eta region
          if(TMath::Abs(jetEta) >= fJetEtaCut) continue; // Cut for jet eta

          // Do jet energy correction for calorimeter jets
          fCaloJetCorrector2018->SetJetPT(jetPt);
          fCaloJetCorrector2018->SetJetEta(jetEta);
          fCaloJetCorrector2018->SetJetPhi(jetPhi);
          
          jetPt = fCaloJetCorrector2018->GetCorrectedPT();

          // After the jet pT can been corrected, apply analysis jet pT cuts
          if(jetPt < fJetMinimumPtCut) continue;
          if(jetPt > fJetMaximumPtCut) continue;

          //************************************************
          //         Fill histograms for all jets
          //************************************************
          
          // Fill the axes in correct order
          fillerJet[0] = jetPt;             // Axis 0 = calorimeter jet pT
          fillerJet[1] = jetPhi;            // Axis 1 = calorimeter jet phi
          fillerJet[2] = jetEta;            // Axis 2 = calorimeter jet eta
          fillerJet[3] = centrality;        // Axis 3 = centrality
          fillerJet[4] = 0;                 // Axis 4 = not used for calorimeter jets
          fillerJet[5] = 0;                 // Axis 5 = not used for calorimeter jets
          
          fHistograms->fhCalorimeterJet->Fill(fillerJet,fTotalEventWeight); // Fill the data point to histogram

          //**********************************************************************
          //      Fill histograms for calorimeter jet - event plane correlation
          //**********************************************************************

          for(int iFlow = 0; iFlow < nFlowComponentsEP; iFlow++){
          
            // Determine the deltaPhi between jet axis and the event plane in the interval [-pi/2,3pi/2]
            jetEventPlaneDeltaPhi = jetPhi - eventPlaneAngle[iFlow];
            while(jetEventPlaneDeltaPhi > (1.5*TMath::Pi())){jetEventPlaneDeltaPhi += -2*TMath::Pi();}
            while(jetEventPlaneDeltaPhi < (-0.5*TMath::Pi())){jetEventPlaneDeltaPhi += 2*TMath::Pi();}

            // Fill the jet - event plane correlation histograms
            fillerEventPlane[0] = jetEventPlaneDeltaPhi;  // Axis 0: DeltaPhi between jet and event plane
            fillerEventPlane[1] = jetPt;                  // Axis 1: Jet pT
            fillerEventPlane[2] = centrality;             // Axis 2: centrality

            fHistograms->fhCalorimeterJetEventPlane[iFlow]->Fill(fillerEventPlane, fTotalEventWeight);
          } // Flow order loop

        } // Calorimeter jet loop
      } // Calorimeter jet if

      //**************************************************
      //       Second jet loop for jet pT closure
      //**************************************************

      // Only fill the jet pT closure plots if selected
      if(!fFillJetPtClosure) continue;

      // Loop over all generator level jets
      nJets = fEventReader->GetNGeneratorJets();
      for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++){

        jetPt = fEventReader->GetGeneratorJetPt(jetIndex);
        jetPhi = fEventReader->GetGeneratorJetPhi(jetIndex);
        jetEta = fEventReader->GetGeneratorJetEta(jetIndex);

        // Kinematic cuts for generator level jets
        if(TMath::Abs(jetEta) >= fJetEtaCut) continue; // Cut for jet eta
        if(jetPt < fJetClosureMinimumPt) continue;     // Cut for jet pT
        if(jetPt > fJetMaximumPtCut) continue;         // Cut for super high pT jets

        // For closure plots, we need to find a matching reconstructed jet
        if(!fEventReader->HasMatchingRecoJet(jetIndex)) continue;

        // Read the reconstructed jet information
        reconstructedJetPt = fEventReader->GetMatchedRecoPt(jetIndex);
        reconstructedJetEta = fEventReader->GetMatchedRecoEta(jetIndex);
        reconstructedJetPhi = fEventReader->GetMatchedRecoPhi(jetIndex);
        partonFlavor = fEventReader->GetGenJetFlavor(jetIndex);

        // Apply jet energy correction for reconstructed jet
        fJetCorrector2018->SetJetPT(reconstructedJetPt);
        fJetCorrector2018->SetJetEta(reconstructedJetEta);
        fJetCorrector2018->SetJetPhi(reconstructedJetPhi);

        reconstructedJetPt = fJetCorrector2018->GetCorrectedPT();
          
        // Apply gaussian smearing to take into account too good jet energy resolution
        if(fSmearResolution){
          smearingFactor = GetSmearingFactor(reconstructedJetPt, reconstructedJetEta, centrality);
          reconstructedJetPt = reconstructedJetPt * fRng->Gaus(1,smearingFactor);
        }

        // Define index for jet flavor using algoritm: [-6,-1] U [1,6] -> kQuark, 21 -> kGluon, anything else -> kUndetermined
        jetFlavor = JetBackgroundHistograms::kUndetermined;
        if(partonFlavor >= -6 && partonFlavor <= 6 && partonFlavor != 0) jetFlavor = JetBackgroundHistograms::kQuark;
        if(partonFlavor == 21) jetFlavor = JetBackgroundHistograms::kGluon;

        //************************************************
        //       Fill histograms for jet pT closure
        //************************************************

        // Fill the different axes for the filler
        fillerClosure[0] = jetPt;                    // Axis 0: pT of the matched generator level jet
        fillerClosure[1] = reconstructedJetPt;       // Axis 1: pT of the matched reconstructed jet
        fillerClosure[2] = jetEta;                   // Axis 2: eta of the jet under consideration
        fillerClosure[3] = centrality;               // Axis 3: Centrality of the event
        fillerClosure[4] = jetFlavor;                // Axis 4: Jet flavor type (quark/gluon)
        fillerClosure[5] = reconstructedJetPt/jetPt; // Axis 5: Reconstructed level jet to generator level jet pT ratio
        fillerClosure[6] = jetPhi;                   // Axis 6: phi of the jet under consideration
  
        // Fill the closure histogram
        fHistograms->fhJetPtClosure->Fill(fillerClosure,fTotalEventWeight);

      } // Jet pT loop for closures
      
    } // Event loop
    
    //************************************************
    //      Cleanup at the end of the file loop
    //************************************************
    
    // Close the input files after the event has been read
    inputFile->Close();
    
  } // File loop
  
}

/*
 * Get the proper vz weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t vz = Vertex z position for the event
 *
 *   return: Multiplicative correction factor for vz
 */
Double_t JetBackgroundAnalyzer::GetVzWeight(const Double_t vz) const{
  return fVzWeightFunction->Eval(vz); // Weight for 2018 MC
}

/*
 * Get the proper centrality weighting depending on analyzed system
 *
 *  Arguments:
 *   const Int_t hiBin = CMS hiBin
 *
 *   return: Multiplicative correction factor for the given CMS hiBin
 */
Double_t JetBackgroundAnalyzer::GetCentralityWeight(const Int_t hiBin) const{

  // No weighting for the most peripheral centrality bins. Different weight function for central and peripheral.
  if(hiBin < 60) return fCentralityWeightFunctionCentral->Eval(hiBin/2.0);
  return (hiBin < 194) ? fCentralityWeightFunctionPeripheral->Eval(hiBin/2.0) : 1;
}

/*
 * Get a smearing factor corresponding to worsening the smearing resolution in MC by 20 %
 * This is obtained by multiplying the MC smearing resolution by 0.666 and using this as additional
 * smearing for the data. Smearing factor depends on jet pT and centrality.
 *
 *  Arguments:
 *   Double_t jetPt = Jet pT
 *   const Double_t centrality = Centrality of the event
 *
 *  return: Additional smearing factor
 */
Double_t JetBackgroundAnalyzer::GetSmearingFactor(Double_t jetPt, Double_t jetEta, const Double_t centrality) {
  
  // For all the jets above 500 GeV, use the resolution for 500 GeV jet
  if(jetPt > 500) jetPt = 500;
  
  // Find the correct centrality bin
  Int_t centralityBin = GetCentralityBin(centrality);
    
  // Parameters for the smearing function
  // Determined using the macro constructJetPtClosures.C
  // Input file: PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_jetPtClosure_finalMcWeight_processed_2023-03-06.root
  Double_t resolutionFit[4][5] = {
    {0.424589, -0.00260826, 8.06713e-06, -1.17528e-08, 6.49404e-12},
    {0.415669, -0.00300066, 1.18443e-05, -2.30594e-08, 1.74213e-11},
    {0.358008, -0.00270528, 1.10689e-05, -2.15565e-08, 1.59491e-11},
    {0.237325, -0.00138461, 4.77259e-06, -7.80495e-09, 4.79538e-12}
  };
    
  for(int iParameter = 0; iParameter < 5; iParameter++){
    // Settings for PbPb
      
    fSmearingFunction->SetParameter(iParameter, resolutionFit[centralityBin][iParameter]);
  }
  
  // Calculation for resolution worsening: we assume the jet energy resolution is a Gaussian distribution with some certain sigma, if you would like to add a Gaussian noise to make it worse, the sigma getting larger, then it obeys the random variable rule that X=Y+Z, where Y~N(y, sigmay) and Z~N(z,sigmaz), then X~N(y+z, sqrt(sigmay^2+sigmaz^2))). In this case, we assume that noise and the resolution are independent.
  // So let assume the sigmay is the jet energy resolution, then you want the sigmax = 1.2sigmay
  // which means that the sigmaz = sigmay * sqrt(1.2^2-1)
  
  // After the smearing function is set, read the value to return
  // Worsening resolution by 20%: 0.663
  // Worsening resolution by 10%: 0.458
  // Worsening resolution by 30%: 0.831

  // We want to worsen resolution in MC by the amount defined by JetMet group. The scaling factor is given by a JetMet manager
  return fSmearingFunction->Eval(jetPt)*fEnergyResolutionSmearingFinder->GetScalingFactor(jetEta);
  
}

/*
 * Check if the event passes all the track cuts
 *
 *  Arguments:
 *   MonteCarloForestReader* eventReader = MonteCarloForestReader containing the event information checked for event cuts
 *   const Bool_t fillHistograms = Flag for filling the event information histograms.
 *
 *   return = True if all event cuts are passed, false otherwise
 */
Bool_t JetBackgroundAnalyzer::PassEventCuts(MonteCarloForestReader* eventReader, const Bool_t fillHistograms){
  
  // Primary vertex has at least two tracks, is within 25 cm in z-rirection and within 2 cm in xy-direction. Only applied for data.
  if(eventReader->GetPrimaryVertexFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(JetBackgroundHistograms::kPrimaryVertex);
  
  // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV. Only applied for PbPb data.
  if(eventReader->GetHfCoincidenceFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(JetBackgroundHistograms::kHfCoincidence);
  
  // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible. Only applied for PbPb data.
  if(eventReader->GetClusterCompatibilityFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(JetBackgroundHistograms::kClusterCompatibility);
  
  // Cut for vertex z-position
  if(TMath::Abs(eventReader->GetVz()) > fVzCut) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(JetBackgroundHistograms::kVzCut);
  
  return true;
  
}

/*
 * Getter for EEC histograms
 */
JetBackgroundHistograms* JetBackgroundAnalyzer::GetHistograms() const{
  return fHistograms;
}

/*
 * Getter for centrality bin
 */
Int_t JetBackgroundAnalyzer::GetCentralityBin(const Double_t centrality) const{
  
  // Find the correct centrality bin
  Int_t centralityBin = 0;
  for(int iCentrality = 1; iCentrality < fCard->GetNBin("CentralityBinEdges"); iCentrality++){
    if(centrality > fCard->Get("CentralityBinEdges",iCentrality)) centralityBin++;
  }

  return centralityBin;
}

/*
 * Get deltaR between two objects
 *
 *  Arguments:
 *   const Double_t eta1 = Eta of the first object
 *   const Double_t phi1 = Phi of the first object
 *   const Double_t eta2 = Eta of the second object
 *   const Double_t phi2 = Phi of the second object
 *
 *  return: DeltaR between the two objects
 */
Double_t JetBackgroundAnalyzer::GetDeltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2) const{
  
  Double_t deltaEta = eta1 - eta2;
  Double_t deltaPhi = phi1 - phi2;
  
  // Transform deltaPhi to interval [-pi,pi]
  while(deltaPhi > TMath::Pi()){deltaPhi += -2*TMath::Pi();}
  while(deltaPhi < -TMath::Pi()){deltaPhi += 2*TMath::Pi();}
  
  // Return the distance between the objects
  return TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
  
}