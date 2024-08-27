// Implementation for MonteCarloForestReader

// Own includes
#include "MonteCarloForestReader.h"

/*
 * Default constructor
 */
MonteCarloForestReader::MonteCarloForestReader() :
  fJetType(0),
  fJetAxis(0),
  fHeavyIonTree(0),
  fSkimTree(0),
  fJetTree(0),
  fTrackTree(0),
  fGenParticleTree(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fPtHatBranch(0),
  fEventWeightBranch(0),
  fPrimaryVertexBranch(0),
  fHfCoincidenceBranch(0),
  fClusterCompatibilityBranch(0),
  fnJetsBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetWTAPhiBranch(0),
  fJetEtaBranch(0),
  fJetWTAEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fJetRefPtBranch(0),
  fJetRefEtaBranch(0),
  fJetRefPhiBranch(0),
  fJetRefFlavorBranch(0),
  fnGenJetsBranch(0),
  fGenJetPtBranch(0),
  fGenJetPhiBranch(0),
  fGenJetWTAPhiBranch(0),
  fGenJetEtaBranch(0),
  fGenJetWTAEtaBranch(0),
  fnCaloJetsBranch(0),
  fCaloJetPtBranch(0),
  fCaloJetPhiBranch(0),
  fCaloJetEtaBranch(0),
  fnTracksBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fHighPurityTrackBranch(0),
  fTrackVertexDistanceZBranch(0),
  fTrackVertexDistanceZErrorBranch(0),
  fTrackVertexDistanceXYBranch(0),
  fTrackVertexDistanceXYErrorBranch(0),
  fTrackChi2Branch(0),
  fnTrackDegreesOfFreedomBranch(0),
  fnHitsTrackerLayerBranch(0),
  fnHitsTrackBranch(0),
  fTrackEnergyEcalBranch(0),
  fTrackEnergyHcalBranch(0),
  fTrackChargeBranch(0),
  fGenParticlePtBranch(0),
  fGenParticlePhiBranch(0),
  fGenParticleEtaBranch(0),
  fGenParticleChargeBranch(0),
  fGenParticleSubeventBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fPtHat(0),
  fPrimaryVertexFilterBit(1),
  fHfCoincidenceFilterBit(1),
  fClusterCompatibilityFilterBit(1),
  fnJets(0),
  fnGenJets(0),
  fnCaloJets(0),
  fEventWeight(1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetWTAPhiArray(),
  fJetEtaArray(),
  fJetWTAEtaArray(),
  fJetRawPtArray(),
  fJetRefPtArray(),
  fJetRefEtaArray(),
  fJetRefPhiArray(),
  fJetRefFlavorArray(),
  fGenJetPtArray(),
  fGenJetPhiArray(),
  fGenJetWTAPhiArray(),
  fGenJetEtaArray(),
  fGenJetWTAEtaArray(),
  fCaloJetPtArray(),
  fCaloJetPhiArray(),
  fCaloJetEtaArray(),
  fnTracks(0),
  fTrackPtVector(0),
  fTrackPtErrorVector(0),
  fTrackPhiVector(0),
  fTrackEtaVector(0),
  fHighPurityTrackVector(0),
  fTrackVertexDistanceZVector(0),
  fTrackVertexDistanceZErrorVector(0),
  fTrackVertexDistanceXYVector(0),
  fTrackVertexDistanceXYErrorVector(0),
  fTrackNormalizedChi2Vector(0),
  fnHitsTrackerLayerVector(0),
  fnHitsTrackVector(0),
  fTrackEnergyEcalVector(0),
  fTrackEnergyHcalVector(0),
  fTrackChargeVector(0),
  fnGenParticles(0),
  fGenParticlePtArray(0),
  fGenParticlePhiArray(0),
  fGenParticleEtaArray(0),
  fGenParticleChargeArray(0),
  fGenParticleSubeventArray(0)
{
  // Default constructor
  
  // Initialize fJetMaxTrackPtArray to -1
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetMaxTrackPtArray[i] = -1;
  }
  
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   Int_t jetType: 0 = Calo jets, 1 = CSPF jets, 2 = Flow subtracted CSPF jets
 *   Int_t jetAxis: 0 = E-scheme axis, 1 = WTA axis
 */
MonteCarloForestReader::MonteCarloForestReader(Int_t jetType, Int_t jetAxis) :
  fJetType(jetType),
  fJetAxis(jetAxis),
  fHeavyIonTree(0),
  fSkimTree(0),
  fJetTree(0),
  fTrackTree(0),
  fGenParticleTree(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fPtHatBranch(0),
  fEventWeightBranch(0),
  fPrimaryVertexBranch(0),
  fHfCoincidenceBranch(0),
  fClusterCompatibilityBranch(0),
  fnJetsBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetWTAPhiBranch(0),
  fJetEtaBranch(0),
  fJetWTAEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fJetRefPtBranch(0),
  fJetRefEtaBranch(0),
  fJetRefPhiBranch(0),
  fJetRefFlavorBranch(0),
  fnGenJetsBranch(0),
  fGenJetPtBranch(0),
  fGenJetPhiBranch(0),
  fGenJetWTAPhiBranch(0),
  fGenJetEtaBranch(0),
  fGenJetWTAEtaBranch(0),
  fnCaloJetsBranch(0),
  fCaloJetPtBranch(0),
  fCaloJetPhiBranch(0),
  fCaloJetEtaBranch(0),
  fnTracksBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fHighPurityTrackBranch(0),
  fTrackVertexDistanceZBranch(0),
  fTrackVertexDistanceZErrorBranch(0),
  fTrackVertexDistanceXYBranch(0),
  fTrackVertexDistanceXYErrorBranch(0),
  fTrackChi2Branch(0),
  fnTrackDegreesOfFreedomBranch(0),
  fnHitsTrackerLayerBranch(0),
  fnHitsTrackBranch(0),
  fTrackEnergyEcalBranch(0),
  fTrackEnergyHcalBranch(0),
  fTrackChargeBranch(0),
  fGenParticlePtBranch(0),
  fGenParticlePhiBranch(0),
  fGenParticleEtaBranch(0),
  fGenParticleChargeBranch(0),
  fGenParticleSubeventBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fPtHat(0),
  fPrimaryVertexFilterBit(1),
  fHfCoincidenceFilterBit(1),
  fClusterCompatibilityFilterBit(1),
  fnJets(0),
  fnGenJets(0),
  fnCaloJets(0),
  fEventWeight(1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetWTAPhiArray(),
  fJetEtaArray(),
  fJetWTAEtaArray(),
  fJetRawPtArray(),
  fJetRefPtArray(),
  fJetRefEtaArray(),
  fJetRefPhiArray(),
  fJetRefFlavorArray(),
  fGenJetPtArray(),
  fGenJetPhiArray(),
  fGenJetWTAPhiArray(),
  fGenJetEtaArray(),
  fGenJetWTAEtaArray(),
  fCaloJetPtArray(),
  fCaloJetPhiArray(),
  fCaloJetEtaArray(),
  fnTracks(0),
  fTrackPtVector(0),
  fTrackPtErrorVector(0),
  fTrackPhiVector(0),
  fTrackEtaVector(0),
  fHighPurityTrackVector(0),
  fTrackVertexDistanceZVector(0),
  fTrackVertexDistanceZErrorVector(0),
  fTrackVertexDistanceXYVector(0),
  fTrackVertexDistanceXYErrorVector(0),
  fTrackNormalizedChi2Vector(0),
  fnHitsTrackerLayerVector(0),
  fnHitsTrackVector(0),
  fTrackEnergyEcalVector(0),
  fTrackEnergyHcalVector(0),
  fTrackChargeVector(0),
  fnGenParticles(0),
  fGenParticlePtArray(0),
  fGenParticlePhiArray(0),
  fGenParticleEtaArray(0),
  fGenParticleChargeArray(0),
  fGenParticleSubeventArray(0)
{
  // Custom constructor
  
  // Initialize fJetMaxTrackPtArray to -1
  for(int i = 0; i < fnMaxJet; i++){
    fJetMaxTrackPtArray[i] = -1;
  }
  
}

/*
 * Copy constructor
 */
MonteCarloForestReader::MonteCarloForestReader(const MonteCarloForestReader& in) :
  fJetType(in.fJetType),
  fJetAxis(in.fJetAxis),
  fHeavyIonTree(in.fHeavyIonTree),
  fSkimTree(in.fSkimTree),
  fJetTree(in.fJetTree),
  fTrackTree(in.fTrackTree),
  fGenParticleTree(in.fGenParticleTree),
  fHiVzBranch(in.fHiVzBranch),
  fHiBinBranch(in.fHiBinBranch),
  fPtHatBranch(in.fPtHatBranch),
  fEventWeightBranch(in.fEventWeightBranch),
  fPrimaryVertexBranch(in.fPrimaryVertexBranch),
  fHfCoincidenceBranch(in.fHfCoincidenceBranch),
  fClusterCompatibilityBranch(in.fClusterCompatibilityBranch),
  fnJetsBranch(in.fnJetsBranch),
  fJetPtBranch(in.fJetPtBranch),
  fJetPhiBranch(in.fJetPhiBranch),
  fJetWTAPhiBranch(in.fJetWTAPhiBranch),
  fJetEtaBranch(in.fJetEtaBranch),
  fJetWTAEtaBranch(in.fJetWTAEtaBranch),
  fJetRawPtBranch(in.fJetRawPtBranch),
  fJetMaxTrackPtBranch(in.fJetMaxTrackPtBranch),
  fJetRefPtBranch(in.fJetRefPtBranch),
  fJetRefEtaBranch(in.fJetRefEtaBranch),
  fJetRefPhiBranch(in.fJetRefPhiBranch),
  fJetRefFlavorBranch(in.fJetRefFlavorBranch),
  fnGenJetsBranch(in.fnGenJetsBranch),
  fGenJetPtBranch(in.fGenJetPtBranch),
  fGenJetPhiBranch(in.fGenJetPhiBranch),
  fGenJetWTAPhiBranch(in.fGenJetWTAPhiBranch),
  fGenJetEtaBranch(in.fGenJetEtaBranch),
  fGenJetWTAEtaBranch(in.fGenJetWTAEtaBranch),
  fnCaloJetsBranch(in.fnCaloJetsBranch),
  fCaloJetPtBranch(in.fCaloJetPtBranch),
  fCaloJetPhiBranch(in.fCaloJetPhiBranch),
  fCaloJetEtaBranch(in.fCaloJetEtaBranch),
  fnTracksBranch(in.fnTracksBranch),
  fTrackPtBranch(in.fTrackPtBranch),
  fTrackPtErrorBranch(in.fTrackPtErrorBranch),
  fTrackPhiBranch(in.fTrackPhiBranch),
  fTrackEtaBranch(in.fTrackEtaBranch),
  fHighPurityTrackBranch(in.fHighPurityTrackBranch),
  fTrackVertexDistanceZBranch(in.fTrackVertexDistanceZBranch),
  fTrackVertexDistanceZErrorBranch(in.fTrackVertexDistanceZErrorBranch),
  fTrackVertexDistanceXYBranch(in.fTrackVertexDistanceXYBranch),
  fTrackVertexDistanceXYErrorBranch(in.fTrackVertexDistanceXYErrorBranch),
  fTrackChi2Branch(in.fTrackChi2Branch),
  fnTrackDegreesOfFreedomBranch(in.fnTrackDegreesOfFreedomBranch),
  fnHitsTrackerLayerBranch(in.fnHitsTrackerLayerBranch),
  fnHitsTrackBranch(in.fnHitsTrackBranch),
  fTrackEnergyEcalBranch(in.fTrackEnergyEcalBranch),
  fTrackEnergyHcalBranch(in.fTrackEnergyHcalBranch),
  fTrackChargeBranch(in.fTrackChargeBranch),
  fGenParticlePtBranch(in.fGenParticlePtBranch),
  fGenParticlePhiBranch(in.fGenParticlePhiBranch),
  fGenParticleEtaBranch(in.fGenParticleEtaBranch),
  fGenParticleChargeBranch(in.fGenParticleChargeBranch),
  fGenParticleSubeventBranch(in.fGenParticleSubeventBranch),
  fVertexZ(in.fVertexZ),
  fHiBin(in.fHiBin),
  fPtHat(in.fPtHat),
  fPrimaryVertexFilterBit(in.fPrimaryVertexFilterBit),
  fHfCoincidenceFilterBit(in.fHfCoincidenceFilterBit),
  fClusterCompatibilityFilterBit(in.fClusterCompatibilityFilterBit),
  fnJets(in.fnJets),
  fnGenJets(in.fnGenJets),
  fnCaloJets(in.fnCaloJets),
  fEventWeight(in.fEventWeight),
  fnTracks(in.fnTracks),
  fTrackPtVector(in.fTrackPtVector),
  fTrackPhiVector(in.fTrackPhiVector),
  fTrackEtaVector(in.fTrackEtaVector),
  fHighPurityTrackVector(in.fHighPurityTrackVector),
  fTrackVertexDistanceZVector(in.fTrackVertexDistanceZVector),
  fTrackVertexDistanceZErrorVector(in.fTrackVertexDistanceZErrorVector),
  fTrackVertexDistanceXYVector(in.fTrackVertexDistanceXYVector),
  fTrackVertexDistanceXYErrorVector(in.fTrackVertexDistanceXYErrorVector),
  fTrackNormalizedChi2Vector(in.fTrackNormalizedChi2Vector),
  fnHitsTrackerLayerVector(in.fnHitsTrackerLayerVector),
  fnHitsTrackVector(in.fnHitsTrackVector),
  fTrackEnergyEcalVector(in.fTrackEnergyEcalVector),
  fTrackEnergyHcalVector(in.fTrackEnergyHcalVector),
  fTrackChargeVector(in.fTrackChargeVector),
  fnGenParticles(in.fnGenParticles),
  fGenParticlePtArray(in.fGenParticlePtArray),
  fGenParticlePhiArray(in.fGenParticlePhiArray),
  fGenParticleEtaArray(in.fGenParticleEtaArray),
  fGenParticleChargeArray(in.fGenParticleChargeArray),
  fGenParticleSubeventArray(in.fGenParticleSubeventArray)
{
  // Copy constructor
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetWTAPhiArray[i] = in.fJetWTAPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetWTAEtaArray[i] = in.fJetWTAEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
    fJetRefPtArray[i] = in.fJetRefPtArray[i];
    fJetRefEtaArray[i] = in.fJetRefEtaArray[i];
    fJetRefPhiArray[i] = in.fJetRefPhiArray[i];
    fJetRefFlavorArray[i] = in.fJetRefFlavorArray[i];
    fGenJetPtArray[i] = in.fGenJetPtArray[i];
    fGenJetPhiArray[i] = in.fGenJetPhiArray[i];
    fGenJetWTAPhiArray[i] = in.fGenJetWTAPhiArray[i];
    fGenJetEtaArray[i] = in.fGenJetEtaArray[i];
    fGenJetWTAEtaArray[i] = in.fGenJetWTAEtaArray[i];
    fCaloJetPtArray[i] = in.fCaloJetPtArray[i];
    fCaloJetPhiArray[i] = in.fCaloJetPhiArray[i];
    fCaloJetEtaArray[i] = in.fCaloJetEtaArray[i];
  }
}

/*
 * Assignment operator
 */
MonteCarloForestReader& MonteCarloForestReader::operator=(const MonteCarloForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  fJetType = in.fJetType;
  fJetAxis = in.fJetAxis;
  fHeavyIonTree = in.fHeavyIonTree;
  fSkimTree = in.fSkimTree;
  fJetTree = in.fJetTree;
  fTrackTree = in.fTrackTree;
  fGenParticleTree = in.fGenParticleTree;
  fHiVzBranch = in.fHiVzBranch;
  fHiBinBranch = in.fHiBinBranch;
  fPtHatBranch = in.fPtHatBranch;
  fEventWeightBranch = in.fEventWeightBranch;
  fPrimaryVertexBranch = in.fPrimaryVertexBranch;
  fHfCoincidenceBranch = in.fHfCoincidenceBranch;
  fClusterCompatibilityBranch = in.fClusterCompatibilityBranch;
  fnJetsBranch = in.fnJetsBranch;
  fJetPtBranch = in.fJetPtBranch;
  fJetPhiBranch = in.fJetPhiBranch;
  fJetWTAPhiBranch = in.fJetWTAPhiBranch;
  fJetEtaBranch = in.fJetEtaBranch;
  fJetWTAEtaBranch = in.fJetWTAEtaBranch;
  fJetRawPtBranch = in.fJetRawPtBranch;
  fJetMaxTrackPtBranch = in.fJetMaxTrackPtBranch;
  fJetRefPtBranch = in.fJetRefPtBranch;
  fJetRefEtaBranch = in.fJetRefEtaBranch;
  fJetRefPhiBranch = in.fJetRefPhiBranch;
  fJetRefFlavorBranch = in.fJetRefFlavorBranch;
  fnGenJetsBranch = in.fnGenJetsBranch;
  fGenJetPtBranch = in.fGenJetPtBranch;
  fGenJetPhiBranch = in.fGenJetPhiBranch;
  fGenJetWTAPhiBranch = in.fGenJetWTAPhiBranch;
  fGenJetEtaBranch = in.fGenJetEtaBranch;
  fGenJetWTAEtaBranch = in.fGenJetWTAEtaBranch;
  fnCaloJetsBranch = in.fnCaloJetsBranch;
  fCaloJetPtBranch = in.fCaloJetPtBranch;
  fCaloJetPhiBranch = in.fCaloJetPhiBranch;
  fCaloJetEtaBranch = in.fCaloJetEtaBranch;
  fnTracksBranch = in.fnTracksBranch;
  fTrackPtBranch = in.fTrackPtBranch;
  fTrackPtErrorBranch = in.fTrackPtErrorBranch;
  fTrackPhiBranch = in.fTrackPhiBranch;
  fTrackEtaBranch = in.fTrackEtaBranch;
  fHighPurityTrackBranch = in.fHighPurityTrackBranch;
  fTrackVertexDistanceZBranch = in.fTrackVertexDistanceZBranch;
  fTrackVertexDistanceZErrorBranch = in.fTrackVertexDistanceZErrorBranch;
  fTrackVertexDistanceXYBranch = in.fTrackVertexDistanceXYBranch;
  fTrackVertexDistanceXYErrorBranch = in.fTrackVertexDistanceXYErrorBranch;
  fTrackChi2Branch = in.fTrackChi2Branch;
  fnTrackDegreesOfFreedomBranch = in.fnTrackDegreesOfFreedomBranch;
  fnHitsTrackerLayerBranch = in.fnHitsTrackerLayerBranch;
  fnHitsTrackBranch = in.fnHitsTrackBranch;
  fTrackEnergyEcalBranch = in.fTrackEnergyEcalBranch;
  fTrackEnergyHcalBranch = in.fTrackEnergyHcalBranch;
  fTrackChargeBranch = in.fTrackChargeBranch;
  fGenParticlePtBranch = in.fGenParticlePtBranch;
  fGenParticlePhiBranch = in.fGenParticlePhiBranch;
  fGenParticleEtaBranch = in.fGenParticleEtaBranch;
  fGenParticleChargeBranch = in.fGenParticleChargeBranch;
  fGenParticleSubeventBranch = in.fGenParticleSubeventBranch;
  fVertexZ = in.fVertexZ;
  fHiBin = in.fHiBin;
  fPtHat = in.fPtHat;
  fPrimaryVertexFilterBit = in.fPrimaryVertexFilterBit;
  fHfCoincidenceFilterBit = in.fHfCoincidenceFilterBit;
  fClusterCompatibilityFilterBit = in.fClusterCompatibilityFilterBit;
  fnJets = in.fnJets;
  fnGenJets = in.fnGenJets;
  fnCaloJets = in.fnCaloJets;
  fEventWeight = in.fEventWeight;
  fnTracks = in.fnTracks;
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetWTAPhiArray[i] = in.fJetWTAPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetWTAEtaArray[i] = in.fJetWTAEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
    fJetRefPtArray[i] = in.fJetRefPtArray[i];
    fJetRefEtaArray[i] = in.fJetRefEtaArray[i];
    fJetRefPhiArray[i] = in.fJetRefPhiArray[i];
    fJetRefFlavorArray[i] = in.fJetRefFlavorArray[i];
    fGenJetPtArray[i] = in.fGenJetPtArray[i];
    fGenJetPhiArray[i] = in.fGenJetPhiArray[i];
    fGenJetWTAPhiArray[i] = in.fGenJetWTAPhiArray[i];
    fGenJetEtaArray[i] = in.fGenJetEtaArray[i];
    fGenJetWTAEtaArray[i] = in.fGenJetWTAEtaArray[i];
    fCaloJetPtArray[i] = in.fCaloJetPtArray[i];
    fCaloJetPhiArray[i] = in.fCaloJetPhiArray[i];
    fCaloJetEtaArray[i] = in.fCaloJetEtaArray[i];
  }
  
  // Copy the track vectors
  fTrackPtVector = in.fTrackPtVector;
  fTrackPtVector = in.fTrackPtVector;
  fTrackPhiVector = in.fTrackPhiVector;
  fTrackEtaVector = in.fTrackEtaVector;
  fHighPurityTrackVector = in.fHighPurityTrackVector;
  fTrackVertexDistanceZVector = in.fTrackVertexDistanceZVector;
  fTrackVertexDistanceZErrorVector = in.fTrackVertexDistanceZErrorVector;
  fTrackVertexDistanceXYVector = in.fTrackVertexDistanceXYVector;
  fTrackVertexDistanceXYErrorVector = in.fTrackVertexDistanceXYErrorVector;
  fTrackNormalizedChi2Vector = in.fTrackNormalizedChi2Vector;
  fnHitsTrackerLayerVector = in.fnHitsTrackerLayerVector;
  fnHitsTrackVector = in.fnHitsTrackVector;
  fTrackEnergyEcalVector = in.fTrackEnergyEcalVector;
  fTrackEnergyHcalVector = in.fTrackEnergyHcalVector;
  fTrackChargeVector = in.fTrackChargeVector;
  
  // Copy the generator level particle vectors
  fnGenParticles = in.fnGenParticles;
  fGenParticlePtArray = in.fGenParticlePtArray;
  fGenParticlePhiArray = in.fGenParticlePhiArray;
  fGenParticleEtaArray = in.fGenParticleEtaArray;
  fGenParticleChargeArray = in.fGenParticleChargeArray;
  fGenParticleSubeventArray = in.fGenParticleSubeventArray;
  
  return *this;
}

/*
 * Destructor
 */
MonteCarloForestReader::~MonteCarloForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void MonteCarloForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("hiBin",1);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  fHeavyIonTree->SetBranchStatus("pthat",1);
  fHeavyIonTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
  fHeavyIonTree->SetBranchStatus("weight",1);
  fHeavyIonTree->SetBranchAddress("weight",&fEventWeight,&fEventWeightBranch); // event weight only for MC

  // Connect the branches of the skim tree
  fSkimTree->SetBranchStatus("*",0);
  fSkimTree->SetBranchStatus("pprimaryVertexFilter",1);
  fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
  fSkimTree->SetBranchStatus("pphfCoincFilter2Th4",1); 
  fSkimTree->SetBranchAddress("pphfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
  fSkimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
  fSkimTree->SetBranchAddress("pclusterCompatibilityFilter", &fClusterCompatibilityFilterBit, &fClusterCompatibilityBranch);
  
  // Connect the branches to the jet tree  
  fJetTree->SetBranchStatus("*",0);
  fJetTree->SetBranchStatus("jtpt",1);
  fJetTree->SetBranchAddress("jtpt",&fJetPtArray,&fJetPtBranch);
  
  // Load jet phi with E-scheme and WTA axes
  fJetTree->SetBranchStatus("jtphi",1);
  fJetTree->SetBranchAddress("jtphi",&fJetPhiArray,&fJetPhiBranch);
  fJetTree->SetBranchStatus("WTAphi",1);
  fJetTree->SetBranchAddress("WTAphi",&fJetWTAPhiArray,&fJetWTAPhiBranch);
  
  // Load jet eta with E-scheme and WTA axes
  fJetTree->SetBranchStatus("jteta",1);
  fJetTree->SetBranchAddress("jteta",&fJetEtaArray,&fJetEtaBranch);
  fJetTree->SetBranchStatus("WTAeta",1);
  fJetTree->SetBranchAddress("WTAeta",&fJetWTAEtaArray,&fJetWTAEtaBranch);
  
  fJetTree->SetBranchStatus("nref",1);
  fJetTree->SetBranchAddress("nref",&fnJets,&fnJetsBranch);
  fJetTree->SetBranchStatus("rawpt",1);
  fJetTree->SetBranchAddress("rawpt",&fJetRawPtArray,&fJetRawPtBranch);
  fJetTree->SetBranchStatus("trackMax",1);
  fJetTree->SetBranchAddress("trackMax",&fJetMaxTrackPtArray,&fJetMaxTrackPtBranch);

  // Connect the reference jet and gen jet arrays
  fJetTree->SetBranchStatus("refpt",1);
  fJetTree->SetBranchAddress("refpt",&fJetRefPtArray,&fJetRefPtBranch);
  fJetTree->SetBranchStatus("refeta",1);
  fJetTree->SetBranchAddress("refeta",&fJetRefEtaArray,&fJetRefEtaBranch);
  fJetTree->SetBranchStatus("refphi",1);
  fJetTree->SetBranchAddress("refphi",&fJetRefPhiArray,&fJetRefPhiBranch);
  fJetTree->SetBranchStatus("matchedPartonFlavor",1);
  fJetTree->SetBranchAddress("matchedPartonFlavor",&fJetRefFlavorArray,&fJetRefFlavorBranch);

  fJetTree->SetBranchStatus("genpt",1);
  fJetTree->SetBranchAddress("genpt",&fGenJetPtArray,&fGenJetPtBranch);
    
  // Load jet phi with E-scheme and WTA axes
  fJetTree->SetBranchStatus("genphi",1);
  fJetTree->SetBranchAddress("genphi",&fGenJetPhiArray,&fGenJetPhiBranch);
  fJetTree->SetBranchStatus("WTAgenphi",1);
  fJetTree->SetBranchAddress("WTAgenphi",&fGenJetWTAPhiArray,&fGenJetWTAPhiBranch);
    
  // Load jet eta with E-scheme and WTA axes
  fJetTree->SetBranchStatus("geneta",1);
  fJetTree->SetBranchAddress("geneta",&fGenJetEtaArray,&fGenJetEtaBranch);
  fJetTree->SetBranchStatus("WTAgeneta",1);
  fJetTree->SetBranchAddress("WTAgeneta",&fGenJetWTAEtaArray,&fGenJetWTAEtaBranch);
    
  fJetTree->SetBranchStatus("ngen",1);
  fJetTree->SetBranchAddress("ngen",&fnGenJets,&fnGenJetsBranch);

  // Load the variables for calo jets
  fJetTree->SetBranchStatus("ncalo", 1);
  fJetTree->SetBranchAddress("ncalo", &fnCaloJets, &fnCaloJetsBranch);
  fJetTree->SetBranchStatus("calopt", 1);
  fJetTree->SetBranchAddress("calopt", &fCaloJetPtArray, &fCaloJetPtBranch);
  fJetTree->SetBranchStatus("calophi", 1);
  fJetTree->SetBranchAddress("calophi", &fCaloJetPhiArray, &fCaloJetPhiBranch);
  fJetTree->SetBranchStatus("caloeta", 1);
  fJetTree->SetBranchAddress("caloeta", &fCaloJetEtaArray, &fCaloJetEtaBranch);
  
  // Connect the branches to the track tree
  /*
  
  fTrackTree->SetBranchStatus("*",0);
  
  // Read the track vectors    
  fTrackTree->SetBranchStatus("trkPt",1);
  fTrackTree->SetBranchAddress("trkPt",&fTrackPtVector,&fTrackPtBranch);
  fTrackTree->SetBranchStatus("trkPtError",1);
  fTrackTree->SetBranchAddress("trkPtError",&fTrackPtErrorVector,&fTrackPtErrorBranch);
  fTrackTree->SetBranchStatus("trkPhi",1);
  fTrackTree->SetBranchAddress("trkPhi",&fTrackPhiVector,&fTrackPhiBranch);
  fTrackTree->SetBranchStatus("trkEta",1);
  fTrackTree->SetBranchAddress("trkEta",&fTrackEtaVector,&fTrackEtaBranch);
  fTrackTree->SetBranchStatus("nTrk",1);
  fTrackTree->SetBranchAddress("nTrk",&fnTracks,&fnTracksBranch);
  fTrackTree->SetBranchStatus("highPurity",1);
  fTrackTree->SetBranchAddress("highPurity",&fHighPurityTrackVector,&fHighPurityTrackBranch);
  fTrackTree->SetBranchStatus("trkDzFirstVtx",1);
  fTrackTree->SetBranchAddress("trkDzFirstVtx",&fTrackVertexDistanceZVector,&fTrackVertexDistanceZBranch);
  fTrackTree->SetBranchStatus("trkDzErrFirstVtx",1);
  fTrackTree->SetBranchAddress("trkDzErrFirstVtx",&fTrackVertexDistanceZErrorVector,&fTrackVertexDistanceZErrorBranch);
  fTrackTree->SetBranchStatus("trkDxyFirstVtx",1);
  fTrackTree->SetBranchAddress("trkDxyFirstVtx",&fTrackVertexDistanceXYVector,&fTrackVertexDistanceXYBranch);
  fTrackTree->SetBranchStatus("trkDxyErrFirstVtx",1);
  fTrackTree->SetBranchAddress("trkDxyErrFirstVtx",&fTrackVertexDistanceXYErrorVector,&fTrackVertexDistanceXYErrorBranch);
  fTrackTree->SetBranchStatus("trkNormChi2",1);
  fTrackTree->SetBranchAddress("trkNormChi2",&fTrackNormalizedChi2Vector,&fTrackChi2Branch);
  fTrackTree->SetBranchStatus("trkNLayers",1);
  fTrackTree->SetBranchAddress("trkNLayers",&fnHitsTrackerLayerVector,&fnHitsTrackerLayerBranch);
  fTrackTree->SetBranchStatus("trkNHits",1);
  fTrackTree->SetBranchAddress("trkNHits",&fnHitsTrackVector,&fnHitsTrackBranch);
  fTrackTree->SetBranchStatus("pfEcal",1);
  fTrackTree->SetBranchAddress("pfEcal",&fTrackEnergyEcalVector,&fTrackEnergyEcalBranch);
  fTrackTree->SetBranchStatus("pfHcal",1);
  fTrackTree->SetBranchAddress("pfHcal",&fTrackEnergyHcalVector,&fTrackEnergyHcalBranch);
  fTrackTree->SetBranchStatus("trkCharge",1);
  fTrackTree->SetBranchAddress("trkCharge",&fTrackChargeVector,&fTrackChargeBranch);
  */
  
  // Connect the branches to the generator level particle tree
  fGenParticleTree->SetBranchStatus("*",0);
  fGenParticleTree->SetBranchStatus("pt",1);
  fGenParticleTree->SetBranchAddress("pt",&fGenParticlePtArray,&fGenParticlePtBranch);
  fGenParticleTree->SetBranchStatus("phi",1);
  fGenParticleTree->SetBranchAddress("phi",&fGenParticlePhiArray,&fGenParticlePhiBranch);
  fGenParticleTree->SetBranchStatus("eta",1);
  fGenParticleTree->SetBranchAddress("eta",&fGenParticleEtaArray,&fGenParticleEtaBranch);
  fGenParticleTree->SetBranchStatus("chg",1);
  fGenParticleTree->SetBranchAddress("chg",&fGenParticleChargeArray,&fGenParticleChargeBranch);
  fGenParticleTree->SetBranchStatus("sube",1);
  fGenParticleTree->SetBranchAddress("sube",&fGenParticleSubeventArray,&fGenParticleSubeventBranch);
  
}


/*
 * Connect a new tree to the reader
 */
void MonteCarloForestReader::ReadForestFromFile(TFile* inputFile){
  
  // Connect a trees from the file to the reader
  fHeavyIonTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
  fSkimTree = (TTree*)inputFile->Get("skimanalysis/HltTree");
  
  // Possible jet trees to be read
  const char *treeName[3] = {"none","none","none"};
  treeName[0] = "akPu4CaloJetAnalyzer/t";     // Tree for calo jets
  treeName[1] = "akCs4PFJetAnalyzer/t";       // Tree for csPF jets
  treeName[2] = "akFlowPuCs4PFJetAnalyzer/t"; // Tree for flow subtracted csPF jets
  
  fJetTree = (TTree*)inputFile->Get(treeName[fJetType]);
  
  // Read track and generator level particle trees
  //fTrackTree = (TTree*)inputFile->Get("PbPbTracks/trackTree");
  fGenParticleTree = (TTree*)inputFile->Get("HiGenParticleAna/hi");
  
  Initialize();
}

/*
 * Connect a new tree to the reader
 */
void MonteCarloForestReader::ReadForestFromFileList(std::vector<TString> fileList){
  TFile *inputFile = TFile::Open(fileList.at(0));
  ReadForestFromFile(inputFile);
}

/*
 * Burn the current forest.
 */
void MonteCarloForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fSkimTree->Delete();
  fJetTree->Delete();
  //fTrackTree->Delete();
  fGenParticleTree->Delete();
}

/*
 * Load an event to memory
 */
void MonteCarloForestReader::GetEvent(Int_t iEvent){
  fHeavyIonTree->GetEntry(iEvent);
  fSkimTree->GetEntry(iEvent);
  fJetTree->GetEntry(iEvent);
  //fTrackTree->GetEntry(iEvent);
  fGenParticleTree->GetEntry(iEvent);
   
  // Read the numbers of generator level particles for this event
  fnGenParticles = fGenParticlePtArray->size();
}

// Getter for number of events in the tree
Int_t MonteCarloForestReader::GetNEvents() const{
  return fJetPtBranch->GetEntries();
}

// Getter for number of jets
Int_t MonteCarloForestReader::GetNJets(Int_t jetType) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetNJets();
    case kGeneratorLevelJet:
      return GetNGeneratorJets();
    default:
      return -999;
  }
}

// Getter for jet pT
Float_t MonteCarloForestReader::GetJetPt(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetJetPt(iJet);
    case kGeneratorLevelJet:
      return GetGeneratorJetPt(iJet);
    default:
      return -999;
  }
} 

// Getter for jet phi
Float_t MonteCarloForestReader::GetJetPhi(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetJetPhi(iJet);
    case kGeneratorLevelJet:
      return GetGeneratorJetPhi(iJet);
    default:
      return -999;
  }
}

// Getter for jet eta
Float_t MonteCarloForestReader::GetJetEta(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetJetEta(iJet);
    case kGeneratorLevelJet:
      return GetGeneratorJetEta(iJet);
    default:
      return -999;
  }
}

// Getter for jet raw pT
Float_t MonteCarloForestReader::GetJetRawPt(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetJetRawPt(iJet);
    case kGeneratorLevelJet:
      return GetGeneratorJetPt(iJet); // Raw pT is meaningless for generator level jets. Just return regular pT.
    default:
      return -999;
  }
} 

// Getter for maximum track pT inside a jet
Float_t MonteCarloForestReader::GetJetMaxTrackPt(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetJetMaxTrackPt(iJet);
    case kGeneratorLevelJet:
      return -999; // Track max is only relevant for reconstructed jet quality cuts. Return -999 as a nonsensical value here.
    default:
      return -999;
  }
}

// Check if there is matching reconstructed/generator level jet
Bool_t MonteCarloForestReader::HasMatchingJet(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return HasMatchingGenJet(iJet);
    case kGeneratorLevelJet:
      return HasMatchingRecoJet(iJet);
    default:
      return false;
  }
} 

// Get the index of the matched jet
Int_t MonteCarloForestReader::GetMatchingIndex(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetMatchingGenIndex(iJet);
    case kGeneratorLevelJet:
      return GetMatchingRecoIndex(iJet);
    default:
      return -999;
  }
}

// Getter for matched jet pT
Float_t MonteCarloForestReader::GetMatchedPt(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetMatchedGenPt(iJet);
    case kGeneratorLevelJet:
      return GetMatchedRecoPt(iJet);
    default:
      return -999;
  }
}

// Getter for matched jet eta
Float_t MonteCarloForestReader::GetMatchedEta(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetMatchedGenEta(iJet);
    case kGeneratorLevelJet:
      return GetMatchedRecoEta(iJet);
    default:
      return -999;
  }
}

// Getter for matched jet phi
Float_t MonteCarloForestReader::GetMatchedPhi(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetMatchedGenPhi(iJet);
    case kGeneratorLevelJet:
      return GetMatchedRecoPhi(iJet);
    default:
      return -999;
  }
}

// Getter for the jet flavor for input jet type
Int_t MonteCarloForestReader::GetJetFlavor(Int_t jetType, Int_t iJet) const{
  switch (jetType) {
    case kReconstructedJet:
      return GetRecoJetFlavor(iJet);
    case kGeneratorLevelJet:
      return GetGenJetFlavor(iJet);
    default:
      return -999;
  }
}

// Getter for number of jets in an event
Int_t MonteCarloForestReader::GetNJets() const{
  return fnJets;
}

// Getter for number of generator level jets in an event
Int_t MonteCarloForestReader::GetNGeneratorJets() const{
  return fnGenJets;
}

// Getter for number of calorimeter jets in an event
Int_t MonteCarloForestReader::GetNCalorimeterJets() const{
  return fnCaloJets;
}


// Getter for jet pT
Float_t MonteCarloForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t MonteCarloForestReader::GetJetPhi(Int_t iJet) const{
  if(fJetAxis == 0) return fJetPhiArray[iJet];
  return fJetWTAPhiArray[iJet];
}

// Getter for jet eta
Float_t MonteCarloForestReader::GetJetEta(Int_t iJet) const{
  if(fJetAxis == 0) return fJetEtaArray[iJet];
  return fJetWTAEtaArray[iJet];
}

// Getter for jet raw pT
Float_t MonteCarloForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetRawPtArray[iJet];
}

// Getter for maximum track pT inside a jet
Float_t MonteCarloForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return fJetMaxTrackPtArray[iJet];
}

// Check if reconstructed jet has a matching generator level jet
Bool_t MonteCarloForestReader::HasMatchingGenJet(Int_t iJet) const{
  
  // For each reconstructed jet there is a reference pT, which tells the the pT of a matched generator level jet
  // If this number is -999, it means that there are no generator level jets matching the reconstructed jet
  if(fJetRefPtArray[iJet] < 0) return false;
  return true;
}

// Get the matching generator level jet index for the given reconstructed jet
Int_t MonteCarloForestReader::GetMatchingGenIndex(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, it is matched to the given reconstructed jet
  for(Int_t iGenJet = 0; iGenJet < fnGenJets; iGenJet++){
    if(TMath::Abs(fJetRefEtaArray[iJet] - fGenJetEtaArray[iGenJet]) < 0.015){
      if(TMath::Abs(fJetRefPhiArray[iJet] - fGenJetPhiArray[iGenJet]) < 0.015){
        if(TMath::Abs(fJetRefPtArray[iJet] - fGenJetPtArray[iGenJet]) < 0.03*fGenJetPtArray[iGenJet]){
          return iGenJet;
        }
      }
    }
  }
  // If a matching index is not found, return -1 to show that
  return -1;  
}

// Getter for matched generator level jet pT
Float_t MonteCarloForestReader::GetMatchedGenPt(Int_t iJet) const{
  
  // Find the index of the matching generator level jet
  Int_t matchedIndex = GetMatchingGenIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;

  // Return matched gen pT
  return fGenJetPtArray[matchedIndex];
}

// Get the eta of the matched generator level jet
Float_t MonteCarloForestReader::GetMatchedGenEta(Int_t iJet) const{

  // Find the index of the matching generator level jet
  Int_t matchedIndex = GetMatchingGenIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  if(fJetAxis == 0) return fGenJetEtaArray[matchedIndex];
  return fGenJetWTAEtaArray[matchedIndex];
  
}

// Get the phi of the matched reconstructed jet
Float_t MonteCarloForestReader::GetMatchedGenPhi(Int_t iJet) const{
  
  // Find the index of the matching generator level jet
  Int_t matchedIndex = GetMatchingGenIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  if(fJetAxis == 0) return fGenJetPhiArray[matchedIndex];
  return fGenJetWTAPhiArray[matchedIndex];
  
}

// Getter for generator level jet pT
Float_t MonteCarloForestReader::GetGeneratorJetPt(Int_t iJet) const{
  return fGenJetPtArray[iJet];
}

// Getter for generator level jet phi
Float_t MonteCarloForestReader::GetGeneratorJetPhi(Int_t iJet) const{
  if(fJetAxis == 0) return fGenJetPhiArray[iJet];
  return fGenJetWTAPhiArray[iJet];
}

// Getter for generator level jet eta
Float_t MonteCarloForestReader::GetGeneratorJetEta(Int_t iJet) const{
  if(fJetAxis == 0) return fGenJetEtaArray[iJet];
  return fGenJetWTAEtaArray[iJet];
}

// Getter for calorimeter jet pT
Float_t MonteCarloForestReader::GetCalorimeterJetPt(Int_t iJet) const{
  return fCaloJetPtArray[iJet];
}

// Getter for calorimeter jet phi
Float_t MonteCarloForestReader::GetCalorimeterJetPhi(Int_t iJet) const{
  return fCaloJetPhiArray[iJet];
}

// Getter for calorimeter jet eta
Float_t MonteCarloForestReader::GetCalorimeterJetEta(Int_t iJet) const{
  return fCaloJetEtaArray[iJet];
}

// Check if generator level jet has a matching reconstructed jet
Bool_t MonteCarloForestReader::HasMatchingRecoJet(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, check also eta and phi
  // If all values are close by, we must have a matching jet
  Double_t jetPt = GetGeneratorJetPt(iJet);
  for(Int_t iRef = 0; iRef < fnJets; iRef++){
    if(TMath::Abs(fGenJetEtaArray[iJet] - fJetRefEtaArray[iRef]) < 0.015){
      if(TMath::Abs(fGenJetPhiArray[iJet] - fJetRefPhiArray[iRef]) < 0.015){
        if(TMath::Abs(jetPt - fJetRefPtArray[iRef]) < 0.03*jetPt) {
          return true;
        }
      }
    }
  }
  
  return false;
}

// Get the index of the matched reconstructed jet
Int_t MonteCarloForestReader::GetMatchingRecoIndex(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, check also eta and phi
  // If all values are close by, we must have a matching jet
  Double_t jetPt = GetGeneratorJetPt(iJet);
  for(Int_t iRef = 0; iRef < fnJets; iRef++){
    if(TMath::Abs(fGenJetEtaArray[iJet] - fJetRefEtaArray[iRef]) < 0.015){
      if(TMath::Abs(fGenJetPhiArray[iJet] - fJetRefPhiArray[iRef]) < 0.015){
        if(TMath::Abs(jetPt - fJetRefPtArray[iRef]) < 0.03*jetPt) {
          return iRef;
        }
      }
    }
  }
  
  // Return -1 to show a match is not found
  return -1;
}

// Get the pT of the matched reconstructed jet
Float_t MonteCarloForestReader::GetMatchedRecoPt(Int_t iJet) const{
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingRecoIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet pT. Need raw pT for reco jets before jet corrections are in the forest
  return fJetRawPtArray[matchingIndex]; // Use this if doing jet pT correction manually
  //return fJetPtArray[matchingIndex]; // Use this if jet pT corrected in the forest
}

// Get the phi of the matched reconstructed jet
Float_t MonteCarloForestReader::GetMatchedRecoPhi(Int_t iJet) const{
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingRecoIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet phi
  if(fJetAxis == 0) return fJetPhiArray[matchingIndex];
  return fJetWTAPhiArray[matchingIndex];
}

// Get the eta of the matched reconstructed jet
Float_t MonteCarloForestReader::GetMatchedRecoEta(Int_t iJet) const{
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingRecoIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet eta
  if(fJetAxis == 0) return fJetEtaArray[matchingIndex];
  return fJetWTAEtaArray[matchingIndex];
}

// Getter for vertex z position
Float_t MonteCarloForestReader::GetVz() const{
  return fVertexZ;
}

// Getter for centrality. CMS has integer centrality bins from 0 to 200, thus division by 2.
Float_t MonteCarloForestReader::GetCentrality() const{
  return fHiBin/2.0;
}

// Getter for hiBin. Return 1 for negative values (for easier handling of tracking efficiency correction)
Int_t MonteCarloForestReader::GetHiBin() const{
  if(fHiBin < 0) return 1;
  return fHiBin;
}

// Getter for pT hat
Float_t MonteCarloForestReader::GetPtHat() const{
  return fPtHat;
}

// Getter for pT hat
Float_t MonteCarloForestReader::GetEventWeight() const{
  return fEventWeight;
}

// Getter for number of tracks in an event
Int_t MonteCarloForestReader::GetNTracks() const{
  return fnTracks;
}

// Getter for track pT
Float_t MonteCarloForestReader::GetTrackPt(Int_t iTrack) const{
  return fTrackPtVector->at(iTrack);
}

// Getter for track pT error
Float_t MonteCarloForestReader::GetTrackPtError(Int_t iTrack) const{
  return fTrackPtErrorVector->at(iTrack);
}

// Getter for track phi
Float_t MonteCarloForestReader::GetTrackPhi(Int_t iTrack) const{
  return fTrackPhiVector->at(iTrack);
}

// Getter for track eta
Float_t MonteCarloForestReader::GetTrackEta(Int_t iTrack) const{
  return fTrackEtaVector->at(iTrack);
}

// Getter for high purity of the track
Bool_t MonteCarloForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return fHighPurityTrackVector->at(iTrack);
}

// Getter for track distance from primary vertex in z-direction
Float_t MonteCarloForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return fTrackVertexDistanceZVector->at(iTrack);
}

// Getter for error of track distance from primary vertex in z-direction
Float_t MonteCarloForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return fTrackVertexDistanceZErrorVector->at(iTrack);
}

// Getter for track distance from primary vertex in xy-direction
Float_t MonteCarloForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return fTrackVertexDistanceXYVector->at(iTrack);
}

// Getter for error of track distance from primary vertex in xy-direction
Float_t MonteCarloForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return fTrackVertexDistanceXYErrorVector->at(iTrack);
}

// Getter for normalized track chi2 value from reconstruction fit
Float_t MonteCarloForestReader::GetTrackNormalizedChi2(Int_t iTrack) const{
  return fTrackNormalizedChi2Vector->at(iTrack);
}

// Getter for number of hits in tracker layers
Int_t MonteCarloForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return fnHitsTrackerLayerVector->at(iTrack);
}

// Getter for number of hits for the track
Int_t MonteCarloForestReader::GetNHitsTrack(Int_t iTrack) const{
  return fnHitsTrackVector->at(iTrack);
}

// Getter for track energy in ECal
Float_t MonteCarloForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return fTrackEnergyEcalVector->at(iTrack);
}

// Getter for track energy in HCal
Float_t MonteCarloForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return fTrackEnergyHcalVector->at(iTrack);
}

// Getter for track charge
Int_t MonteCarloForestReader::GetTrackCharge(Int_t iTrack) const{
  return fTrackChargeVector->at(iTrack);
}

// Getter for number of generator level particles
Int_t MonteCarloForestReader::GetNGenParticles() const{
  return fnGenParticles;
}

// Getter for generator level particle pT
Float_t MonteCarloForestReader::GetGenParticlePt(Int_t iTrack) const{
  return fGenParticlePtArray->at(iTrack);
}

// Getter for generator level particle phi
Float_t MonteCarloForestReader::GetGenParticlePhi(Int_t iTrack) const{
  return fGenParticlePhiArray->at(iTrack);
}

// Getter for generator level particle eta
Float_t MonteCarloForestReader::GetGenParticleEta(Int_t iTrack) const{
  return fGenParticleEtaArray->at(iTrack);
}

// Getter for generator level particle charge
Int_t MonteCarloForestReader::GetGenParticleCharge(Int_t iTrack) const{
  return fGenParticleChargeArray->at(iTrack);
}

// Getter for generator level particle subevent index
Int_t MonteCarloForestReader::GetGenParticleSubevent(Int_t iTrack) const{
  return fGenParticleSubeventArray->at(iTrack);
}

// Getter for reconstructed jet flavor
Int_t MonteCarloForestReader::GetRecoJetFlavor(Int_t iJet) const{
  return fJetRefFlavorArray[iJet];
}

// Getter for generator level jet flavor
Int_t MonteCarloForestReader::GetGenJetFlavor(Int_t iJet) const{

  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingRecoIndex(iJet);
  
  // If we did not find match, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching parton flavor
  return fJetRefFlavorArray[matchingIndex];

}      

// Getter for primary vertex filter bit
Int_t MonteCarloForestReader::GetPrimaryVertexFilterBit() const{
  return fPrimaryVertexFilterBit;
}

// Getter for HF energy coincidence filter bit
Int_t MonteCarloForestReader::GetHfCoincidenceFilterBit() const{
  return fHfCoincidenceFilterBit;
}

// Getter for cluster compatibility filter bit
Int_t MonteCarloForestReader::GetClusterCompatibilityFilterBit() const{
  return fClusterCompatibilityFilterBit;
}