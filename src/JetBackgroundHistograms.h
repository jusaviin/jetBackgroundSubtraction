// Class for histograms needed in the jet background analysis

#ifndef JETBACKGROUNDHISTOGRAMS_H
#define JETBACKGROUNDHISTOGRAMS_H

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

// Own includes
#include "ConfigurationCard.h"

class JetBackgroundHistograms{
  
public:
  
  // Enumeration for event types to event histogram and track cuts for track cut histogram
  enum enumEventTypes {kAll, kPrimaryVertex, kHfCoincidence, kClusterCompatibility, kVzCut, knEventTypes};
  enum enumClosureParticleType {kQuark,kGluon,knClosureParticleTypes};
  enum enumEventPlaneOrder {kSecondOrderEventPlane, kThirdOrderEventPlane, kFourthOrderEventPlane, knEventPlanes};
    
  // Constructors and destructor
  JetBackgroundHistograms(); // Default constructor
  JetBackgroundHistograms(ConfigurationCard* newCard); // Custom constructor
  JetBackgroundHistograms(const JetBackgroundHistograms& in); // Copy constructor
  virtual ~JetBackgroundHistograms(); // Destructor
  JetBackgroundHistograms& operator=(const JetBackgroundHistograms& obj); // Equal sign operator
  
  // Methods
  void CreateHistograms();                   // Create all histograms
  void Write() const;                        // Write the histograms to a file that is opened somewhere else
  void Write(TString outputFileName) const;  // Write the histograms to a file
  void SetCard(ConfigurationCard* newCard);  // Set a new configuration card for the histogram class
  
  // Histograms defined public to allow easier access to them. Should not be abused
  TH1F* fhVertexZ;                 // Vertex z-position
  TH1F* fhVertexZWeighted;         // Weighted vertex z-position (only meaningfull for MC)
  TH1F* fhEvents;                  // Number of events. For binning see enumEventTypes.
  TH1F* fhCentrality;              // Centrality information. -0.5 for pp or PYTHIA.
  TH1F* fhCentralityWeighted;      // Weighted centrality distribution (only meaningful for MC)
  TH1F* fhPtHat;                   // pT hat for MC events (only meaningful for MC)
  TH1F* fhPtHatWeighted;           // Weighted pT hat distribution
  THnSparseF* fhInclusiveJet; // Inclusive jet information
  THnSparseF* fhLeadingJet;   // Leading jet information
  THnSparseF* fhJetPtClosure; // Jet pT closure histograms. Also information for response matrix.
  THnSparseF *fhInclusiveJetEventPlane[knEventPlanes];  // Correlation between jets and event plane angles
  THnSparseF *fhLeadingJetEventPlane[knEventPlanes];    // Correlation between leading jets and event plane angles

private:
  
  ConfigurationCard* fCard;    // Card for binning info
  const TString kEventTypeStrings[knEventTypes] = {"All", "PrimVertex", "HfCoin2Th4", "ClustCompt", "v_{z} cut"}; // Strings corresponding to event types
  
};

#endif
