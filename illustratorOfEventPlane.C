/*
 * Macro to illustrate the event plane direction
 */
void illustratorOfEventPlane(){

	// Define the partices used to determine the event plane
	vector<double> particlePhiVector = {TMath::Pi()/2.0+0.5};

  // Calculate the second order event plane direction
  double eventPlaneQx = 0;
  double eventPlaneQy = 0;
  for(double particlePhi : particlePhiVector){
    eventPlaneQx += TMath::Cos(2.0*particlePhi);
    eventPlaneQy += TMath::Sin(2.0*particlePhi);
  }
  double eventPlaneAngle = 0.5 * TMath::ATan2(eventPlaneQy, eventPlaneQx);

  // Illustrate the event plane angle with a simple drawing
  cout << "Event plane angle is: " << eventPlaneAngle << endl;

  // Create a histogram for jet-event plane correlation
  // DeltaPhi in [-pi/2,3pi/2]
  const Double_t minDeltaPhiJetEventPlane = -TMath::Pi()/2.0;    // Minimum deltaPhi for jet-event plane correlations
  const Double_t maxDeltaPhiJetEventPlane = 3.0*TMath::Pi()/2.0; // Maximum deltaPhi for jet-event plane correlations
  const Int_t nDeltaPhiBinsJetEventPlane = 200;                  // Number of deltaPhi bins for jet-event plane correlations
  TH1D* hJetEventPlane = new TH1D("jetEventPlane", "jetEventPlane", nDeltaPhiBinsJetEventPlane, minDeltaPhiJetEventPlane, maxDeltaPhiJetEventPlane);
  hJetEventPlane->Sumw2();

  // Simulate jet-event plane correlation with a gap in jet acceptance
  TRandom3* myRandom = new TRandom3();
  myRandom->SetSeed(0);
  TF1* eventPlaneCorrelation = new TF1("correlation", "0.2 * cos(2 * x - [0]) + 1", -TMath::Pi(), TMath::Pi());
  eventPlaneCorrelation->SetParameter(0,0);

  // Need to introduce correlation here and then see if acceptance makes it warped
  const int nEvents = 300000;
  double jetPhi;
  double eventPlanePhi;
  double jetEventPlaneDeltaPhi;
  for(int iEvent = 0; iEvent < nEvents; iEvent++){
    if(iEvent % 100000 == 0) cout << "Doing event " << iEvent << endl;

    // Get random phi for the jet and skip events whete this is between 0 and 1
    eventPlanePhi = myRandom->Uniform(-TMath::Pi(), TMath::Pi());
    eventPlaneCorrelation->SetParameter(0, 2*eventPlanePhi);
    jetPhi = eventPlaneCorrelation->GetRandom();
    if(jetPhi > 0){
      if(jetPhi < 1){
        continue;
      }
    }

    // Calculate the deltaPhi between jet and event plane
    jetEventPlaneDeltaPhi = jetPhi - eventPlanePhi;
    while(jetEventPlaneDeltaPhi > (1.5*TMath::Pi())){jetEventPlaneDeltaPhi += -2*TMath::Pi();}
    while(jetEventPlaneDeltaPhi < (-0.5*TMath::Pi())){jetEventPlaneDeltaPhi += 2*TMath::Pi();}

    hJetEventPlane->Fill(jetEventPlaneDeltaPhi);

  }

  // Draw the correlation histogram between jet and event plane
  hJetEventPlane->Draw();
  

}
