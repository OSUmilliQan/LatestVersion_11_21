#include "globals.hh"
#include "MilliQAnalysis.hh"
#include "G4SystemOfUnits.hh"
#include <math.h>
#include <limits>
#include <numeric>
#include <algorithm>
#include "MilliQWaveform.hh"
#include <tgmath.h>

MilliQAnalysis::MilliQAnalysis(std::vector< std::vector<G4double> > ppmtTime, std::vector< std::vector<G4double> > pscintTime,std::vector< std::vector<G4double> > pscintEn, const boost::property_tree::ptree pt) : fIsActive(false), fPTree(pt), fVerbose(false) {

  fpmtTime = ppmtTime;
  fscintTime = pscintTime;
  fscintEn = pscintEn;



  CAENID();
  OnlineTrigger();
  Waveform();
  if(IsActive() ){
	  CalculateUsefulInfo();
  }

  NearestN();

}

// This class assumes that a trigger Module is a 2 x 2 x Nstack object
// and that groups are neighbouring PMTs along "z"

void MilliQAnalysis::CAENID() {

	  try {
	    boost::property_tree::ini_parser::read_ini(fPTree.get<std::string>("Configuration.GeometryConfigFile"), fGeometryPTree);
	    boost::property_tree::ini_parser::read_ini(fPTree.get<std::string>("Configuration.PMTConfigFile"), fPMTPTree);
	    boost::property_tree::ini_parser::read_ini(fPTree.get<std::string>("Configuration.ParticleConfigFile"), fParticlePTree);
	  }
	  catch(boost::property_tree::ptree_error &e) {
	    G4ExceptionDescription msg;
	    msg << G4endl << "Configuration file " << e.what() << G4endl;
	    G4Exception("MilliQAnalysis::ReadConfiguration()", "MilliQAnalysis::ConfigFileReadError", FatalException, msg);
	  }


	// Calculates the CAEN IDs for which we need to check if there has been activation
	  fNblock = G4ThreeVector(fGeometryPTree.get<G4int>("DetectorGeometry.NBlocks_X"),
	                          fGeometryPTree.get<G4int>("DetectorGeometry.NBlocks_Y"),
	                          fGeometryPTree.get<G4int>("DetectorGeometry.NBlocks_Z"));
	  fNstack = fGeometryPTree.get<G4int>("DetectorGeometry.NStacks");

	  fDetectorSize = fNblock[0]*fNblock[1]*fNblock[2]*fNstack;


	  fcoincidenceThreshold = fParticlePTree.get<G4double>("ParticleProperties.coincidenceThreshold")*ns;

	  G4int Nblocks = fNblock[0]*fNblock[1]*fNblock[2];


	  // Suppose we have a 1x2x3 (x,y,z) layer, the numbering goes as:
	  // 	+y
	  // 	^
	  // 	| 3 4 5
	  // 	| 0 1 2
	  // 	------> +z
	  //   /
	  //  /
	  // \/
	  // -x
	  // The rest of the detector goes in the +x direction

	  // Collects IDs of each module
	  for(G4int iy = 0; iy<fNblock[1]/2; iy++){
		  for(G4int iz = 0; iz<fNblock[2]/2; iz++){

			  std::vector<G4int> SubModules;
			  for(G4int istack = 0; istack < fNstack; istack++){
				  SubModules.push_back(0+2*iz+iy*2*fNblock[2]+Nblocks*istack);
				  SubModules.push_back(1+2*iz+iy*2*fNblock[2]+Nblocks*istack);
				  SubModules.push_back(0+fNblock[2]+2*iz+iy*2*fNblock[2]+Nblocks*istack);
				  SubModules.push_back(1+fNblock[2]+2*iz+iy*2*fNblock[2]+Nblocks*istack);
			  }
			  fModules.push_back(SubModules);
		  }
	  }
}

void MilliQAnalysis::OnlineTrigger() {

	//Loops over all modules
	for(G4int i = 0; i < fNblock[1]*fNblock[2]/4; i++){

		//Loops over all combinations  of PMTs inside module
		for(G4int j=1; j<4*fNstack; j++){

			G4int jind=fModules[i][j];
			if(fpmtTime[jind].size() == 0)
				continue;

			for(G4int k=0; k<j; k++){

    			G4int kind=fModules[i][k];
    			if(fpmtTime[kind].size() == 0)
    				continue;

    			//Skip over neighbouring PMTs inside a group
    			if ( std::abs(jind - kind) == 1)
    				continue;

    			if(fIsActive) break;

   				// Loops over all hit times in 2 candidate PMTs
   				for(unsigned int a = 0; a < fpmtTime[jind].size(); a++){
   					for(unsigned int b = 0; b < fpmtTime[kind].size(); b++){

   						//Compare all combinations of hit times
   						if(fabs(fpmtTime[jind][a]/ns - fpmtTime[kind][b]/ns) <= fcoincidenceThreshold/ns) {
   							fIsActive = true;
   						}
   					}
   				}
    		}
    	}
		if(fIsActive) break;
    }


}

void MilliQAnalysis::Waveform() {

	G4double R = fPMTPTree.get<G4double>("PMT.Resistance");
	G4double G = fPMTPTree.get<G4double>("PMT.Gain");
	G4double tcurrent = fPMTPTree.get<G4double>("PMT.CurrentTime");
	G4double nBins = 20.;

	TriggerCount.resize(fDetectorSize, 0);
	TimeCount.resize(fDetectorSize, 0);
	EventID.resize(fDetectorSize, 0);
	TDC.resize(fDetectorSize, 0);
	WaveMax.resize(fDetectorSize, 0);
	WaveMin.resize(fDetectorSize, 0);
	ChargeIntegral.resize(fDetectorSize, 0);

	MilliQWaveform* MyWaveform;

	for ( unsigned int i = 0; i < fpmtTime.size(); i++){

		if( fpmtTime[i].size() == 0){
			//Create an empty waveform
			MyWaveform = new MilliQWaveform(
					0.0, // fTriggerCount
					0.0, // fTimeCount
					0, // fEventID
					0, // fTDC
					0.0, // fWaveMax
					0.0, // fWaveMin
					0.0); // fChargeIntegral

			continue;
		}


		std::vector<G4double > voltage((G4int)nBins,0);
		G4double firstTime = fpmtTime[i][0];
		G4double lastTime = fpmtTime[i].back();
		G4double Qtotal = 0;
		G4double deltaTime = 0;
		G4double nindex = 0.;
		G4double range = 0;


		if( fpmtTime[i].size() == 1 ){
			deltaTime = tcurrent; // ***** IS THIS RIGHT? THE TOTAL Charge might be bad
			voltage[0]++;
		}

		if ( fpmtTime[i].size() > 1 ){
			deltaTime = (lastTime - firstTime)/nBins;
			range = lastTime- firstTime;

			for (unsigned int j = 0; j < fpmtTime[i].size(); j++){

				nindex = floor( nBins*(fpmtTime[i][j] - firstTime)/range );
				if(lastTime == fpmtTime[i][j]){
					nindex = nindex - 1.;
				}
				voltage[(G4int)nindex]++;
			}
		}


		for(unsigned int k = 0; k < voltage.size(); k++){
			voltage[k] = voltage[k]*R*G/tcurrent;
			Qtotal += voltage[k]*deltaTime/R;

		}

		auto EventWaveMax = std::max_element(voltage.begin(), voltage.end());
		auto EventWaveMin = std::min_element(voltage.begin(),voltage.end());


		MyWaveform = new MilliQWaveform(
				0.0, // fTriggerCount
				0.0, // fTimeCount
				0, // fEventID
				0, // fTDC
				*EventWaveMax, // fWaveMax
				*EventWaveMin, // fWaveMin
				Qtotal); // fChargeIntegral

		//Debugging
		G4cout<<"Charge "<< Qtotal << " EventWaveMax " << *EventWaveMax << " EventWaveMin "
				<< *EventWaveMin << G4endl;


		// Output Info:
		// Output MilliQCAEN info
		TriggerCount[i] = MyWaveform->GetTriggerCount();
		TimeCount[i] = MyWaveform->GetTimeCount();
		EventID[i] = MyWaveform->GetEventID();
		TDC[i] = MyWaveform->GetTDC();
		WaveMax[i] = MyWaveform->GetWaveMax();
		WaveMin[i] = MyWaveform->GetWaveMin();
		ChargeIntegral[i] = MyWaveform->GetChargeIntegral();

	}

	delete(MyWaveform);
}

void MilliQAnalysis::CalculateUsefulInfo(){

	G4int median;
    G4double FirstHitTime = std::numeric_limits<double>::max() ;
    FirstHitScintillator = std::numeric_limits<int>::max();

    PmtMedianHitTimes.resize(fDetectorSize, 0);
    TotalEnergyDeposit.resize(fDetectorSize, 0);
    NumberPMTHits.resize(fDetectorSize, 0);
    TimeOfFlight.resize(fDetectorSize, 0);

	for (unsigned int j = 0; j < fDetectorSize; j++){


		  if( fscintEn[j].size() > 0  && fscintTime[j].size() > 0){

			  // Figures out total Energy deposit for each Scintillator
		  	  TotalEnergyDeposit[j] = std::accumulate(fscintEn[j].begin(),fscintEn[j].end(),0)/MeV;

		  	  // Figures out Scintillator ID of first Hit
		 	  if( fscintTime[j][0] < FirstHitTime ){
		 		  FirstHitTime = fscintTime[j][0];
			 	  FirstHitScintillator = j;
		 	  }


		  }

	  	  // Record All PMT Hit Times
	  	  for(unsigned int k = 0; k < fpmtTime[j].size();k++) {
		  	  PmtAllHitTimes.push_back( fpmtTime[j][k] /ns );
	  	  }


		  // Figures out which PMTs are Active
		  if( fpmtTime[j].size() > 0 ){

			  ActivePMT.push_back(j);

			  // Figures out Median hit time for each PMT
			  median = (fpmtTime[j].size() - 1)/2;
		  	  PmtMedianHitTimes[j] = fpmtTime[j][median]/ns;

		  	  // Figures out Number of PMT Hits
		  	  NumberPMTHits[j] = fpmtTime[j].size() ;

		  	  // Figures out time of flight
		  	  TimeOfFlight[j] = fpmtTime[j][median]/ns - fscintTime[j][0]/ns;

		  }

	}

}


void MilliQAnalysis::NearestN() {

	/*

  //Figures out for which sequence all 3 pmt and Scint light up
  if(fVerbose) G4cout << "//////// " << G4endl;

  std::vector<G4int> activePMT;
  std::vector< std::vector<G4double> > activePMTTimes (fNstack);

  //The following identifies events in which 3 consecutive layers light up

  bool recordEvent = true;

  // For all stack layers
  for(G4int j = 0; j < fNstack; j++) {

    // For all blocks in that stack
    bool layerHit = false;
    for(G4int i = 0; i < fNblock; i++) {
      if(fpmtTime[i + fNblock*j].size() > 0) {
	       layerHit = true;

	       // store PMT numbers with hits
	       activePMT.push_back(i + fNblock*j);

	       // store the times of this cell's hits in the list of hit times in the stack layer
	       for(unsigned int k = 0; k < fpmtTime[i + fNblock*j].size(); k++) activePMTTimes[j].push_back(fpmtTime[i + fNblock*j][k]);
      }
    }

    // if any layer has no hits anywhere, don't record this event
    if(!layerHit) {
      recordEvent = false;
      //break;
    }
    else if(fVerbose) G4cout << " Passed Layer " << j << G4endl;
  }

  if(fVerbose) G4cout << "record Event " << recordEvent << G4endl;

  // If all three layers have hits (anywhere)
  if(recordEvent) {

    // loop through the last 2 stacks/layers
    for(G4int j = 1; j < fNstack; j++) {

      for(G4int k = 0; k < j; k++) {

	       bool tFlag = false;

	      // if any two hits, in different layers, are within 15ns of each other, record the event
      	for(unsigned int a = 0; a < activePMTTimes[k].size(); a++) {
	         for(unsigned int b = 0; b < activePMTTimes[j].size(); b++) {

	            if(fabs(activePMTTimes[k][a]/ns - activePMTTimes[j][b]/ns) < 15/ns) {
	               tFlag = true;
	               break;
	            }
	         }

	         if(tFlag) break;
	      }

	      if(!tFlag) {
	         recordEvent = false;
	         break;
	      }

      }

      if(!recordEvent) break;
    }
  }

  // The next part sets all zeroes three times per event for pmtTimes, timeOfFlight, and activeEvent???

  if(recordEvent) {
    if(fVerbose) G4cout << "Recorded Event" << G4endl;
    fIsActive = true;
    //		ComputeTandESummed();
    G4double ensum = 0;
    for(G4int i = 0; i < fNstack; i++) {
      pmtTimes.push_back(ensum);
      timeOfFlight.push_back(ensum);
      activeEvent.push_back(0);
    }
    totalEdep.push_back(ensum);
  }

  */


  //The following identifies events in which 3 consecutive PMTs light up
  /*
    for(G4int i=0; i < fNblock; i++){
      if(fpmtTime[i].size()>0){
        activePMT.push_back(i);
      }
    }

    for(unsigned int i=0; i < activePMT.size(); i++) {
      bool recordEvent = true;
      for(G4int j=1; j < fNstack; j++){
        if(fpmtTime[activePMT[i]+fNblock*j].size()==0){
          recordEvent = false;
          break;
        }
    }
    if(recordEvent == true){

      for(G4int j=1; j < fNstack; j++){
        for(G4int k = 0; k < j; k++){

          bool tFlag = false;

          for(unsigned int a = 0; a < fpmtTime[activePMT[i]+fNblock*k].size(); a++){
            for(unsigned int b = 0; b < fpmtTime[activePMT[i]+fNblock*j].size(); b++){

              if(fabs(fpmtTime[activePMT[i]+fNblock*k][a]/ns-fpmtTime[activePMT[i]+fNblock*j][b]/ns)<15/ns){
                tFlag = true;
                break;
              }
            }

            if(tFlag == true) break;
         }
         if(tFlag == false){
          recordEvent = false;
          break;
         }

       }

    if(recordEvent == false) break;

    }
  }

    if(recordEvent == true) activeEvent.push_back( activePMT[i] );

  }
    //If passes this, means that the same 3 scint and pmt light up!
    if(activeEvent.size()==1){
      G4cout<<"Recorded Event"<<G4endl;
      for(G4int j = 1; j < fNstack; j++)
        activeEvent.push_back(activeEvent[0]+j*fNblock);
        fIsActive = true;
        ComputeTandE();
    }
  */
}
