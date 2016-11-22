/*
 * MilliQWaveform.cpp
 *
 *  Created on: Aug 23, 2016
 *      Author: chaffy
 */

#include "../include/MilliQWaveform.hh"

MilliQWaveform::MilliQWaveform(
		G4float fTriggerCount,
		G4float fTimeCount,
		G4int fEventID,
		G4int fTDC,
		G4float fWaveMax,
		G4float fWaveMin,
		G4float fChargeIntegral){

	TriggerCount = fTriggerCount;
	TimeCount = fTimeCount;
	EventID = fEventID;
	TDC = fTDC;
	WaveMax = fWaveMax;
	WaveMin = fWaveMin;
	ChargeIntegral = fChargeIntegral;

}

MilliQWaveform::~MilliQWaveform() {
	// TODO Auto-generated destructor stub
}
