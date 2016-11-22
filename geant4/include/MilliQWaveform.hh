/*
 * MilliQWaveform.h
 *
 *  Created on: Aug 23, 2016
 *      Author: chaffy
 */
#include "globals.hh"


#ifndef SOURCE_DIRECTORY__SRC_MILLIQWAVEFORM_H_
#define SOURCE_DIRECTORY__SRC_MILLIQWAVEFORM_H_

class MilliQWaveform {
public:
	MilliQWaveform(
			G4float fTriggerCount,
			G4float fTimeCount,
			G4int fEventID,
			G4int fTDC,
			G4float fWaveMax,
			G4float fWaveMin,
			G4float fChargeIntegral);

	virtual ~MilliQWaveform();

	G4float GetTriggerCount() { return TriggerCount;};
	G4float GetTimeCount() { return TimeCount; };
	G4int GetEventID() { return EventID; };
	G4int GetTDC() { return TDC; };
	G4float GetWaveMax() { return WaveMax; };
	G4float GetWaveMin() { return WaveMin; };
	G4float GetChargeIntegral() {return ChargeIntegral; };

private:

	G4float TriggerCount;
	G4float TimeCount;
	G4int EventID;
	G4int TDC;
	G4float WaveMax;
	G4float WaveMin;
	G4float ChargeIntegral;

};

#endif /* SOURCE_DIRECTORY__SRC_MILLIQWAVEFORM_H_ */
