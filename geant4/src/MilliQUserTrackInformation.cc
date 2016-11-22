
#include "MilliQUserTrackInformation.hh"


MilliQUserTrackInformation::MilliQUserTrackInformation()
  : fStatus(active),fReflections(0),fForcedraw(false) {}


MilliQUserTrackInformation::~MilliQUserTrackInformation() {}


void MilliQUserTrackInformation::AddTrackStatusFlag(int s)
{
  if(s&active) //track is now active
    fStatus&=~inactive; //remove any flags indicating it is inactive
  else if(s&inactive) //track is now inactive
    fStatus&=~active; //remove any flags indicating it is active
  fStatus|=s; //add new flags
}
