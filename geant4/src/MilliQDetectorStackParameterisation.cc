
#include "MilliQDetectorStackParameterisation.hh"

MilliQDetectorStackParameterisation::MilliQDetectorStackParameterisation(G4int  pN,
                                                                         G4ThreeVector  pStackDimensions,
                                                                         G4ThreeVector  pDistributionUnitVector,
																		 G4ThreeVector	pOffset,
                                                                         G4double       pStartDepth,
                                                                         G4double       pEndDepth)
{
    fStartCenterDepth = pStartDepth+(pStackDimensions.x()/2.);
    fN = pN;
    
    fDepthGap = ( pEndDepth-pStartDepth-pStackDimensions.x() )/ (pN==1 ? pN : pN-1);
    
    fStackDimensions = pStackDimensions;
    fDistributionUnitVector = pDistributionUnitVector;
    fOffset = pOffset;
}


MilliQDetectorStackParameterisation::~MilliQDetectorStackParameterisation(){}


void MilliQDetectorStackParameterisation::ComputeTransformation(const G4int pId,
                                                                G4VPhysicalVolume*  pPV) const
{
    G4ThreeVector translation = (fStartCenterDepth + pId*fDepthGap)*fDistributionUnitVector;
    if(pId == 1)
    	translation += fOffset;
    if(pId == 2)
    	translation += -fOffset;
    pPV->SetTranslation(translation);
    //pPV->SetRotation( new G4RotationMatrix( fDistributionUnitVector.getPhi() , fDistributionUnitVector.getTheta() , fDistributionUnitVector.getRho() ) );
}


G4ThreeVector MilliQDetectorStackParameterisation::GetStackDimensions()
{
    G4ThreeVector stackDimensions;
    stackDimensions.setX(fStackDimensions.x());
    stackDimensions.setY(fStackDimensions.y());
    stackDimensions.setZ((fStackDimensions.z()*fN) + (fDepthGap*(fN-1)));
    
    return stackDimensions;
}
