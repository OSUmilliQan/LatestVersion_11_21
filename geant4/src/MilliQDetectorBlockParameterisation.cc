
#include "MilliQDetectorBlockParameterisation.hh"

MilliQDetectorBlockParameterisation::MilliQDetectorBlockParameterisation(G4ThreeVector  pN,
                                                         G4ThreeVector  pBlockDimensions,
                                                         G4ThreeVector  pGapDimensions)
{
    fN = pN;
    fBlockDimensions = pBlockDimensions;
    fGapDimensions = pGapDimensions;
}


MilliQDetectorBlockParameterisation::~MilliQDetectorBlockParameterisation(){}


void MilliQDetectorBlockParameterisation::ComputeTransformation(const G4int pId,
                                                                G4VPhysicalVolume* pPV) const
{
    G4ThreeVector idDimensions;
    GetIdDimensionsFromId(pId,idDimensions);
    
    // origin = (idDimensions - ((fN-1.)/2)) * (fBlockDimensions+fGapDimensions)
    //
    G4ThreeVector origin;
    G4ThreeVector one(1.,1.,1.);
    G4ThreeVector blockIdPosition =  idDimensions - ((fN-one)/2.);
    G4ThreeVector individualBlockVolume = fBlockDimensions+fGapDimensions;
    origin.setX( blockIdPosition.x() * individualBlockVolume.x() );
    origin.setY( blockIdPosition.y() * individualBlockVolume.y() );
    origin.setZ( blockIdPosition.z() * individualBlockVolume.z() );

    pPV->SetTranslation(origin);
    pPV->SetRotation(0);
}


void MilliQDetectorBlockParameterisation::GetIdDimensionsFromId(G4int           pId,
                                                        G4ThreeVector&  pIdDimensions) const
{
    pIdDimensions.setX( 0. );
    pIdDimensions.setZ( (pId % (G4int)fN.z()) );
    pIdDimensions.setY( (pId-pIdDimensions.z())/(G4int)fN.z() );

}


G4ThreeVector MilliQDetectorBlockParameterisation::GetStackDimensions()
{
    // stackDimensions = fN*(fBlockDimensions+fGapDimensions)
    //
    G4ThreeVector stackDimensions;
    G4ThreeVector individualBlockVolume = fBlockDimensions+fGapDimensions;
    stackDimensions.setX(fN.x() * fBlockDimensions.x());//individualBlockVolume.x());
    stackDimensions.setY(fN.y() * individualBlockVolume.y());
    stackDimensions.setZ(fN.z() * individualBlockVolume.z());
    
    return stackDimensions;
}
