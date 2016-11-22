setenv USERBASE `pwd`
set ARCH=x86_64-slc6-gcc49-opt

source /afs/cern.ch/sw/lcg/contrib/gcc/4.9/${ARCH}/setup.csh 

setenv QTHOME /afs/cern.ch/sw/lcg/external/qt/4.8.4/${ARCH}/
setenv G4BASE /afs/cern.ch/sw/lcg/external/geant4
setenv DAWNHOME /afs/cern.ch/sw/lcg/external/dawn/3_88a/x86_64-slc5-gcc43-opt/
setenv G4DAWNFILE_DEST_DIR ${USERBASE}/DawnFiles/
setenv HEPMC_DIR /afs/cern.ch/sw/lcg/external/HepMC/2.06.08/${ARCH}/
setenv FASTJET_INSTALL /afs/cern.ch/sw/lcg/external/fastjet/3.0.3/${ARCH}/

source $G4BASE/10.2/${ARCH}/GNUmake-setup.csh

#for boost latest version
setenv BOOSTSYS /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/boost/1.57.0-cms
setenv LD_LIBRARY_PATH ${XERCESCROOT}/bin:${XERCESCROOT}/lib:${HEPMC_DIR}/lib:${USERBASE}/userlib/lib:${USERBASE}/analysis/lib:${FASTJET_INSTALL}/lib:${BOOSTSYS}/lib:/afs/cern.ch/sw/lcg/contrib/gcc/4.9/${ARCH}/lib64:/afs/cern.ch/sw/lcg/app/releases/GEANT4/clhep/2.3.1.0/x86_64-slc6-gcc49-opt/lib:/afs/cern.ch/sw/lcg/external/geant4/10.2/x86_64-slc6-gcc49-opt/lib64:${LD_LIBRARY_PATH}

cd /afs/cern.ch/sw/lcg/app/releases/ROOT/6.02.12/${ARCH}/root/
source bin/thisroot.csh

cd -
setenv PATH $DAWNHOME/bin:${PATH}:${FASTJET_INSTALL}/bin:${XERCESCROOT}/bin
