set USERBASE=`pwd`
set ARCH=x86_64-slc6-gcc49-opt

alias root $USERBASE/myROOT.sh
    
# GCC 4.9.3
set GCC_ROOT="/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6"
set GCC_VERSION="4.9.1"
set GCC_REVISION="1"
set GCC_CATEGORY="lcg"
setenv PATH "/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/bin:$PATH"
setenv LD_LIBRARY_PATH "/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/lib64:$LD_LIBRARY_PATH"

# geant4
source /cvmfs/geant4.cern.ch/geant4/10.2/x86_64-slc6-gcc49-opt/CMake-setup.csh

# QT
setenv QTHOME /cvmfs/sft.cern.ch/lcg/external/qt/4.8.4/x86_64-slc6-gcc48-opt
setenv QTLIBPATH ${QTHOME}/lib

# CMake
set CMAKE_ROOT="/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.3.2/Linux-x86_64"
set CMAKE_VERSION="3.3.2"
set CMAKE_REVISION="1"
set CMAKE_CATEGORY="lcg"
setenv PATH "/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.3.2/Linux-x86_64/bin:$PATH"
setenv LD_LIBRARY_PATH "/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.3.2/Linux-x86_64/lib:$LD_LIBRARY_PATH"

# expat
set EXPAT_ROOT="/cvmfs/sft.cern.ch/lcg/external/expat/2.0.1/x86_64-slc6-gcc48-opt"
set EXPAT_VERSION="2.0.1"
set EXPAT_REVISION="1"
set EXPAT_CATEGORY="lcg"
setenv PATH "/cvmfs/sft.cern.ch/lcg/external/expat/2.0.1/x86_64-slc6-gcc48-opt/bin:$PATH"
setenv LD_LIBRARY_PATH "/cvmfs/sft.cern.ch/lcg/external/expat/2.0.1/x86_64-slc6-gcc48-opt/lib:$LD_LIBRARY_PATH"

# bz2, zlib
source /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/bz2lib/1.0.5-cms/etc/profile.d/init.csh
source /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/zlib/1.2.8-cms/etc/profile.d/init.csh

# ROOT
source /cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/6.02.12/etc/profile.d/init.csh
