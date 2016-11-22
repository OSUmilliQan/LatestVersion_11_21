#!/bin/bash

GCC_ROOT="/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6"
GCC_VERSION="4.9.1"
GCC_REVISION="1"
GCC_CATEGORY="lcg"
export PATH=/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/bin:$PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/lib64:$LD_LIBRARY_PATH

source /cvmfs/sft.cern.ch/lcg/external/ROOT/5.34.23/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

root "$@"
