#!/bin/bash

# Set Paths
source paths.mCP.sh

# Set arguments
charge=$1
mass=$2
nEv=$3
proc=$4
configFile=$5
nCores=$6

#Pick correct distribution files
cnew="$(echo $charge | sed 's/0.//')"
if [ "$cnew" -lt "0100" ]; then
        sourcecharge=0.001
elif [ "$cnew" -lt "1000" ]; then
        sourcecharge=0.01
elif [ "$cnew" -ge "1000" ]; then
        sourcecharge=0.1
fi

#Define names
outputname="$proc"."$mass"GeV."$charge"Q."$configFile"Config
sourcename="$proc"/"$mass"/"$sourcecharge"/hit_4_vecs.txt
JOB=$SCRATCH/Scratchy."$outputname"
SRC=$JOB/geant4/src
CONFIG=$JOB/geant4/config
echo "$outputname" "$nEv"

# Prepare configuration
rm -r $JOB
mkdir $JOB
cp -r $GEANT $JOB
cd $JOB
cmake -DGeant4_DIR=$G4DIR geant4/
cd $JOB
sed -i'.bak' -e 's/ElectricCharge =.*/ElectricCharge = '"$charge"'/g' $CONFIG/particles.ini
sed -i'.bak' -e 's/.*fMonopoleMass =.*/fMonopoleMass = '"$mass"'*GeV;/g' $CONFIG/particles.ini
sed -i'.bak' -e 's~FileName.*~FileName = '$sourcename'~g' $CONFIG/particles.ini
sed -i'.bak' -e 's~PathName.*~PathName = '$DATA'/~g' $CONFIG/particles.ini
sed -i'.bak' -e 's/.*beamOn.*/\/run\/beamOn '"$nEv"'/g' $CONFIG/mcp.mac

#Run Program
make -j $nCores MilliQ
./MilliQ $CONFIG/mcp.mac $CONFIG/"$configFile".ini

# Collect Output
cp MilliQ.root $RESULTS/"$outputname".root
echo $mass $charge $nEv >> $RESULTS/NEventsInitial."$proc"."$configFile"Config.dat
rm -r $JOB
