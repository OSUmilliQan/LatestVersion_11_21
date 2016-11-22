
NOTE: 
	All path are given starting from the milliq folder, unless otherwise specified. 
	All information on the simulation and its output is located in milliQanDetector.pdf
	The procedure to calculate limits is specified in Set_Limits_miliQ_notes.pdf
	KnownBugs.txt indicates what currently doesn't work in the simulation

INSTRUCTIONS:

A) Make sure the following is installed:
	1) Install Geant4 and ROOT following instructions in ./recipes/
	2) Download acceptances kinematic distributions from GITHUB. 
	   On OSUT3, copy /home/bfrancis/public/acceptances_v3.zip to someplace in your working directory, and unzip it. 
	   If running on a Windows machine, ask for a web link.
	3) Create Build Directory for the simulation according to the instructions in ./recipes/, type "make".
	   Information for how to create a general MilliQBuild directory in order for all the steps in C) to execute is shown in recipes/compileCommands.txt 

B) Make sure this information is set:
	1) Path to kinematic distributions in A.1: geant4/config/particles.ini  
	2) Mass, Charge, and nanosecond trigger threshold: geant4/config/particles.ini  
	3) Detector specifications: see geant4/config/{default.ini, Geometry, PMT, Scintillator, particles.ini} 
	4) Number of events in a run: geant4/config/mcp.mac
	5) Paths in BASHScript/paths.mCP.sh

C) To run the simulation interactively, via run card or via bash script:
	1) Go in the BuildDirectory described above in A)3) and type "./MilliQ PATH_to_milliq/geant4/config/mcp.mac PATH_to_milliq/geant4/config/default.ini". 
	   The second argument is the location of the mcp.mac file which tells it how may particles to shoot in one run. 
	   The 3rd argument sets all the detector configurations (dimensions, number of PMTs, etc) from the chosen file.
	   If you omit writting the 3rd argument, it chooses geant4/config/default.ini by default.
	2) Alternatively, type "./MilliQ". This will open the interactive Qt visualizer. Make sure Qt is installed.
	3) Alternatively, the simulation can be run via BASHScript/{run.mCP.sh, submit.mcp.exclusions.sh} by typing bash run.mCP.sh
	   In run.mCP.sh, you will set which charges, masses, configuration file, and processes that are scanned over.
	   Make sure BASHScript/paths.mCP.sh is accurate.

D) To analyse the output
	1) To run noWaveformSensitivities.cc, make sure DataFiles/ contains CrossSections, Efficiencies, and Acceptances in the correct format, 
	   for each of the physics processes (JPsi, MilliQ_UFO, Y1S, Y2S, Y3S)
	2) Need to install gmp.hi (GNU Multiple Precision Arithmetic Library)
	3) Go in OfflineAnalysis, all the analysis and plotting files are there, along with the makefiles
	4) Point the Makefile to noWaveformSensitivities.cc, type make to compile and ./noWaveformSensitivities to run the executable.
	5) To plot the output of noWaveformSensitivities according to what we put in the paper, run the Mathematica file. 
	   It just does the plotting (all calculations are in the c++ / ROOT programmes). 


The README/ folder contains additional resources that should be helpful if you want to make modifications to the simulation	
