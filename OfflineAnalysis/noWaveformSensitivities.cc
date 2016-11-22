//Standard c++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include "sstream"
#include <stdio.h>
#include <unistd.h>
#include <string.h>

//Library for big numbers
#include <gmp.h>

//ROOT Libraries
#include <TRandom1.h>
#include <TH1F.h>
#include <TApplication.h>
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TStyle.h"
#include "Math/WrappedTF1.h"
#include "Math/RootFinder.h"
#include "Math/BrentRootFinder.h"

using namespace std;

//Set constants in GeV
const double massMuon = 0.105658;

const double widthJPsi = 5.93 / 100.;
const double massJPsi = 3.097;

const double widthY1S = 2.48 / 100.;
const double massY1S = 9.460;

const double widthY2S = 1.93 / 100.;
const double massY2S = 10.023;

const double widthY3S = 2.18 / 100.;
const double massY3S = 10.355;

//Cross section k factor
const double kFactor = 1.3;


//Cross section correction factor for Onium decays in mCP 
double phaseSpaceRatio(double mmcp, double mOnia, double muBR) {

	double a = mmcp * mmcp / (mOnia * mOnia);
	double b = massMuon * massMuon / (mOnia * mOnia);

	double numerator = TMath::Sqrt(1. - 4 * a) * (1. + 2 * a);
	double denominator = TMath::Sqrt(1. - 4 * b) * (1. + 2 * b);
	return muBR*numerator / denominator;
}


//Reads data from cross section, acceptance and efficiency, and adds them to myVec
void setVector(string fFilename, vector< vector<double>> &myVec){

	//Picks the directory containing README, CC
	//Assumes DataFiles is inside it
	char the_path[256];
	string fPath2 = "/../DataFiles/";
	string fPath1 = getcwd(the_path, 255);
	string fPath = fPath1.append(fPath2).c_str();

	std::ifstream infile;
	infile.open(fPath.append(fFilename).c_str());

	string line;
	double a, b, c, d;
	double mOld = -1;
	vector<double> mySubVec;

	while (std::getline(infile, line)) {

		std::istringstream iss(line);
		iss.clear();
		iss.str(line);

		// a always corresponds to mass, b, c, d depend on file
		// Assumes all 46 mass points for mCP_UFO, 14 for Y1,2,3S and 7 for JPsi
		iss >> a >> b >> c >> d;

		//Sets mOld to the first mass entry
		if(mOld == -1)
			mOld = a;


		//Read N-1 masses
		if(mOld != a){
			mOld = a;
			//Adds mySubVec to myVec
			myVec.push_back(mySubVec);
			mySubVec.clear();
		}

		//Cross Sections:
		//Assumes ordering by mass
		if(fFilename == "XSection.MilliQUFO.14TeV.dat") mySubVec.push_back(b);
		else if(fFilename == "XSection.JPsi.14TeV.dat")	mySubVec.push_back(b*phaseSpaceRatio(a, massJPsi, widthJPsi ));
		else if(fFilename == "XSection.Y1S.14TeV.dat")	mySubVec.push_back(b*phaseSpaceRatio(a, massY1S, widthY1S ));
		else if(fFilename == "XSection.Y2S.14TeV.dat")	mySubVec.push_back(b*phaseSpaceRatio(a, massY2S, widthY2S ));
		else if(fFilename == "XSection.Y3S.14TeV.dat")	mySubVec.push_back(b*phaseSpaceRatio(a, massY3S, widthY3S ));

		//Efficiencies
		//Assumes ordering by mass, and by charge
		if(fFilename.find("Acceptances") != std::string::npos ){
			mySubVec.push_back(b);
			mySubVec.push_back(d/c);
		}

		//Acceptances
		//Assumes ordering by mass, and then by charge, assumes 3 charges 0.001, 0.01, 0.1
		if(fFilename.find("acceptances.geometric.v3") != std::string::npos){
			mySubVec.push_back(c);
		}

	}

	// Read out last mass info
	if( infile.eof() ){
		myVec.push_back(mySubVec);
		mySubVec.clear();
	}


}


//Efficiency formula for 3 scintillator layers
Double_t EfficiencyFunction(Double_t *x, Double_t *par)
{
   Double_t xx =x[0];
   Double_t f1 = TMath::Power(par[0],3);
   Double_t f2 = TMath::Power(xx/par[1],2);
   Double_t f = f1*TMath::Power(1-TMath::Exp(-f2), 3);
   return f;
}

//Performs fit of G4 efficiency data as a function of charge based on EfficiencyFunction
void setEfficiencyFit(vector<vector<double>> &myVec, vector<vector<double>> myData){
	// Fit Efficiencies
	for(unsigned int i = 0; i < myData.size(); i++){

		Double_t x[myData[i].size() / 2];
		Double_t y[myData[i].size() / 2];
		vector<double> mySubVec;

		for(unsigned j = 0; j < myData[i].size() / 2; j++){
			x[j] = myData[i][2*j];
			y[j] = myData[i][2*j+1];
		}

		TGraph *gr1 = new TGraph(sizeof(x)/sizeof(x[0]),x,y);
		TF1 *fit = new TF1("fit",EfficiencyFunction,0,2,2);

		//Define starting point for fit
		fit->SetParameter(0,0.95);
		fit->SetParameter(1,0.002);
		gr1->Fit("fit","q");

		Double_t p0 = fit->GetParameter(0);
		Double_t p1 = fit->GetParameter(1);

//		cout<<" p0 "<< p0 << " p1 "<<p1<<endl;

		mySubVec.push_back(p0);
		mySubVec.push_back(p1);

		myVec.push_back(mySubVec);
		mySubVec.clear();

	}
}

//Interpolates acceptance as a function of charge
//Assumes 3 charges in acceptance data: 0.001, 0.01, 0.1
Double_t AcceptanceFunction(Double_t *x, Double_t *par)
{
	Double_t xv =x[0];
	Double_t f =
			(xv >= 0.000 && xv < 0.001)*(0.000000000001)+
			(xv >= 0.001 && xv < 0.010)*(par[0]+(xv-0.001)*(par[1]-par[0])/(0.01-0.001)) +
			(xv >= 0.010 && xv < 0.100)*(par[1]+(xv-0.01)*(par[2]-par[1])/(0.1-0.01)) +
			(xv >= 0.100)*(par[2]);
	return f;
}


// Outputs number of signal events minus number of events consistent with X% C.L. from ratio likelihood analysis as a function of charge
// given efficiency, acceptance, cross section, luminosity
Double_t NEventsFunction(Double_t *x, Double_t *par){

	Double_t xv = x[0];
	vector<Double_t> myXSection;
	vector<Double_t> myEffFit;
	vector<Double_t> myAccept;

	//Set NSig for 0
	Double_t f0 = par[0];
	// Set Lumi for 31
	Double_t f1 = par[31];
	//Set Cross-Section kFactor
	Double_t f2 = kFactor;
	//Set epsilon^2
	Double_t f3 = TMath::Power(xv,2);
	int NProcess = par[32];


	//Set XSections for 1 - 5
	for(unsigned int i = 0; i < NProcess; i++){
			myXSection.push_back( par[1+i] ) ;
	}
	// Set EffFit for 6 - 15
	for(unsigned int i = 0; i < NProcess; i++){
			myEffFit.push_back( par[ 6 + 2*i] );
			myEffFit.push_back( par[6 + 2*i + 1] );
	}
	// Set Accept for 16 - 30
	for(unsigned int i = 0; i < NProcess; i++){
			myAccept.push_back( par[16 + 3*i] );
			myAccept.push_back( par[16 + 3*i + 1] );
			myAccept.push_back( par[16 + 3*i + 2] );
	}


	Double_t output = 0;
	for(unsigned int i = 0; i < NProcess; i++){
		Double_t f4 = myXSection[i];
		Double_t dEff[2] = {myEffFit[2*i],myEffFit[2*i+1]};
		Double_t dAccept[3] = {myAccept[3*i], myAccept[3*i+1], myAccept[3*i+2] };
		Double_t f5 = EfficiencyFunction(x, dEff);
		Double_t f6 = AcceptanceFunction(x, dAccept);

		output += f4*f5*f6;

	}

	return f1*f2*f3*output-f0;

}

//Does factorial using the big number library
void factorial(unsigned long long n, mpf_t &result) {
    mpf_set_ui(result, 1);

    while (n > 1) {
        mpf_mul_ui(result, result, n);
        n = n-1;
    }
}


//Applies likelihood ratio analysis to determine number of signal events consistent 
// with X% C.L., given number of background events 
// Uses big library because c++ can't handle 150 factorial.
Double_t LikelihoodRatio(Double_t *x, Double_t *par){

	Double_t xx = x[0];
	Double_t Nbg = par[0];
	Double_t alpha = par[1];
	Double_t f1 = TMath::Exp(-xx);


	mpf_t a, b, f2, fhelper, f3 ;

	mpf_init (f2);
	mpf_init (f3);
	mpf_init (a);
	mpf_init (b);
	mpf_init (fhelper);

	mpf_set_d (f2, 0);
	mpf_set_d (f3, 0);


	for(int i = 0; i <= Nbg; i++){
		//Set Things
		mpf_set_d (a, xx+Nbg);
		factorial(i, b);
		mpf_set_d (fhelper,0);
		//Do algebra
		mpf_add(fhelper,fhelper,a);
		mpf_pow_ui(fhelper, fhelper, i);
		mpf_div(fhelper, fhelper, b);
		mpf_add(f2, f2, fhelper);
	}

	for(int i = 0; i <= Nbg; i++){
		//Set things
		mpf_set_d (a, Nbg);
		factorial(i, b);
		mpf_set_d (fhelper,0);
		//Do algebra
		mpf_add(fhelper,fhelper,a);
		mpf_pow_ui(fhelper, fhelper, i);
		mpf_div(fhelper, fhelper, b);
		mpf_add(f3, f3, fhelper);
	}

	return alpha - f1*mpf_get_d(f2)/mpf_get_d(f3);

}

//Implements LikelihoodRatio function given C.L. criteria and number bkg events
void setNSig(Double_t myNBkg, Double_t myAlpha, vector<double> &myVec){

	Double_t lowBound = 0;
	Double_t upBound = 500;

	// Create the function and wrap it
	TF1 *fNSig = new TF1("fNSig", LikelihoodRatio, lowBound, upBound, 2 );
	fNSig->SetParameter(0, myNBkg);
	fNSig->SetParameter(1, myAlpha);

	ROOT::Math::WrappedTF1 wfNSig(*fNSig);
	// Create the Integrator & Set parameters of the method
	ROOT::Math::BrentRootFinder brfNSig;
	brfNSig.SetFunction( wfNSig, lowBound, upBound );
	brfNSig.Solve();
	myVec.push_back(brfNSig.Root());

}

int main(int argc, char* argv[]) {

	//Sets Masses
	const int nMasses = 46;
	vector<double> massListUFO, massListJPsi, massListYS; 

	for(int i = 0; i < nMasses; i++) 
		massListUFO.push_back( 0.1 * TMath::Power(10., TMath::Sqrt(0.2 * i)));

	for(unsigned int i = 0; i < 14; i++) {
		massListYS.push_back(massListUFO[i]);
		if(i < 7) massListJPsi.push_back(massListUFO[i]);
	}




	//Calculates NSig from likelihood ratio
	vector<double> NSig;
	vector<double> Alpha {0.05, 0.0027}; //95% and 3sigma C.L.
	vector<double> NBkg {150, 300};

	setNSig(NBkg[0], Alpha[0], NSig); // 150Bkg, 95%
	setNSig(NBkg[0], Alpha[1], NSig); // 150Bkg, 3sigma
	setNSig(NBkg[1], Alpha[0], NSig); // 300Bkg, 95%
	setNSig(NBkg[1], Alpha[1], NSig); // 300Bkg, 3sigma

	//Set Luminosity to match the NSig vector
	vector<double> Lumi {300*1000, 300*1000, 3000*1000, 3000*1000};


	// Stores Cross Section
	vector<vector <double>> XSectionUFO, XSectionJPsi, XSectionY1S, XSectionY2S, XSectionY3S;
	// Stores Charges and Efficiencies
	vector<vector <double>> EffUFO, EffJPsi, EffY1S, EffY2S, EffY3S;
	// Stores Charges and Acceptances
	vector<vector <double>> AcceptUFO, AcceptJPsi, AcceptY1S, AcceptY2S, AcceptY3S;
	// Stores Efficiency Fit Parameters
	vector<vector<double>> EffFitUFO, EffFitJPsi, EffFitY1S, EffFitY2S, EffFitY3S;


	//Import Madgraph CrossSections
	setVector("XSection.MilliQUFO.14TeV.dat", XSectionUFO);
	setVector("XSection.JPsi.14TeV.dat", XSectionJPsi);
	setVector("XSection.Y1S.14TeV.dat", XSectionY1S);
	setVector("XSection.Y2S.14TeV.dat", XSectionY2S);
	setVector("XSection.Y3S.14TeV.dat", XSectionY3S);

	//Import G4 Efficiencies
	setVector("Acceptances.mCP_UFO.FullSettings.5b5.refl0.98.LG.dat", EffUFO);
	setVector("Acceptances.JPsi.FullSettings.5b5.refl0.98.LG.dat", EffJPsi);
	setVector("Acceptances.Y1S.FullSettings.5b5.refl0.98.LG.dat", EffY1S);
	setVector("Acceptances.Y2S.FullSettings.5b5.refl0.98.LG.dat", EffY2S);
	setVector("Acceptances.Y3S.FullSettings.5b5.refl0.98.LG.dat", EffY3S);

	//Perform G4 Efficiency Fits
	setEfficiencyFit(EffFitUFO, EffUFO);
	setEfficiencyFit(EffFitJPsi, EffJPsi);
	setEfficiencyFit(EffFitY1S, EffY1S);
	setEfficiencyFit(EffFitY2S, EffY2S);
	setEfficiencyFit(EffFitY3S, EffY3S);

	//Import Geometric Acceptances
	setVector("acceptances.geometric.v3.mCP_UFO.dat", AcceptUFO);
	setVector("acceptances.geometric.v3.JPsi.dat", AcceptJPsi);
	setVector("acceptances.geometric.v3.Y1S.dat", AcceptY1S);
	setVector("acceptances.geometric.v3.Y2S.dat", AcceptY2S);
	setVector("acceptances.geometric.v3.Y3S.dat", AcceptY3S);

	Double_t lowRange = 0.001;
	Double_t upRange = 3;

	vector<vector<vector<double>>> XSection { XSectionUFO, XSectionY1S, XSectionY2S, XSectionY3S, XSectionJPsi};
	vector<vector<vector<double>>> EffFit { EffFitUFO, EffFitY1S, EffFitY2S, EffFitY3S, EffFitJPsi};
	vector<vector<vector<double>>> Accept { AcceptUFO, AcceptY1S, AcceptY2S, AcceptY3S, AcceptJPsi};

	vector<vector<double>> Sensitivities;
	vector<double> sensitivityPre;

	//Calculate the Sensitivities
	int NProcess = 5; //We can through the masses, this tells us how many processes contribute given a mCP mass
	for(unsigned int j = 0; j < XSectionUFO.size(); j++){

		if( j >= XSectionJPsi.size() )
			NProcess = 4;
		if( j >= XSectionY1S.size() )
			NProcess = 1;


		// Create the function and wrap it
		TF1 *fNEv = new TF1("fNEv", NEventsFunction, lowRange, upRange, 33 );

		//Set XSections for 1 - 5
		for(unsigned int i = 0; i < NProcess; i++){
			fNEv->SetParameter(1+i,XSection[i][j][0]);
		}

		// Set EffFit for 6 - 15
		for(unsigned int i = 0; i < NProcess; i++){
			fNEv->SetParameter(6 + 2*i, EffFit[i][j][0]);
			fNEv->SetParameter(6 + 2*i + 1, EffFit[i][j][1]);
		}
		// Set Accept for 16 - 30
		for(unsigned int i = 0; i < NProcess; i++){
			fNEv->SetParameter(16 + 3*i, Accept[i][j][0]);
			fNEv->SetParameter(16 + 3*i + 1, Accept[i][j][1]);
			fNEv->SetParameter(16 + 3*i + 2, Accept[i][j][2]);
		}
		// Set Number of Processes
		fNEv->SetParameter(32, NProcess);
		sensitivityPre.clear();
		for(unsigned int a = 0; a < NSig.size(); a++){

			// Set Lumi for 31
			fNEv->SetParameter(31, Lumi[a]);
			//Set NSig for 0
			fNEv->SetParameter(0,NSig[a]);

			ROOT::Math::WrappedTF1 wf1(*fNEv);
			// Create the Integrator & Set parameters of the method
			ROOT::Math::BrentRootFinder brf;
			brf.SetFunction( wf1, lowRange, upRange );
			brf.Solve();

			sensitivityPre.push_back(brf.Root());

		}
		Sensitivities.push_back(sensitivityPre);

	}

	//Output in a form mathematica can read for doing easy plotting
	ofstream output;
	output.open ("noWaveformSensitivities.dat");
	output<<"{";
	for(unsigned int j = 0; j < NSig.size(); j++){
		output<<"{";
		for(unsigned int i = 0; i < Sensitivities.size(); i ++){
			if(i!= Sensitivities.size() -1)
				output << "{"<<massListUFO[i] << "," << Sensitivities[i][j]<<"},";
			else
				output << "{"<<massListUFO[i] << "," << Sensitivities[i][j]<<"}";
		}
		if(j!=NSig.size()-1)
			output<<"},";
		else
			output<<"}";

	}
	output << "}";
	output.close();


	// ROOT
	// Write the plotting scripts for ROOT, currently easier to simply do them in mathematica


	return 0;

}
