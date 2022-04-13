#include "../common/function.C"


//-----------------CMS Data----------------------------------------
// const int     nPots_CMS2022 = 2;
// const double  Rap_CMS2022[nPots_CMS2022]        	= { -2.2,  -1.8 };
// const double  RapErrLow_CMS2022[nPots_CMS2022]  	= { 0.20,  0.20 };
// const double  RapErrHig_CMS2022[nPots_CMS2022]  	= { 0.20,  0.20 };

const int     nPots_CMS2022 = 3;
const double  Rap_CMS2022[nPots_CMS2022]        	= { -2.25,  -2.0,  -1.75 };
//-----------------------------------------------------------------

std::vector<double> energy, flux, fluxTable;
std::vector<double> yTable, NyTable;
std::vector<double>  PhotonEnergy_CMS2022, PhotonFlux_CMS2022, Raps_CMS2022;
const int NnCases = 6;
const int nstep = 100;
const double Emin = 1.000000e-05, Emax = 4.777493e+02;
const double JpsiMass = 3.096916;
const TString nCasesName[NnCases] = {"AnAn", "OnOn", "0nXn", "Xn0n", "OnXnSum", "XnXn"};
TString Name, Case;

double Jpsi_rap2EgammCM(double rap);
double interpolateFlux(const double Egamma);
void plotFlux();
double Derivative(const double x);
void nk2Ny();

void getPhotonFlux()
{
	TString inFileDir 		  = "flux/Flux_XnXn.txt";
	Name = TString(inFileDir);
	Name.ReplaceAll("flux/","");
	Name.ReplaceAll(".txt", "");
	Case = TString(Name);
	Case.ReplaceAll("Flux_","");
	cout<<Name<<endl;

	ifstream myfile(inFileDir);

	//Calculate photon energy from the rap
	for (int i = 0; i < nPots_CMS2022; ++i)
	{
		PhotonEnergy_CMS2022.push_back( JpsiMass/2 * exp(Rap_CMS2022[i]) );
		PhotonEnergy_CMS2022.push_back( JpsiMass/2 * exp(-Rap_CMS2022[i]) );
		Raps_CMS2022		.push_back(Rap_CMS2022[i]);
		Raps_CMS2022		.push_back(-Rap_CMS2022[i]);
	}

	if (myfile.is_open())
	{
		std::string line;
	    // Read one line at a time into the variable line:
	    while(std::getline(myfile, line))
	    {
	        std::vector<double>   	lineData;
	        std::stringstream  		lineStream(line);

	        double value;
	        // Read an integer at a time from the line
	        while(lineStream >> value)
	        {
	            // Add the integers from a line to a 1D array (vector)
	            lineData.push_back(value);
	        }
	        cout<<line<<endl;
	        // When all the integers have been read, add the 1D array
	        energy		.push_back(lineData[0]);
	        flux 		.push_back(lineData[1]);
	        fluxTable	.push_back(lineData[2]);
	    }
	}
	else cout << "ERROR!!! Unable to open file!!!";

	plotFlux();
}

void plotFlux()
{	
	nk2Ny();
	TCanvas *c = new TCanvas();
	c->SetLogy();
	TGraph* gr = new TGraph(yTable.size(), & yTable[0], & NyTable[0]);

	// TF1 *f = new TF1("f","pol9",0,100);
	// gr->Fit("f","R");

	TGraph* points = new TGraph();
	points->SetMarkerColor(kRed);

	for (int i = 0; i < PhotonEnergy_CMS2022.size(); ++i)
	{
		PhotonFlux_CMS2022.push_back(	PhotonEnergy_CMS2022[i] * interpolateFlux(PhotonEnergy_CMS2022[i])	);
		points->SetPoint(i,Raps_CMS2022[i],PhotonFlux_CMS2022[i]);
	}

	auto mg = new TMultiGraph();
	mg->Add(gr);
	mg->Add(points);
	

	mg->SetTitle(";y;dN/dy");
	mg->Draw("AP");
	mg->GetXaxis()->SetRangeUser(-4,4);
	mg->GetYaxis()->SetRangeUser(1e-2,1e4);
	drawLatex(0.3, 0.85, "UPC Pb+Pb #sqrt{s_{NN}} = 5.02 TeV (" + Case +")",      42,       0.05,      1);
	// drawLatex(0.3, 0.80, Form("Emin = %.0e GeV to Emax = %.f GeV (CM frame)", Emin, Emax),      42,       0.04,      1);
	for (int i = 0; i < PhotonEnergy_CMS2022.size(); ++i)
	{
		drawLatex(0.15, 0.35-i*0.04, Form("y = %.2f, E = %.2f GeV, dN/dy = %.3f ",Raps_CMS2022[i],PhotonEnergy_CMS2022[i],PhotonFlux_CMS2022[i]),      42,       0.04,      1);
	}

	c->SaveAs( "out4Flux/" + Name + "_dNdy.png" );
	c->SaveAs( "out4Flux/" + Name + "_dNdy.pdf" );	
}

double Jpsi_rap2EgammCM(double rap)
{
	double PbMass = 0.938 * 207.2;
	double Ebeam = 2510*207.2;
	double EgammLab = JpsiMass/2.0 * exp(rap);
	double M2 = pow((EgammLab + Ebeam),2) + pow(( EgammLab - sqrt(Ebeam*Ebeam-PbMass*PbMass) ), 2);
	double EgammaCM = ( M2-PbMass*PbMass ) / ( 2*sqrt(M2) );

	return EgammaCM;
}

double interpolateFlux(const double Egamma)
{

	double lnEmin=log(Emin);
	double lnEmax=log(Emax);
	double lEgamma=log(Egamma);
	double dlnE=(lnEmax-lnEmin)/nstep; 

	//       >> Egamma between Ilt and Ilt+1
	int Ilt = int((lEgamma-lnEmin)/dlnE);
	//       >> ln(Egamma) for first point 
	double lnElt = lnEmin + Ilt*dlnE; 
	//       >> Interpolate
	double flux_r = fluxTable[Ilt] + ((lEgamma-lnElt)/dlnE)*(fluxTable[Ilt+1]- fluxTable[Ilt]);
	flux_r = flux_r/Egamma;

	return flux_r;
}

double Derivative(const double x)
{
	const double dx = 1e-6;
	// return - (interpolateFlux(x + dx) - interpolateFlux(x)) / dx;
	return -(-interpolateFlux(x + 2*dx) + 8*interpolateFlux(x + dx) - 8*interpolateFlux(x - dx) + interpolateFlux(x - 2*dx))/(12*dx);
}

void nk2Ny()
{
	cout<<flux.size()<<endl;
	for (int i = 0; i < flux.size(); ++i)
	{
		yTable 		.push_back( log( 2 / JpsiMass * energy[i] )	);
		NyTable 	.push_back(	energy[i] * interpolateFlux(energy[i])	);
	}
}


	