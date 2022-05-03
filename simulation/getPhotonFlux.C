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
std::vector<double> biter_0n0n, PofPhotonB_0n0n, PofHadronB_0n0n, PofB_0n0n;
std::vector<double> biter_0nXnSum, PofPhotonB_0nXnSum, PofHadronB_0nXnSum, PofB_0nXnSum;
std::vector<double> biter_XnXn, PofPhotonB_XnXn, PofHadronB_XnXn, PofB_XnXn;
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
void plotPofB();
double Derivative(const double x);
void nk2Ny();

void getPhotonFlux()
{
	TString inFileDir 		  = "flux/Flux_AnAn.txt";
	TString inFileDir1 		  = "flux/PofB_0n0n.txt";
	TString inFileDir2 		  = "flux/PofB_0nXnSum.txt";
	TString inFileDir3 		  = "flux/PofB_XnXn.txt";
	Name = TString(inFileDir);
	Name.ReplaceAll("flux/","");
	Name.ReplaceAll(".txt", "");
	Case = TString(Name);
	Case.ReplaceAll("Flux_","");
	cout<<Name<<endl;

	ifstream myfile(inFileDir);
	ifstream myfile1(inFileDir1);
	ifstream myfile2(inFileDir2);
	ifstream myfile3(inFileDir3);

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
	else cout << "ERROR!!! Unable to open Flux file!!!";


	if (myfile1.is_open())
	{
		std::string line;
	    // Read one line at a time into the variable line:
	    while(std::getline(myfile1, line))
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
	        biter_0n0n			.push_back(lineData[0]);
	        PofPhotonB_0n0n 	.push_back(lineData[1]);
	        PofHadronB_0n0n 	.push_back(lineData[2]);
	        PofB_0n0n 	 		.push_back(lineData[3]);
	    }
	}
	else cout << "ERROR!!! Unable to open PofB file!!!";


	if (myfile2.is_open())
	{
		std::string line;
	    // Read one line at a time into the variable line:
	    while(std::getline(myfile2, line))
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
	        biter_0nXnSum			.push_back(lineData[0]);
	        PofPhotonB_0nXnSum 		.push_back(lineData[1]);
	        PofHadronB_0nXnSum 		.push_back(lineData[2]);
	        PofB_0nXnSum 	 		.push_back(lineData[3]);
	    }
	}
	else cout << "ERROR!!! Unable to open PofB file!!!";

	if (myfile3.is_open())
	{
		std::string line;
	    // Read one line at a time into the variable line:
	    while(std::getline(myfile3, line))
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
	        biter_XnXn			.push_back(lineData[0]);
	        PofPhotonB_XnXn 	.push_back(lineData[1]);
	        PofHadronB_XnXn 	.push_back(lineData[2]);
	        PofB_XnXn 	 		.push_back(lineData[3]);
	    }
	}
	else cout << "ERROR!!! Unable to open PofB file!!!";	
		/* code */

	plotFlux();
	plotPofB();
}

void plotFlux()
{	
	nk2Ny();
	TCanvas *c = new TCanvas();
	c->SetLogy();

	TH2D* htem2d = new TH2D("htem2d", ";y;dN/dy", 10,-4,4, 10, 1e-2, 1e3);

	TGraph* gr = new TGraph(yTable.size(), & yTable[0], & NyTable[0]);
	TGraph* points = new TGraph();

	for (int i = 0; i < PhotonEnergy_CMS2022.size(); ++i)
	{
		PhotonFlux_CMS2022.push_back(	PhotonEnergy_CMS2022[i] * interpolateFlux(PhotonEnergy_CMS2022[i])	);
		points->SetPoint(i,Raps_CMS2022[i],PhotonFlux_CMS2022[i]);
	}

	gr->SetLineColor(kBlack);
	gr->SetLineWidth(3);
	points->SetMarkerStyle(kFullCircle);
	points->SetMarkerColor(kRed);

	htem2d->Draw();
	gr->Draw("SAME");
	// points->Draw("SAME P");

	drawLatex(0.3, 0.85, "UPC Pb+Pb #sqrt{s_{NN}} = 5.02 TeV (" + Case +")",      42,       0.05,      1);
	// drawLatex(0.3, 0.80, Form("Emin = %.0e GeV to Emax = %.f GeV (CM frame)", Emin, Emax),      42,       0.04,      1);
	for (int i = 0; i < PhotonEnergy_CMS2022.size(); ++i)
	{
		// drawLatex(0.15, 0.35-i*0.04, Form("y = %.2f, E = %.2f GeV, dN/dy = %.3f ",Raps_CMS2022[i],PhotonEnergy_CMS2022[i],PhotonFlux_CMS2022[i]),      42,       0.04,      1);
	}

	c->SaveAs( "out4Flux/" + Name + "_dNdy.png" );
	c->SaveAs( "out4Flux/" + Name + "_dNdy.pdf" );
	delete c;
	delete gr;
	delete points;
	delete htem2d;
}

void plotPofB()
{
	TCanvas *c = new TCanvas();
	c->SetLogx();

	TH2D* FluxNeuConfig = new TH2D("FluxNeuConfig", ";b (fm);P_{fn}(b)", 10,6,1e3, 10, 0, 1.1);

	FluxNeuConfig->GetYaxis()->SetTitleSize(0.06);
	FluxNeuConfig->GetYaxis()->SetTitleOffset(0.85);
	FluxNeuConfig->GetYaxis()->SetLabelSize(0.04);
	FluxNeuConfig->GetXaxis()->SetTitleSize(0.06);
	FluxNeuConfig->GetXaxis()->SetTitleOffset(0.95);
	FluxNeuConfig->GetXaxis()->SetLabelSize(0.04);
	FluxNeuConfig->SetTickLength(0.04);
	FluxNeuConfig->Draw();

	TGraph* gr1 = new TGraph(biter_0n0n.size(), & biter_0n0n[0], & PofPhotonB_0n0n[0]);
	TGraph* gr2 = new TGraph(biter_0nXnSum.size(), & biter_0nXnSum[0], & PofPhotonB_0nXnSum[0]);
	TGraph* gr3 = new TGraph(biter_XnXn.size(), & biter_XnXn[0], & PofPhotonB_XnXn[0]);
	gr1->SetTitle("0n0n");
	gr1->SetLineColor(kBlue);
	gr1->SetMarkerColor(kBlue);
	gr1->SetLineWidth(3);
	gr2->SetTitle("0nXn");
	gr2->SetLineColor(kBlack);
	gr2->SetMarkerColor(kBlack);
	gr2->SetLineWidth(3);
	gr3->SetTitle("XnXn");
	gr3->SetLineColor(kRed);
	gr3->SetMarkerColor(kRed);
	gr3->SetLineWidth(3);

	gr1->Draw("same");
	gr2->Draw("same");
	gr3->Draw("same");
	// gPad->BuildLegend();

	auto legend = new TLegend(0.55, 0.25, 0.88, 0.50);
	legend->AddEntry(gr1);
	legend->AddEntry(gr2);
	legend->AddEntry(gr3);
	legend->Draw();

	c->SaveAs( "out4Flux/PofB.png" );
	c->SaveAs( "out4Flux/PofB.pdf" );

	delete c;
	delete FluxNeuConfig;
	delete gr1; delete gr2; delete gr3;
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


	