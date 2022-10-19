#include "../common/function.C"

std::map<TString, std::vector<double>> PhotonFluxMap	=
{
	{"Energy_Table_AnAn",{}},	{"Raps_Table_AnAn",{}}	,{"dNdk_Table_AnAn",	{}},{"dNdy_Table_AnAn",	{}},
	{"Energy_Table_0n0n",{}},	{"Raps_Table_0n0n",{}}	,{"dNdk_Table_0n0n",	{}},{"dNdy_Table_0n0n",	{}},
	{"Energy_Table_0nXnSum",{}},{"Raps_Table_0nXnSum",{}},{"dNdk_Table_0nXnSum",{}},{"dNdy_Table_0nXnSum",	{}},
	{"Energy_Table_XnXn",{}},{"Raps_Table_XnXn",{}}	,	{"dNdk_Table_XnXn",		{}},{"dNdy_Table_XnXn",	{}},
	{"biter_0n0n",	{}},	{"PofPhotonB_0n0n",	{}},	{"PofHadronB_0n0n",	{}},	{"PofB_0n0n",	{}},
	{"biter_0nXnSum",{}},	{"PofPhotonB_0nXnSum",{}},	{"PofHadronB_0nXnSum",{}},	{"PofB_0nXnSum",{}},
	{"biter_XnXn",	{}},	{"PofPhotonB_XnXn",	{}},	{"PofHadronB_XnXn",	{}},	{"PofB_XnXn",	{}}
};

auto PhotonFluxMapTemp = PhotonFluxMap;

double Interpolate(const double Egamma, TString Case);
void plotFlux(std::map<TString, std::vector<double>> Map, TString Case);
void plotPofB();

void loadPhotonFlux(TString inFileDir	= "flux/",	TString subCase = "")
{
	std::vector<TString> CasesName = {"AnAn", "0n0n", "0nXnSum", "XnXn"};
	const double JpsiMass = 3.096916;

	for (int i = 0; i < CasesName.size(); ++i)
	{
		TString FluxFileName  = Form("Flux_%s%s.txt", CasesName[i].Data(), subCase.Data());
		TString PofBFileName  = Form("PofB_%s%s.txt", CasesName[i].Data(), subCase.Data());
		cout<<"loadPhotonFlux-------->Loading Photon Flux From " + FluxFileName + " in the Dir:"<<inFileDir<<endl;

		ifstream myfile(Form("%s%s", inFileDir.Data(), FluxFileName.Data()));

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
				// cout<<line<<endl;
				// When all the integers have been read, add the 1D array
				PhotonFluxMap.at("Energy_Table_" + CasesName[i]).push_back(lineData[0]);
				PhotonFluxMap.at("dNdk_Table_" + CasesName[i]) 	.push_back(lineData[1]);
				PhotonFluxMap.at("dNdy_Table_" + CasesName[i])	.push_back(lineData[2]);

				PhotonFluxMap.at("Raps_Table_" + CasesName[i])	.push_back( log( 2 / JpsiMass * lineData[0] )	);
			}
		}
		else throw std::runtime_error( "ERROR!!! Unable to open Flux file!!!");

		if (i == 0) continue;
		ifstream myfile1(Form("%s%s", inFileDir.Data(), PofBFileName.Data()));
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
				// cout<<line<<endl;
				// When all the integers have been read, add the 1D array
				PhotonFluxMap.at("biter_" + CasesName[i])		.push_back(lineData[0]);
				PhotonFluxMap.at("PofPhotonB_" + CasesName[i]) 	.push_back(lineData[1]);
				PhotonFluxMap.at("PofHadronB_" + CasesName[i]) 	.push_back(lineData[2]);
				PhotonFluxMap.at("PofB_" + CasesName[i])  		.push_back(lineData[3]);
			}
		}
		else cout << "ERROR!!! Unable to open PofB file!!!";
	}
	cout<<"loadPhotonFlux-------->DONE "<<endl;
}

void plotFlux(std::map<TString, std::vector<double>> Map, TString Case)
{	
	TCanvas *c = new TCanvas();
	c->SetLogy();

	TH2D* htem2d = new TH2D("htem2d_"+Case, ";y;dN/dy", 10,-4,4, 10, 1e-2, 1e3);

	TGraph* gr = new TGraph(PhotonFluxMap.at("Raps_Table_"+Case).size(), & PhotonFluxMap.at("Raps_Table_"+Case)[0], & PhotonFluxMap.at("dNdy_Table_"+Case)[0]);
	TGraph* points = new TGraph();

	for (int i = 0; i < Map.at("Raps").size(); ++i)
	{
		points->SetPoint(i,Map.at("Raps")[i],Map.at("dNdy_"+Case)[i]);
	}
	gr->SetLineColor(kBlack);
	gr->SetLineWidth(3);
	points->SetMarkerStyle(kFullCircle);
	points->SetMarkerColor(kRed);

	htem2d->Draw();
	gr->Draw("SAME l");
	points->Draw("SAME P");

	drawLatex(0.3, 0.85, "UPC Pb+Pb #sqrt{s_{NN}} = 5.02 TeV (" + Case +")",      42,       0.05,      1);
	// drawLatex(0.3, 0.80, Form("Emin = %.0e GeV to Emax = %.f GeV (CM frame)", Emin, Emax),      42,       0.04,      1);
	for (int i = 0; i < Map.at("Raps").size(); ++i)
	{
		drawLatex(0.15, 0.35-i*0.04, Form("y = %.2f, E = %.2f GeV, dN/dy = %.3f ",Map.at("Raps")[i],Map.at("Energy")[i],Map.at("dNdy_"+Case)[i]),      42,       0.04,      1);
	}

	c->SaveAs( "out4Flux/Flux_" + Case + "_dNdy.png" );
	c->SaveAs( "out4Flux/Flux_" + Case + "_dNdy.pdf" );
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

	TGraph* gr1 = new TGraph(PhotonFluxMap.at("biter_0n0n").size(), & PhotonFluxMap.at("biter_0n0n")[0], & PhotonFluxMap.at("PofPhotonB_0n0n")[0]);
	TGraph* gr2 = new TGraph(PhotonFluxMap.at("biter_0nXnSum").size(), & PhotonFluxMap.at("biter_0nXnSum")[0], & PhotonFluxMap.at("PofPhotonB_0nXnSum")[0]);
	TGraph* gr3 = new TGraph(PhotonFluxMap.at("biter_XnXn").size(), & PhotonFluxMap.at("biter_XnXn")[0], & PhotonFluxMap.at("PofPhotonB_XnXn")[0]);
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

	// double min = 6.67*1, max = 10000;
	// TF1 f2("f2",[&](double*x, double*){ return gr2->Eval(x[0]); }, 0, max, 0); 
	// TF1 f3("f3",[&](double*x, double*){ return gr3->Eval(x[0]); }, 0, max, 0); 
	// auto n2 = f2.Integral(min, max, 1e-5);
	// auto n3 = f3.Integral(min, max, 1e-5);
	// cout<<"Ratio 0nXnSum/XnXn:	"<< n2/n3 <<endl;

	// c->SaveAs( "out4Flux/PofB.png" );
	// c->SaveAs( "out4Flux/PofB.pdf" );

	delete c;
	delete FluxNeuConfig;
	delete gr1; delete gr2; delete gr3;
}

double Interpolate(const double Egamma, TString Case)
{
	const int nstep = 100;
	const double Emin = 1.000000e-05, Emax = 4.777493e+02; //Emax = 2.627452e+02 for 2.76TeV

	double lnEmin=log(Emin);
	double lnEmax=log(Emax);
	double lEgamma=log(Egamma);
	double dlnE=(lnEmax-lnEmin)/nstep; 

	//       >> Egamma between Ilt and Ilt+1
	int Ilt = int((lEgamma-lnEmin)/dlnE);
	//       >> ln(Egamma) for first point 
	double lnElt = lnEmin + Ilt*dlnE; 
	//       >> Interpolate
	double flux_r = PhotonFluxMap.at("dNdy_Table_" + Case)[Ilt-1] + ((lEgamma-lnElt)/dlnE)*(PhotonFluxMap.at("dNdy_Table_" + Case)[Ilt]- PhotonFluxMap.at("dNdy_Table_" + Case)[Ilt-1]);
	flux_r = flux_r/Egamma;

	return flux_r;
}

void InterpolateFlux(std::map<TString, std::vector<double>> &Map, TString inFileDir = "flux/", TString subCase = "")
{
	cout<<"InterpolateFlux-------->Raps:	";
	for (int i = 0; i < Map.at("Raps").size(); ++i)
	{
		cout<<Map.at("Raps")[i]<<"	";
	}
	cout<<endl;

	const double JpsiMass = 3.096916;
	std::vector<TString> CasesName = {"AnAn",	"0n0n", "0nXnSum", "XnXn"};

	// if (PhotonFluxMap.at("Energy_Table_AnAn").size() == 0){	loadPhotonFlux(inFileDir, subCase);	}
	PhotonFluxMap = PhotonFluxMapTemp;
	loadPhotonFlux(inFileDir, subCase);

	//Calculate photon energy from the rap
	for (int i = 0; i < CasesName.size(); ++i)
	{
		TString Case = CasesName[i];
		for (int j = 0; j < Map.at("Raps").size(); ++j)
		{
			double k = JpsiMass/2 * exp(Map.at("Raps")[j]);
			Map["Energy"].push_back( k );
			Map["dNdy_" + Case].push_back(	k * Interpolate(k, Case)	);
		}
	}

	cout<<"InterpolateFlux-------->DONE"<<endl;
}

std::vector<double> getFluxUncer(std::vector<double> Default, std::vector<double> Diff)
{
	std::vector<double> Uncer{};
	for (int i = 0; i < Default.size(); ++i)
	{
		Uncer.push_back( abs( Default[i]-Diff[i] ) / Default[i] * 100);
	}
	return Uncer;
}

void CompareFlux(TString Case,	std::vector<TString> subCases, std::vector<TString> legendName, TString inFileDir = "flux/")
{
	loadPhotonFlux(inFileDir, "_SigNN68p3R6p67a0p56");
	auto PhotonFluxMap_Default = PhotonFluxMap;
	PhotonFluxMap = PhotonFluxMapTemp;

	std::vector<int> colors{1,2,3,4};

	TCanvas *c = new TCanvas();
	// c->SetLogy();
	TH2D* htem2d = new TH2D("htem2d_"+Case, "_Uncer;y;dN/dy Uncer. (%)", 10,-4,4, 10, 0, 20);
	htem2d->Draw();
	auto legend = new TLegend(0.2, 0.6, 0.5, 0.8);

	for (int i = 0; i < subCases.size(); ++i)
	{
		auto subCase = subCases[i];
		loadPhotonFlux(inFileDir, subCase);
		auto PhotonFluxMap_subCase = PhotonFluxMap;
		PhotonFluxMap = PhotonFluxMapTemp;

		auto Uncer_subCase = getFluxUncer(PhotonFluxMap_Default.at("dNdy_Table_"+Case),	PhotonFluxMap_subCase.at("dNdy_Table_"+Case));
		TGraph* gr_subCase = new TGraph(PhotonFluxMap_subCase.at("Raps_Table_"+Case).size(), & PhotonFluxMap_subCase.at("Raps_Table_"+Case)[0], & Uncer_subCase[0]);
		
		gr_subCase->SetLineColor(colors[i]);
		gr_subCase->SetLineWidth(3);

		gr_subCase->Draw("SAME l");
		legend->AddEntry(gr_subCase, 	legendName[i],	"l");
	}

	drawLatex(0.3, 0.85, "UPC Pb+Pb #sqrt{s_{NN}} = 5.02 TeV (" + Case +")",      42,       0.05,      1);
	legend->Draw();

	c->SaveAs( "out4Flux/CompareFlux_" + Case + "_dNdy.png" );
	c->SaveAs( "out4Flux/CompareFlux_" + Case + "_dNdy.pdf" );
	delete c;
	delete htem2d;
}

void getPhotonFlux()
{
	//------------------Standard Test-------------------------------------------
	std::map<TString, std::vector<double>>	TestMap = 
	{
		{"Raps",	{ 1.75, -1.75, 2, -2, 2.25, -2.25 }},
		{"dNdy",	{}}
	};
	InterpolateFlux(TestMap, "flux/", "_SigNN68p3R6p67a0p56");
	plotFlux(TestMap,	"0n0n");
	plotFlux(TestMap,	"0nXnSum");
	plotFlux(TestMap,	"XnXn");
	plotFlux(TestMap,	"AnAn");
	// plotPofB();
	//--------------------------------------------------------------------------

	//--------------------------------------------------------------------------
	//------------------Flux For Drawing ALICE and LHCb values------------------
	// std::map<TString, std::vector<double>>	TestMap = 
	// {
	// 	{"Raps",	{ 0, -3.875, -3.625, -4.25, -3.75 }},
	// 	{"dNdy",	{}}
	// };

	// auto TestMap1 = TestMap;
	// auto TestMap2 = TestMap;
	// InterpolateFlux(TestMap, "flux/", "_SigNN68p3R6p67a0p56");
	// plotFlux(TestMap,	"0n0n");
	// plotFlux(TestMap,	"0nXnSum");
	// plotFlux(TestMap,	"XnXn");
	// plotFlux(TestMap,	"AnAn");

	// InterpolateFlux(TestMap1, "flux/", "_SigNN68p3R6p64a0p53");
	// InterpolateFlux(TestMap2, "flux/", "_SigNN68p3R6p70a0p59");
	// auto v1 = getFluxUncer(TestMap.at("dNdy_AnAn"),	TestMap1.at("dNdy_AnAn"));
	// auto v2 = getFluxUncer(TestMap.at("dNdy_AnAn"),	TestMap2.at("dNdy_AnAn"));
	// for (int i = 0; i < v1.size(); ++i)
	// {
	// 	auto uncer = (v1[i] > v2[i]) ? v1[i]: v2[i];
	// 	cout << "Raps:" << TestMap.at("Raps")[i] << "	Flux:" << TestMap.at("dNdy_AnAn")[i]; cout<< "	Uncer:" << uncer <<endl;
	// }

	// CompareFlux("AnAn",		{"_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59", "_SigNN68p3R6p64a0p59", "_SigNN68p3R6p70a0p53"},	{"#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59", "#sigma_{NN} = 68.3, R = 6.64, a = 0.59", "#sigma_{NN} = 68.3, R = 6.70, a = 0.53"});
	// CompareFlux("0n0n",		{"_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59", "_SigNN68p3R6p64a0p59", "_SigNN68p3R6p70a0p53"},	{"#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59", "#sigma_{NN} = 68.3, R = 6.64, a = 0.59", "#sigma_{NN} = 68.3, R = 6.70, a = 0.53"});
	// CompareFlux("0nXnSum",	{"_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59", "_SigNN68p3R6p64a0p59", "_SigNN68p3R6p70a0p53"},	{"#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59", "#sigma_{NN} = 68.3, R = 6.64, a = 0.59", "#sigma_{NN} = 68.3, R = 6.70, a = 0.53"});
	// CompareFlux("XnXn",		{"_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59", "_SigNN68p3R6p64a0p59", "_SigNN68p3R6p70a0p53"},	{"#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59", "#sigma_{NN} = 68.3, R = 6.64, a = 0.59", "#sigma_{NN} = 68.3, R = 6.70, a = 0.53"});

	// CompareFlux("AnAn",		{"_SigNN68p3R6p67a0p56", "_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59"},	{"#sigma_{NN} = 68.3, R = 6.67, a = 0.56", "#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59"});
	// CompareFlux("0n0n",		{"_SigNN68p3R6p67a0p56", "_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59"},	{"#sigma_{NN} = 68.3, R = 6.67, a = 0.56", "#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59"});
	// CompareFlux("0nXnSum",	{"_SigNN68p3R6p67a0p56", "_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59"},	{"#sigma_{NN} = 68.3, R = 6.67, a = 0.56", "#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59"});
	// CompareFlux("XnXn",		{"_SigNN68p3R6p67a0p56", "_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59"},	{"#sigma_{NN} = 68.3, R = 6.67, a = 0.56", "#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59"});
	//--------------------------------------------------------------------------
}
