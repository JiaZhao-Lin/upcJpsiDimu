#include "Map_IO.C"
#include "PhysParameters.h"
#include "../simulation/getPhotonFlux.C"


double y2W(const double y)
{
	const double JpsiMass   = 3.096916;
	const double Gamma_beam = 2672.9;
	const double Mass_N     = (0.93827+0.93957) / 2;
	return sqrt( (2 * Gamma_beam * Mass_N * JpsiMass) * exp(y) );
}


TF1 * fitSigmas(const int data_option = 3, const int func_option = 0)
{
	auto ShadowRatio_ParamsMap 	= readMap("rootfiles/Results_Map.root");
	auto TotalSysUncer_Map 		= readMap("rootfiles/TotalSysUncer_Map.root");

	auto Ws_ToFit 			= ShadowRatio_ParamsMap["Ws"];
	auto Sigmas_ToFit 		= ShadowRatio_ParamsMap["Sigmas"];
	auto Sigmas_Err_ToFit 	= ShadowRatio_ParamsMap["Sigmas_Err"];

	// std::vector<double> Ws_ToFit 			= {};
	// std::vector<double> Sigmas_ToFit 		= {};
	// std::vector<double> Sigmas_Err_ToFit 	= {};

	TF1* funcSigmas;

	switch(func_option) {
		case 1:
			funcSigmas = new TF1("SigmasPrediction","[0] * x/([1]+[2]*x)");
			funcSigmas->SetParameters(0.001, 10, 1);
			break;
		case 2:
			funcSigmas = new TF1("SigmasPrediction","[0] * pow(x,[1]) / (pow(x,[1]) + [2])");
			funcSigmas->SetParameters(0.01, 1, 10);
			break;
		case 3: //Asymptotic regression
			funcSigmas = new TF1("SigmasPrediction","[0] - ( [0]-[1] ) * exp(-[2]*x)");
			funcSigmas->SetParameters(0.01, 1, 0.01);
			break;
		case 4:	//Logistic curve
			funcSigmas = new TF1("SigmasPrediction","[0] + ([1]-[0])/(1+exp([2]*(x-[3])))");
			funcSigmas->SetParameters(-0.1, 0.01,-0.01,100);
			break;
		case 5:	//Gompertz Curve
			funcSigmas = new TF1("SigmasPrediction","[0] + ([1]-[0])*(1-exp(-exp([2]*(x-[3]))))");
			funcSigmas->SetParameters(-0.1, 0.01,-0.01,100);
			break;
		case 6:	//Alt. Gompertz Curve
			funcSigmas = new TF1("SigmasPrediction","[0] + ([1]-[0])*exp(-exp([2]*(x-[3])))");
			funcSigmas->SetParameters(-0.1, 0.01,-0.01,100);
			break;
		case 7:	//Log-logistic curve
			funcSigmas = new TF1("SigmasPrediction","[0] + ([1]-[0])/( 1+exp( [2]*(log(x)-log([3])) ) )");
			funcSigmas->SetParameters(0.01, 0.01,-1,10);
			break;
		default: // Err Funcion
			funcSigmas = new TF1("SigmasPrediction","[0]*TMath::Erf((x-[1])/[2])");
			funcSigmas->SetParameters(0.1,10,10);
	}
	cout<<endl;
	funcSigmas->Print();
	cout<<endl;

	for (int i = 0; i < ALICE_Run2_FwdRap_y.size(); ++i)
	{
		//Calculating the values
		ALICE_Run2_FwdRap_W        		.push_back( y2W(ALICE_Run2_FwdRap_y[i]) );
		ALICE_Run2_FwdRap_Sigma    		.push_back( ALICE_Run2_FwdRap_dSigmady[i] / ALICE_Run2_FwdRap_Flux[i] );
		ALICE_Run2_FwdRap_Sigma_StatErr	.push_back( ALICE_Run2_FwdRap_Sigma[i] * ALICE_Run2_FwdRap_dSigmady_StatErr[i] / ALICE_Run2_FwdRap_dSigmady[i]		);
	}
	TGraphErrors* ge_ALICE_Run2_MidRap	= new TGraphErrors(ALICE_Run2_MidRap_W.size(),	&ALICE_Run2_MidRap_W[0],	&ALICE_Run2_MidRap_Sigma[0],	0,	&ALICE_Run2_MidRap_Sigma_StatErr[0]	);
	TGraphErrors* ge_ALICE_Run2_FwdRap	= new TGraphErrors(ALICE_Run2_FwdRap_W.size(),	&ALICE_Run2_FwdRap_W[0],	&ALICE_Run2_FwdRap_Sigma[0],	0,	&ALICE_Run2_FwdRap_Sigma_StatErr[0]	);


	switch(data_option) {
		case 1:
			Ws_ToFit.insert(Ws_ToFit.end(),	ALICE_Run2_FwdRap_W.begin(),	ALICE_Run2_FwdRap_W.end());
			Sigmas_ToFit.insert(Sigmas_ToFit.end(),	ALICE_Run2_FwdRap_Sigma.begin(),	ALICE_Run2_FwdRap_Sigma.end());
			Sigmas_Err_ToFit.insert(Sigmas_Err_ToFit.end(),	ALICE_Run2_FwdRap_Sigma_StatErr.begin(),	ALICE_Run2_FwdRap_Sigma_StatErr.end());
			break;
		case 2:
			Ws_ToFit.insert(Ws_ToFit.end(),	ALICE_Run2_MidRap_W.begin(),	ALICE_Run2_MidRap_W.end());
			Sigmas_ToFit.insert(Sigmas_ToFit.end(),	ALICE_Run2_MidRap_Sigma.begin(),	ALICE_Run2_MidRap_Sigma.end());
			Sigmas_Err_ToFit.insert(Sigmas_Err_ToFit.end(),	ALICE_Run2_MidRap_Sigma_StatErr.begin(),	ALICE_Run2_MidRap_Sigma_StatErr.end());
			break;
		case 3:
			Ws_ToFit.insert(Ws_ToFit.end(),	ALICE_Run2_FwdRap_W.begin(),	ALICE_Run2_FwdRap_W.end());
			Ws_ToFit.insert(Ws_ToFit.end(),	ALICE_Run2_MidRap_W.begin(),	ALICE_Run2_MidRap_W.end());
			Sigmas_ToFit.insert(Sigmas_ToFit.end(),	ALICE_Run2_FwdRap_Sigma.begin(),	ALICE_Run2_FwdRap_Sigma.end());
			Sigmas_ToFit.insert(Sigmas_ToFit.end(),	ALICE_Run2_MidRap_Sigma.begin(),	ALICE_Run2_MidRap_Sigma.end());
			Sigmas_Err_ToFit.insert(Sigmas_Err_ToFit.end(),	ALICE_Run2_FwdRap_Sigma_StatErr.begin(),	ALICE_Run2_FwdRap_Sigma_StatErr.end());
			Sigmas_Err_ToFit.insert(Sigmas_Err_ToFit.end(),	ALICE_Run2_MidRap_Sigma_StatErr.begin(),	ALICE_Run2_MidRap_Sigma_StatErr.end());
			break;
		default:;
	}

	auto c = new TCanvas();
	c->SetLogy();
	auto SigmaVsW = new TH2D("SigmaVsW", ";W_{#gammaPb} (GeV);#sigma(#gamma A #rightarrow J/#psi A) (mb);", 10,0,420, 10, 0.005, 0.20);

	SigmaVsW->GetYaxis()->CenterTitle();
	SigmaVsW->GetYaxis()->SetTitleSize(0.065);
	SigmaVsW->GetYaxis()->SetTitleSize(0.065);
	SigmaVsW->GetYaxis()->SetTitleOffset(0.85);
	SigmaVsW->GetYaxis()->SetLabelSize(0.04);
	SigmaVsW->GetXaxis()->SetTitleSize(0.05);
	SigmaVsW->GetXaxis()->SetTitleOffset(0.98);
	SigmaVsW->GetXaxis()->SetLabelSize(0.04);
	SigmaVsW->SetTickLength(0.04);
	SigmaVsW->Draw();
	drawLatex(0.15, 0.85, "UPC Pb+Pb #sqrt{s_{NN}} = 5.02 TeV ()",      42,       0.05,      1);
	drawLatex(0.15, 0.78, Form("Fit Func: %s",funcSigmas->GetFormula()->GetExpFormula().Data()),      42,       0.03,      1);

	cout<<"fitSigmas----------------->Proceessing<----------------"<<endl;
	TGraphErrors* gr = new TGraphErrors(Ws_ToFit.size(), &Ws_ToFit[0], &Sigmas_ToFit[0],0, &Sigmas_Err_ToFit[0]);
	



	gr->Fit(funcSigmas);
	gr->SetTitle("TGraph1D TF1 Fit; dN_{1}/dy; dN_{2}/dy; d#sigma/dy");
	gr->SetMarkerColor(kBlue);
	// gr->SetMarkerSize(0.8);
	gr->SetLineWidth(5);
	gr->Draw("same pez");

	// ge_ALICE_Run2_MidRap->SetMarkerStyle(24);
	// ge_ALICE_Run2_MidRap->SetMarkerColor(4);
	// ge_ALICE_Run2_MidRap->SetLineColor(4);
	// ge_ALICE_Run2_MidRap->SetLineWidth(2);
	// ge_ALICE_Run2_MidRap->Draw("pezsame");

	// ge_ALICE_Run2_FwdRap->SetMarkerStyle(25);
	// ge_ALICE_Run2_FwdRap->SetMarkerColor(4);
	// ge_ALICE_Run2_FwdRap->SetLineColor(4);
	// ge_ALICE_Run2_FwdRap->SetLineWidth(2);
	// ge_ALICE_Run2_FwdRap->Draw("pezsame");

	// funcSigmas->SetMarkerColor(kBlue);
	// funcSigmas->Draw("same surf");

	// drawLatex(0.05, 0.95, Form("(#sigma(y = %.2f), #sigma(y = %.2f)) = (%.4f #pm %.4f, %.4f #pm %.4f)", Rap, -Rap,
	//                       f->GetParameter(0), f->GetParError(0), f->GetParameter(1), f->GetParError(1)),      42,       0.05,      1);

	// c->SaveAs(Form("outplots/fit2D_%d.png",i));
	cout<<endl;

	// Sigmas      .push_back(f->GetParameter(0));  Sigmas      .push_back(f->GetParameter(1));
	// Sigmas_Err  .push_back(f->GetParError(0));   Sigmas_Err  .push_back(f->GetParError(1));

	// delete gr; delete c; delete f;
	cout<<"fitSigmas----------------->DONE<----------------"<<endl;

	return funcSigmas;
}

void predictDSigmaDy()
{
	int nSteps = 40;
	double y_step = 0.1;

	TF1 * funcSigmas	=	fitSigmas();

	std::vector<double> Raps_Prediction = {};
	std::vector<double> Ws_Prediction = {};
	std::vector<double> Xsec_AnAn_Prediction = {};
	std::vector<double> Sigmas_Prediction = {};

	for (int i = 0; i < nSteps; ++i)
	{
		auto y_plus	 = i*y_step;
		auto y_minus = - y_plus;

		auto W_plus		= y2W(y_plus);
		auto W_minus 	= y2W(y_minus);

		Raps_Prediction	.push_back(	y_plus	);
		Raps_Prediction	.insert(Raps_Prediction.begin(),	y_minus	);
		Ws_Prediction	.push_back(	W_plus	);
		Ws_Prediction	.insert(Ws_Prediction.begin(),	W_minus	);

		Sigmas_Prediction.push_back(	funcSigmas->Eval(W_plus)	);
		Sigmas_Prediction.insert(Sigmas_Prediction.begin(),		funcSigmas->Eval(W_minus)	);

	}

	std::map<TString, std::vector<double>> Temp_Map	=	
	{
		{"Raps",	Raps_Prediction},
		{"Ws",		Ws_Prediction}
	};

	InterpolateFlux(Temp_Map,	"../simulation/flux/");

	for (int i = 0; i < nSteps; ++i)
	{
		auto Sigmas_plus	= Sigmas_Prediction[i];
		auto Sigmas_minus 	= Sigmas_Prediction[Sigmas_Prediction.size()-1-i];
		auto dNdy_plus		= Temp_Map.at("dNdy_AnAn")[i];
		auto dNdy_minus		= Temp_Map.at("dNdy_AnAn")[Sigmas_Prediction.size()-1-i];

		auto Xsec = dNdy_minus*Sigmas_minus + dNdy_plus*Sigmas_plus;
		Xsec_AnAn_Prediction.push_back(Xsec);
	}
	Temp_Map["Xsec_AnAn_Prediction"] = Xsec_AnAn_Prediction;

	// auto c1 = new TCanvas();
	// TH2D* htem2d = new TH2D("AnAn", "AnAn;y;d#sigma_{J/#psi}/dy (mb);", 10, -4.1, 0, 10, 0, 8.0);
	// htem2d ->GetYaxis()->CenterTitle();
	// //htem2d ->GetYaxis()->SetNdivisions(6);
	// htem2d ->GetYaxis()->SetTitleSize(0.07);
	// htem2d ->GetYaxis()->SetTitleOffset(0.63);
	// htem2d ->GetYaxis()->SetLabelSize(0.06);
	// htem2d ->GetXaxis()->CenterTitle();
	// htem2d ->GetXaxis()->SetTitleSize(0.07);
	// htem2d ->GetXaxis()->SetTitleOffset(0.75);
	// htem2d ->GetXaxis()->SetLabelSize(0.06);
	// //htem2d ->SetTickLength(0.08);

	// htem2d->Draw();

	// auto gr = new TGraph(Xsec_AnAn_Prediction.size(),Raps_Prediction.data(),Xsec_AnAn_Prediction.data());
	// gr->Draw("Same");

	cout<<"predictDSigmaDy----------------->DSigmaDy<----------------"<<endl;
	saveMap(Temp_Map,	"rootfiles/DataDrivenPrediction.root");
	cout<<"predictDSigmaDy----------------->DONE Saving<----------------"<<endl;
}

