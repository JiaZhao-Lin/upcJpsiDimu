#include "Map_IO.C"
#include "PhysParameters.h"
#include "../simulation/getPhotonFlux.C"

	const int template_option  = 1; //0:OldCohJpsi; 1:NewCohJpsi (w R+1fm);
	const TString template_Name[2] = {"", "_NewCohJpsi"};

double y2W(const double y)
{
	const double JpsiMass   = 3.096916;
	const double Gamma_beam = 2672.9;
	const double Mass_N     = (0.93827+0.93957) / 2;
	return sqrt( (2 * Gamma_beam * Mass_N * JpsiMass) * exp(y) );
}

double fPieceWise4Fit(double *x, double *par)
{
	if((x[0]>42))
	{
		return exp( par[2]*x[0] + par[3]);
	}

	return exp( par[0]*x[0] + par[1]);

}

TF1 * fitSigmas(const int data_option = 0, const int func_option = 9)
{
	auto ShadowRatio_ParamsMap 	= readMap(Form("rootfiles/Results%s_Map.root", template_Name[template_option].Data()));
	auto TotalSysUncer_Map 		= readMap(Form("rootfiles/TotalSysUncer%s_Map.root", template_Name[template_option].Data()));

	auto Ws_ToFit 			= ShadowRatio_ParamsMap["Ws"];
	auto Sigmas_ToFit 		= ShadowRatio_ParamsMap["Sigmas"];
	auto Sigmas_Err_ToFit 	= ShadowRatio_ParamsMap["Sigmas_Err"];

	// --------------------------For moving points up and down-----------------------------------
	// for (int i = 0; i < Ws_ToFit.size(); i++)
	// {
	// 	// cout<<Ws_ToFit[i]<<endl;
	// 	// if (i%2==1)
	// 	// {
	// 	// 	Sigmas_ToFit[i]-= 0.01 * TotalSysUncer_Map.at("Sigmas")[i]*TotalSysUncer_Map.at("Sigmas_TotalSysUncer")[i];
	// 	// }
	// 	// else
	// 	// {
	// 	// 	Sigmas_ToFit[i]+= 0.01 * TotalSysUncer_Map.at("Sigmas")[i]*TotalSysUncer_Map.at("Sigmas_TotalSysUncer")[i];
	// 	// }
	// 	Sigmas_ToFit[i]-= 0.01 * TotalSysUncer_Map.at("Sigmas")[i]*TotalSysUncer_Map.at("Sigmas_TotalSysUncer")[i];
	// }
	// ------------------------------------------------------------------------------------------


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
		case 8:	//piecewise
			funcSigmas = new TF1("SigmasPrediction", fPieceWise4Fit, 0, 400,4);
			funcSigmas->SetParameters(0.01, 0.01,0.01,0.01);
			funcSigmas->SetParNames("a1", "b1", "a2", "b2");
			// fJpsiPeak->SetParameter(1,1);
			// fJpsiPeak->SetParameter(2,1);
			// fJpsiPeak->SetParameter(3,1);
			// fJpsiPeak->SetParameter(4,3.1);
			// funcSigmas = new TF1("SigmasPrediction","[0] + ([1]-[0])/( 1+exp( [2]*(log(x)-log([3])) ) )");
			// funcSigmas->SetParameters(0.01, 0.01,-1,10);
			break;
		case 9:	//Linear
			funcSigmas = new TF1("SigmasPrediction", "[0] * x + [1]", 42,400);
			funcSigmas->SetParameters(0.001, 1);
			funcSigmas->SetParNames("a", "b");
			break;
		case 10: //Log
			funcSigmas = new TF1("SigmasPrediction", "exp( [0]*x + [1])", 42,400);
			funcSigmas->SetParameters(0.001, 0.1);
			funcSigmas->SetParNames("a", "b");
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
		ALICE_Run2_FwdRap_W        		.push_back( y2W(ALICE_Run2_FwdRap_y[i]) );
		ALICE_Run2_FwdRap_Sigma    		.push_back( ALICE_Run2_FwdRap_dSigmady[i] / ALICE_Run2_FwdRap_Flux[i] );
		ALICE_Run2_FwdRap_Sigma_StatErr	.push_back( ALICE_Run2_FwdRap_Sigma[i] * ALICE_Run2_FwdRap_dSigmady_StatErr[i] / ALICE_Run2_FwdRap_dSigmady[i]		);
	}
	for (int i = 0; i < LHCb_Run2_FwdRap_y.size(); ++i)
	{
		LHCb_Run2_FwdRap_W        		.push_back( y2W(LHCb_Run2_FwdRap_y[i]) );
		LHCb_Run2_FwdRap_Sigma    		.push_back( LHCb_Run2_FwdRap_dSigmady[i] / LHCb_Run2_FwdRap_Flux[i] );
		LHCb_Run2_FwdRap_Sigma_StatErr	.push_back( LHCb_Run2_FwdRap_Sigma[i] * LHCb_Run2_FwdRap_dSigmady_StatErr[i] / LHCb_Run2_FwdRap_dSigmady[i]		);
	}

	switch(data_option) {
		case 1:	// use ALICE LHCb Fwd and CMS Data points 
			Ws_ToFit.insert(Ws_ToFit.end(),	ALICE_Run2_FwdRap_W.begin(),	ALICE_Run2_FwdRap_W.end());
			Ws_ToFit.insert(Ws_ToFit.end(),	LHCb_Run2_FwdRap_W.begin(),	LHCb_Run2_FwdRap_W.end());
			Sigmas_ToFit.insert(Sigmas_ToFit.end(),	ALICE_Run2_FwdRap_Sigma.begin(),	ALICE_Run2_FwdRap_Sigma.end());
			Sigmas_ToFit.insert(Sigmas_ToFit.end(),	LHCb_Run2_FwdRap_Sigma.begin(),	LHCb_Run2_FwdRap_Sigma.end());
			Sigmas_Err_ToFit.insert(Sigmas_Err_ToFit.end(),	ALICE_Run2_FwdRap_Sigma_StatErr.begin(),	ALICE_Run2_FwdRap_Sigma_StatErr.end());
			Sigmas_Err_ToFit.insert(Sigmas_Err_ToFit.end(),	LHCb_Run2_FwdRap_Sigma_StatErr.begin(),	LHCb_Run2_FwdRap_Sigma_StatErr.end());
			break;
		case 2:	// use ALICE Mid and CMS Data points
			Ws_ToFit.insert(Ws_ToFit.end(),	ALICE_Run2_MidRap_W.begin(),	ALICE_Run2_MidRap_W.end());
			Sigmas_ToFit.insert(Sigmas_ToFit.end(),	ALICE_Run2_MidRap_Sigma.begin(),	ALICE_Run2_MidRap_Sigma.end());
			Sigmas_Err_ToFit.insert(Sigmas_Err_ToFit.end(),	ALICE_Run2_MidRap_Sigma_StatErr.begin(),	ALICE_Run2_MidRap_Sigma_StatErr.end());
			break;
		case 3:	// use ALICE Fwd, Mid and CMS Data points
			Ws_ToFit.insert(Ws_ToFit.end(),	ALICE_Run2_FwdRap_W.begin(),	ALICE_Run2_FwdRap_W.end());
			Ws_ToFit.insert(Ws_ToFit.end(),	ALICE_Run2_MidRap_W.begin(),	ALICE_Run2_MidRap_W.end());
			Sigmas_ToFit.insert(Sigmas_ToFit.end(),	ALICE_Run2_FwdRap_Sigma.begin(),	ALICE_Run2_FwdRap_Sigma.end());
			Sigmas_ToFit.insert(Sigmas_ToFit.end(),	ALICE_Run2_MidRap_Sigma.begin(),	ALICE_Run2_MidRap_Sigma.end());
			Sigmas_Err_ToFit.insert(Sigmas_Err_ToFit.end(),	ALICE_Run2_FwdRap_Sigma_StatErr.begin(),	ALICE_Run2_FwdRap_Sigma_StatErr.end());
			Sigmas_Err_ToFit.insert(Sigmas_Err_ToFit.end(),	ALICE_Run2_MidRap_Sigma_StatErr.begin(),	ALICE_Run2_MidRap_Sigma_StatErr.end());
			break;
		default: // use CMS Data points ONLY
			break;
	}

	auto c = new TCanvas();
	c->SetLogy();
	auto SigmaVsW = new TH2D("SigmaVsW", ";W_{#gammaN}^{Pb} (GeV);#sigma(#gamma Pb #rightarrow J/#psi Pb) (mb);", 10,0,420, 10, 0.004, 0.13);

	SigmaVsW->GetYaxis()->CenterTitle();
    SigmaVsW->GetXaxis()->CenterTitle();
	SigmaVsW->GetYaxis()->SetTitleSize(0.05);
	SigmaVsW->GetYaxis()->SetTitleSize(0.05);
	SigmaVsW->GetYaxis()->SetTitleOffset(0.99);
	SigmaVsW->GetYaxis()->SetLabelSize(0.04);
	SigmaVsW->GetXaxis()->SetTitleSize(0.05);
	SigmaVsW->GetXaxis()->SetTitleOffset(0.98);
	SigmaVsW->GetXaxis()->SetLabelSize(0.04);
	SigmaVsW->SetTickLength(0.04);
	SigmaVsW->Draw();

    drawLatex(0.12,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);	
	// drawLatex(0.15, 0.78, Form("Fit Func: %s",funcSigmas->GetFormula()->GetExpFormula().Data()),      42,       0.03,      1);
	drawLatex(0.15, 0.80, "Fit Func: y = a*x + b",      42,       0.04,      1);
	drawLatex(0.15, 0.70, "Fitting to the data moved to lower syst. error",      42,       0.03,      1);

	cout<<"fitSigmas----------------->Proceessing<----------------"<<endl;
	TGraphErrors* gr = new TGraphErrors(Ws_ToFit.size(), Ws_ToFit.data(), Sigmas_ToFit.data(),0, Sigmas_Err_ToFit.data());
	
	funcSigmas->SetLineStyle(2);
	gr->Fit(funcSigmas, "R");
	gr->SetTitle("TGraph1D TF1 Fit; dN_{1}/dy; dN_{2}/dy; d#sigma/dy");
	gr->SetMarkerColor(kBlue);
	// gr->SetMarkerSize(0.8);
	gr->SetLineWidth(2);
	gr->Draw("same pez");
	cout<<endl;

	// TPaveStats *st = (TPaveStats*)gr->FindObject("stats");
	// st->SetX1NDC(0.8); //new x start position
	// st->SetX2NDC(0.9); //new x end position

	c->SaveAs("outplots/SigmasVsW_Fit.pdf");
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

	InterpolateFlux(Temp_Map,	"../simulation/flux/", "_SigNN68p3R6p67a0p56");

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

	// cout<<"predictDSigmaDy----------------->DSigmaDy<----------------"<<endl;
	// saveMap(Temp_Map,	Form("rootfiles/DataDrivenPrediction%s.root", template_Name[template_option].Data() ));
	// cout<<"predictDSigmaDy----------------->DONE Saving<----------------"<<endl;
}

