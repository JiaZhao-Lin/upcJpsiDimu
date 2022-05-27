const double JpsiMass   = 3.096916;
const double Sqrt_s     = 5020;
const double Gamma_beam = 2672.9;
const double Mass_N     = (0.93827+0.93957) / 2;

const std::vector<double>  CGC_JpsiNoFluct_W	= {37, 47, 59,75, 95, 121, 153, 194, 246, 311, 394, 599  },
						   CGC_JpsiNoFluct_CohXsec 		= {0.01759987558438179, 0.020489530220598236, 0.023324292435666014, 0.02714433695284619, 0.03137306416640146, 0.036389263389293006, 0.041829124040095714, 0.04793443973914736, 0.05473231107873076, 0.06212140529983476, 0.07032820925862308, 0.0874426817784331},
						   CGC_JpsiNoFluct_InCohXsec	= {0.005314082049208781, 0.006081113144220922, 0.00684657722834769, 0.008057272329685066, 0.009032600408812817, 0.010363958837189501, 0.01173848797350525, 0.013306012198701258, 0.014398401192654236, 0.016650908354141573, 0.01846420961715733, 0.022516649102744427};

//------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------
//-------------------------------------------Run2 ALICE MidRap------------------------------------------------
const std::vector<double> ALICE_Run2_MidRap_W  				= 	{ 124.67672		},
						ALICE_Run2_MidRap_Sigma	    		= 	{ 4.07/2/81.026 },
						ALICE_Run2_MidRap_Sigma_StatErr		=	{ ALICE_Run2_MidRap_Sigma[0] * 0.11/4.07 },
						ALICE_Run2_MidRap_Sigma_SysErrLow	=	{ ALICE_Run2_MidRap_Sigma[0] * TMath::Hypot( 0.23/4.07,	 0.05)	}, // 5% error for flux
						ALICE_Run2_MidRap_Sigma_SysErrHig	=	{ ALICE_Run2_MidRap_Sigma[0] * TMath::Hypot( 0.23/4.07,	 0.05)	}, // 5% error for flux
						ALICE_Run2_MidRap_x	    			= 	{ 0.00061691554 },
						ALICE_Run2_MidRap_Sigma_IA	    	= 	{ 0.069926461 	},
						ALICE_Run2_MidRap_Sigma_IA_Err		= 	{ 0.0037760289 	},
						ALICE_Run2_MidRap_R	    			= 	{ sqrt(ALICE_Run2_MidRap_Sigma[0] / ALICE_Run2_MidRap_Sigma_IA[0]) },
						ALICE_Run2_MidRap_R_StatErr	    	= 	{ 0.5 * ALICE_Run2_MidRap_R[0] * ALICE_Run2_MidRap_Sigma_StatErr[0]/ALICE_Run2_MidRap_Sigma[0] },
						ALICE_Run2_MidRap_R_SysErrLow  		= 	{ 0.5 * ALICE_Run2_MidRap_R[0] * TMath::Hypot( ALICE_Run2_MidRap_Sigma_SysErrLow[0]/ALICE_Run2_MidRap_Sigma[0],
																								ALICE_Run2_MidRap_Sigma_IA_Err[0]/ALICE_Run2_MidRap_Sigma_IA[0]) },
						ALICE_Run2_MidRap_R_SysErrHig  		= 	{ 0.5 * ALICE_Run2_MidRap_R[0] * TMath::Hypot( ALICE_Run2_MidRap_Sigma_SysErrHig[0]/ALICE_Run2_MidRap_Sigma[0],
																								ALICE_Run2_MidRap_Sigma_IA_Err[0]/ALICE_Run2_MidRap_Sigma_IA[0]) };
//-------------------------------------------Run2 ALICE FwdRap------------------------------------------------
const std::vector<double> ALICE_Run2_FwdRap_y  				= 	{ -3.875,		-3.625 	},
						ALICE_Run2_FwdRap_dSigmady			=	{ 1.615,		1.938 	},
						ALICE_Run2_FwdRap_dSigmady_StatErr	=	{ 0.060,		0.0420 	},
						ALICE_Run2_FwdRap_dSigmady_SysErrLow=	{ 0.147,		0.190 	},
						ALICE_Run2_FwdRap_dSigmady_SysErrHig=	{ 0.135,		0.166 	},
						ALICE_Run2_FwdRap_Flux 				=	{ 201.308,		193.547 },
						ALICE_Run2_FwdRap_Flux_Err 			=	{ 201.308*0.05,	193.547*0.05}, // 5% error for flux
						ALICE_Run2_FwdRap_Sigma_IA 			=	{ 0.011413847,	0.01347654		},
						ALICE_Run2_FwdRap_Sigma_IA_Err 		=	{ 0.00061634772,0.00072773330	};

std::vector<double> 	ALICE_Run2_FwdRap_W  				= 	{},	
						ALICE_Run2_FwdRap_Sigma				=	{},	
						ALICE_Run2_FwdRap_Sigma_StatErr		=	{},
						ALICE_Run2_FwdRap_Sigma_SysErrLow	=	{},
						ALICE_Run2_FwdRap_Sigma_SysErrHig	=	{},
						ALICE_Run2_FwdRap_x 				=	{},
						ALICE_Run2_FwdRap_R 				=	{},
						ALICE_Run2_FwdRap_R_StatErr 		=	{},
						ALICE_Run2_FwdRap_R_SysErrLow 		=	{},
						ALICE_Run2_FwdRap_R_SysErrHig 		=	{};
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

const TString outDir = "outfigure4Paper";


std::map<TString, std::vector<double>> readMap(TString fileName = "rootfiles/Results_Map.root")
{
	TFile *file = TFile::Open(fileName, "READ");
	std::map<TString, std::vector<double>> map;

	std::vector<TString> *keys;
	file->GetObject("keys", keys); // I try to retrieve the vector
	for(auto it = keys->begin(); it != keys->end(); ++it) {
		std::vector<Double_t> *temp;

		file->GetObject(*it, temp);
		std::cout << "Retrieving key to map:	" << *it << endl;
		map[*it] = *temp;
	}
	return map;
}
double y2x(const double y)
{
	return (JpsiMass / Sqrt_s) * exp(-y);
}
double y2W(const double y)
{
	return sqrt( (2 * Gamma_beam * Mass_N * JpsiMass) * exp(y) );
}
TLatex* drawLatex(double x, double y, TString text, int textFont, double textSize, int colorIndex, double textAngle=0)
{
	TLatex *latex = new TLatex(x,y,text.Data());
	latex->SetNDC();
	latex->SetTextFont(textFont);
	latex->SetTextSize(textSize);
	latex->SetTextColor(colorIndex);
	latex->SetTextAngle(textAngle);
	latex->Draw("same");
	return latex;
}
void plotSigmaVsW( std::map<TString, std::vector<double>> &ShadowRatio_ParamsMap, const std::map<TString, std::vector<double>> &TotalSysUncer_Map, const int flag4Axis = 0 ) //0: logY only, 1: logX and LogY
{
	system("mkdir -p outfigure4Paper");

	for (int i = 0; i < ALICE_Run2_FwdRap_y.size(); ++i)
	{
		//Calculating the values
		ALICE_Run2_FwdRap_W        		.push_back( y2W(ALICE_Run2_FwdRap_y[i]) );
		ALICE_Run2_FwdRap_Sigma    		.push_back( ALICE_Run2_FwdRap_dSigmady[i] / ALICE_Run2_FwdRap_Flux[i] );
		ALICE_Run2_FwdRap_Sigma_StatErr	.push_back( ALICE_Run2_FwdRap_Sigma[i] * ALICE_Run2_FwdRap_dSigmady_StatErr[i] / ALICE_Run2_FwdRap_dSigmady[i]		);
		ALICE_Run2_FwdRap_Sigma_SysErrLow.push_back( ALICE_Run2_FwdRap_Sigma[i] * TMath::Hypot(ALICE_Run2_FwdRap_dSigmady_SysErrLow[i] / ALICE_Run2_FwdRap_dSigmady[i],
																							 ALICE_Run2_FwdRap_Flux_Err[i] / ALICE_Run2_FwdRap_Flux[i])		);
		ALICE_Run2_FwdRap_Sigma_SysErrHig.push_back( ALICE_Run2_FwdRap_Sigma[i] * TMath::Hypot(ALICE_Run2_FwdRap_dSigmady_SysErrHig[i] / ALICE_Run2_FwdRap_dSigmady[i],
																							 ALICE_Run2_FwdRap_Flux_Err[i] / ALICE_Run2_FwdRap_Flux[i])		);
		
		ALICE_Run2_FwdRap_x        		.push_back( y2x(ALICE_Run2_FwdRap_y[i]) );
		ALICE_Run2_FwdRap_R 			.push_back( sqrt(ALICE_Run2_FwdRap_Sigma[i]/ALICE_Run2_FwdRap_Sigma_IA[i])	);
		ALICE_Run2_FwdRap_R_StatErr		.push_back( 0.5 * ALICE_Run2_FwdRap_R[i] * ALICE_Run2_FwdRap_Sigma_StatErr[i]/ALICE_Run2_FwdRap_Sigma[i]	);
		ALICE_Run2_FwdRap_R_SysErrLow	.push_back( 0.5 * ALICE_Run2_FwdRap_R[i] * TMath::Hypot( ALICE_Run2_FwdRap_Sigma_SysErrLow[i]/ALICE_Run2_FwdRap_Sigma[i],
																								ALICE_Run2_FwdRap_Sigma_IA_Err[i]/ALICE_Run2_FwdRap_Sigma_IA[i])	);
		ALICE_Run2_FwdRap_R_SysErrHig	.push_back( 0.5 * ALICE_Run2_FwdRap_R[i] * TMath::Hypot( ALICE_Run2_FwdRap_Sigma_SysErrHig[i]/ALICE_Run2_FwdRap_Sigma[i],
																								ALICE_Run2_FwdRap_Sigma_IA_Err[i]/ALICE_Run2_FwdRap_Sigma_IA[i])	);
	}


	auto c = new TCanvas();
	// c->SetLogx();

	//gPad->SetPad(0.0,0.25,1.0,0.96);
	gPad->SetTopMargin(0.05);
	gPad->SetBottomMargin(0.12);
	//gPad->SetLeftMargin(0.05);
	gPad->SetRightMargin(0.05);
	
	TH2D* SigmaVsW;
	
	if(flag4Axis==0) //logY only
	{
		c->SetLogy();
		SigmaVsW = new TH2D("SigmaVsW", ";W_{#gammaPb} (GeV);#sigma(#gamma A #rightarrow J/#psi A) (mb);", 10,0,420, 10, 0.005, 0.20);
	}
	else //logX and logY as default one
	{
		c->SetLogx();
		c->SetLogy();
		//SigmaVsW = new TH2D("SigmaVsW", ";W_{#gammaPb} (GeV);#sigma(#gamma A #rightarrow J/#psi A) (mb)", 10,28,520, 10, 0.015, 0.11);
		SigmaVsW = new TH2D("SigmaVsW", ";W_{#gammaPb} (GeV);#sigma(#gamma A #rightarrow J/#psi A) (mb);", 10,14,520, 10, 0.006, 0.10);
	}

	SigmaVsW->GetYaxis()->CenterTitle();
	SigmaVsW->GetYaxis()->SetTitleSize(0.065);
	SigmaVsW->GetYaxis()->SetTitleSize(0.065);
	SigmaVsW->GetYaxis()->SetTitleOffset(0.85);
	SigmaVsW->GetYaxis()->SetLabelSize(0.04);
	SigmaVsW->GetXaxis()->SetTitleSize(0.05);
	SigmaVsW->GetXaxis()->SetTitleOffset(0.98);
	SigmaVsW->GetXaxis()->SetLabelSize(0.04);
	SigmaVsW->SetTickLength(0.04);
	SigmaVsW->Draw("");

	auto X_AXIS_ERR = std::vector<double>(ShadowRatio_ParamsMap["Ws"].size(), 3);
	auto Sigmas_TotalSysErr = TotalSysUncer_Map.at("Sigmas_TotalSysUncer");
	for (int i = 0; i < Sigmas_TotalSysErr.size(); ++i)
	{
		Sigmas_TotalSysErr[i] *= 0.01 * TotalSysUncer_Map.at("Sigmas")[i];
	}
	TGraphErrors* ge_CMS        			 = new TGraphErrors(ShadowRatio_ParamsMap["Ws"].size(),	&ShadowRatio_ParamsMap["Ws"][0],	&ShadowRatio_ParamsMap["Sigmas"][0],	0,	&ShadowRatio_ParamsMap["Sigmas_Err"][0]	);
	TGraphAsymmErrors* gae_CMS 	= new TGraphAsymmErrors(
												ShadowRatio_ParamsMap["Ws"].size(),
												&ShadowRatio_ParamsMap["Ws"][0],        &ShadowRatio_ParamsMap["Sigmas"][0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&Sigmas_TotalSysErr[0], &Sigmas_TotalSysErr[0]);
	
	TGraphErrors* ge_ALICE_Run2_MidRap	= new TGraphErrors(ALICE_Run2_MidRap_W.size(),	&ALICE_Run2_MidRap_W[0],	&ALICE_Run2_MidRap_Sigma[0],	0,	&ALICE_Run2_MidRap_Sigma_StatErr[0]	);
	TGraphErrors* ge_ALICE_Run2_FwdRap	= new TGraphErrors(ALICE_Run2_FwdRap_W.size(),	&ALICE_Run2_FwdRap_W[0],	&ALICE_Run2_FwdRap_Sigma[0],	0,	&ALICE_Run2_FwdRap_Sigma_StatErr[0]	);
	TGraphAsymmErrors* gae_ALICE_Run2_MidRap 	= new TGraphAsymmErrors(
												ALICE_Run2_MidRap_W.size(),
												&ALICE_Run2_MidRap_W[0],        &ALICE_Run2_MidRap_Sigma[0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&ALICE_Run2_MidRap_Sigma_SysErrLow[0], &ALICE_Run2_MidRap_Sigma_SysErrHig[0]);
	TGraphAsymmErrors* gae_ALICE_Run2_FwdRap 	= new TGraphAsymmErrors(
												ALICE_Run2_FwdRap_W.size(),
												&ALICE_Run2_FwdRap_W[0],        &ALICE_Run2_FwdRap_Sigma[0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&ALICE_Run2_FwdRap_Sigma_SysErrLow[0], &ALICE_Run2_FwdRap_Sigma_SysErrHig[0]);

	TGraphErrors* ge_CGCnoFluct = new TGraphErrors(CGC_JpsiNoFluct_W.size(),	&CGC_JpsiNoFluct_W[0],	&CGC_JpsiNoFluct_CohXsec[0],	0,	0);
	TGraphErrors* ge_IA        	= new TGraphErrors(ShadowRatio_ParamsMap["Ws_FitIA"].size(),	&ShadowRatio_ParamsMap["Ws_FitIA"][0],	&ShadowRatio_ParamsMap["Sigmas_FitIA"][0],	0,	0);

	gae_CMS ->SetMarkerStyle(24);
	gae_CMS ->SetFillColorAlpha(16, 0.7);
	gae_CMS ->SetFillStyle(1001);
	gae_CMS ->Draw("2same");
	ge_CMS->SetMarkerStyle(20);
	ge_CMS->SetMarkerSize(1.1);
	ge_CMS->SetMarkerColor(1);
	ge_CMS->SetLineColor(1);
	ge_CMS->SetLineWidth(2);
	ge_CMS->Draw("pezsame");

	gae_ALICE_Run2_MidRap ->SetMarkerStyle(24);
	gae_ALICE_Run2_MidRap ->SetFillColorAlpha(16, 0.7);
	gae_ALICE_Run2_MidRap ->SetFillStyle(1001);
	gae_ALICE_Run2_MidRap ->Draw("2same");
	ge_ALICE_Run2_MidRap->SetMarkerStyle(24);
	ge_ALICE_Run2_MidRap->SetMarkerColor(4);
	ge_ALICE_Run2_MidRap->SetLineColor(4);
	ge_ALICE_Run2_MidRap->SetLineWidth(2);
	ge_ALICE_Run2_MidRap->Draw("pezsame");

	gae_ALICE_Run2_FwdRap ->SetMarkerStyle(24);
	gae_ALICE_Run2_FwdRap ->SetFillColorAlpha(16, 0.7);
	gae_ALICE_Run2_FwdRap ->SetFillStyle(1001);
	gae_ALICE_Run2_FwdRap ->Draw("2same");
	ge_ALICE_Run2_FwdRap->SetMarkerStyle(25);
	ge_ALICE_Run2_FwdRap->SetMarkerColor(4);
	ge_ALICE_Run2_FwdRap->SetLineColor(4);
	ge_ALICE_Run2_FwdRap->SetLineWidth(2);
	ge_ALICE_Run2_FwdRap->Draw("pezsame");

	ge_CGCnoFluct->SetMarkerColor(1);
	ge_CGCnoFluct->SetLineColor(4);
	ge_CGCnoFluct->SetLineStyle(2);
	ge_CGCnoFluct->SetLineWidth(2);
	ge_CGCnoFluct->Draw("lsame");

	ge_IA->SetMarkerColor(2);
	ge_IA->SetLineStyle(5);
	ge_IA->SetLineColor(2);
	ge_IA->SetLineWidth(2);
	ge_IA->Draw("lsame");

	//drawLatex(0.15, 0.86, "Pb+Pb #rightarrow Pb+Pb+J/#psi #sqrt{s_{NN}} = 5.02 TeV",  42,        0.05,      1 );

	if(flag4Axis == 1)
	{
		drawLatex(0.15, 0.87, "Pb+Pb UPC #sqrt{s_{NN}} = 5.02 TeV",  42,        0.05,      1 );
		
		TLegend  *legData =  new TLegend(0.13, 0.66, 0.35, 0.82);
		legData->SetFillStyle(0);
		legData->SetFillColor(0);
		legData->SetTextSize(0.045);
		legData->AddEntry(ge_CMS,               "CMS",               "p");
		legData->AddEntry(ge_ALICE_Run2_FwdRap, "ALICE (-4 < y < -3.5)", "p");
		legData->AddEntry(ge_ALICE_Run2_MidRap, "ALICE (|y| < 0.15)",  "p");
		legData->Draw("same");

		TLegend  *legTheory =  new TLegend(0.50, 0.20, 0.80, 0.30);
		legTheory->SetFillStyle(0);
		legTheory->SetFillColor(0);
		legTheory->SetTextSize(0.045);
		legTheory->AddEntry(ge_CGCnoFluct,        "CGC",                   "l");
		legTheory->AddEntry(ge_IA,                "Impulse Approximation", "l");
		legTheory->Draw("same");

		c->SaveAs("outplots/SigmaVsW_logXY.png");
		c->SaveAs("outplots/SigmaVsW_logXY.pdf");
	}
	else
	{
		drawLatex(0.15, 0.87, "Pb+Pb UPC #sqrt{s_{NN}} = 5.02 TeV",  42,        0.06,      1 );

		TLegend  *legData =  new TLegend(0.30, 0.20, 0.55, 0.35);
		legData->SetFillStyle(0);
		legData->SetFillColor(0);
		legData->SetTextSize(0.045);
		legData->AddEntry(ge_CMS,               "CMS",                   "p");
		legData->AddEntry(ge_ALICE_Run2_FwdRap, "ALICE (-4 < y < -3.5)", "p");
		legData->AddEntry(ge_ALICE_Run2_MidRap, "ALICE (|y| < 0.15)",    "p");
		legData->Draw("same");

		TLegend  *legTheory =  new TLegend(0.75, 0.20, 0.95, 0.30);
		legTheory->SetFillStyle(0);
		legTheory->SetFillColor(0);
		legTheory->SetTextSize(0.045);
		legTheory->AddEntry(ge_CGCnoFluct,        "CGC",                   "l");
		legTheory->AddEntry(ge_IA,                "IA", "l");
		legTheory->Draw("same");

		c->SaveAs(outDir+"/SigmaVsW_logY.png");
		c->SaveAs(outDir+"/SigmaVsW_logY.pdf");
	}

	delete c;
}

void plotRvsX( std::map<TString, std::vector<double>> &ShadowRatio_ParamsMap, const std::map<TString, std::vector<double>> &TotalSysUncer_Map )
{	
	auto c = new TCanvas();

	gPad->SetTopMargin(0.05);
	gPad->SetBottomMargin(0.12);
	//gPad->SetLeftMargin(0.05);
	gPad->SetRightMargin(0.05);


	c->SetLogx();

	//TH2D* htem2d = new TH2D("htem2d", "", 10,3.0e-5,5.0e-2, 10, 0, 1);
	TH2D* htem2d = new TH2D("htem2d", "", 10,4.0e-5,5e-2, 10, 0.2, 1.1);
	htem2d->SetYTitle("R^{Pb}_{g}(x, #mu^{2}=2.4 GeV^{2})");
	htem2d->SetXTitle("x");

	htem2d->GetYaxis()->SetTitleSize(0.06);
	htem2d->GetYaxis()->SetTitleOffset(0.85);
	htem2d->GetYaxis()->SetLabelSize(0.04);
	htem2d->GetXaxis()->CenterTitle();
	htem2d->GetXaxis()->SetTitleSize(0.07);
	htem2d->GetXaxis()->SetTitleOffset(0.69);
	htem2d->GetXaxis()->SetLabelSize(0.04);
	htem2d->SetTickLength(0.04);
	htem2d->Draw();

	std::vector<double> X_AXIS_ERR = {5e-6,1.5e-4,4e-6,2e-4,3e-6,2.5e-4};
	auto R_TotalSysErr = TotalSysUncer_Map.at("R_TotalSysUncer");
	for (int i = 0; i < R_TotalSysErr.size(); ++i)
	{
		R_TotalSysErr[i] *= 0.01 * TotalSysUncer_Map.at("R")[i];
	}
	TGraphErrors* ge_CMS        	= new TGraphErrors(ShadowRatio_ParamsMap["Xs"].size(),    &ShadowRatio_ParamsMap["Xs"][0], 
					&ShadowRatio_ParamsMap["R"][0],	&ShadowRatio_ParamsMap["Xs_Err"][0],	&ShadowRatio_ParamsMap["R_Err"][0]);
	TGraphAsymmErrors* gae_CMS 	= new TGraphAsymmErrors(
												ShadowRatio_ParamsMap["Xs"].size(),
												&ShadowRatio_ParamsMap["Xs"][0],        &ShadowRatio_ParamsMap["R"][0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&R_TotalSysErr[0], &R_TotalSysErr[0]);

	TGraphErrors* ge_ALICE_Run2_MidRap	= new TGraphErrors(ALICE_Run2_MidRap_x.size(),	&ALICE_Run2_MidRap_x[0],	&ALICE_Run2_MidRap_R[0],	0,	&ALICE_Run2_MidRap_R_StatErr[0]	);
	TGraphErrors* ge_ALICE_Run2_FwdRap	= new TGraphErrors(ALICE_Run2_FwdRap_x.size(),	&ALICE_Run2_FwdRap_x[0],	&ALICE_Run2_FwdRap_R[0],	0,	&ALICE_Run2_FwdRap_R_StatErr[0]	);
	X_AXIS_ERR = {3e-5};
	TGraphAsymmErrors* gae_ALICE_Run2_MidRap 	= new TGraphAsymmErrors(
												ALICE_Run2_MidRap_x.size(),
												&ALICE_Run2_MidRap_x[0],        &ALICE_Run2_MidRap_R[0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&ALICE_Run2_MidRap_R_SysErrLow[0], &ALICE_Run2_MidRap_R_SysErrHig[0]);
	X_AXIS_ERR = {1.6e-3,1.2e-3};
	TGraphAsymmErrors* gae_ALICE_Run2_FwdRap 	= new TGraphAsymmErrors(
												ALICE_Run2_FwdRap_x.size(),
												&ALICE_Run2_FwdRap_x[0],        &ALICE_Run2_FwdRap_R[0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&ALICE_Run2_FwdRap_R_SysErrLow[0], &ALICE_Run2_FwdRap_R_SysErrHig[0]);
	gae_CMS ->SetMarkerStyle(24);
	gae_CMS ->SetFillColorAlpha(16, 0.7);
	gae_CMS ->SetFillStyle(1001);
	gae_CMS ->Draw("2same");
	ge_CMS->SetMarkerStyle(20);
	ge_CMS->SetMarkerColor(1);
	ge_CMS->SetLineColor(1);
	ge_CMS->SetLineWidth(2);
	ge_CMS->Draw("pezsame");

	gae_ALICE_Run2_MidRap ->SetMarkerStyle(24);
	gae_ALICE_Run2_MidRap ->SetFillColorAlpha(16, 0.7);
	gae_ALICE_Run2_MidRap ->SetFillStyle(1001);
	gae_ALICE_Run2_MidRap ->Draw("2same");
	ge_ALICE_Run2_MidRap->SetMarkerStyle(24);
	ge_ALICE_Run2_MidRap->SetMarkerColor(4);
	ge_ALICE_Run2_MidRap->SetLineColor(4);
	ge_ALICE_Run2_MidRap->SetLineWidth(2);
	ge_ALICE_Run2_MidRap->Draw("pezsame");

	gae_ALICE_Run2_FwdRap ->SetMarkerStyle(24);
	gae_ALICE_Run2_FwdRap ->SetFillColorAlpha(16, 0.7);
	gae_ALICE_Run2_FwdRap ->SetFillStyle(1001);
	gae_ALICE_Run2_FwdRap ->Draw("2same");
	ge_ALICE_Run2_FwdRap->SetMarkerStyle(25);
	ge_ALICE_Run2_FwdRap->SetMarkerColor(4);
	ge_ALICE_Run2_FwdRap->SetLineColor(4);
	ge_ALICE_Run2_FwdRap->SetLineWidth(2);
	ge_ALICE_Run2_FwdRap->Draw("pezsame");

	TLegend  *leg =  new TLegend(0.13, 0.60, 0.40, 0.80);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.050);
	leg->AddEntry(ge_CMS,         "CMS",        "p");
	leg->AddEntry(ge_ALICE_Run2_FwdRap, "ALICE (-4 < y < -3.5)", "p");
	leg->AddEntry(ge_ALICE_Run2_MidRap, "ALICE (|y| < 0.15)",  "p");
	leg->Draw("same");

	drawLatex(0.15, 0.86, "Pb+Pb UPC #sqrt{s_{NN}} = 5.02 TeV",  42,        0.05,      1 );
	drawLatex(0.42, 0.22, "Nuclear suppression factor", 42, 0.06, 1);

	c->SaveAs(outDir+"/ShadowingRatiovsX.png");
	c->SaveAs(outDir+"/ShadowingRatiovsX.pdf");

	delete c;

	return ShadowRatio_ParamsMap;
}

void PlottingForPaper()
{
	auto ShadowRatio_ParamsMap 	= readMap("rootfiles/Results_Map.root");
	auto TotalSysUncer_Map 		= readMap("rootfiles/TotalSysUncer_Map.root");
	plotSigmaVsW(ShadowRatio_ParamsMap,TotalSysUncer_Map);
	plotRvsX(ShadowRatio_ParamsMap,TotalSysUncer_Map);
}
