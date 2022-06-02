#include "plotShadowingRatio.C"


std::map<TString, std::vector<double>> TotalSysUncer_Map = {
	{"Xs",			{}},	{"Ws",					{}},
	{"Sigmas",		{}},	{"R",					{}},
	{"R_TotalSysUncer",	{}},	{"Sigmas_TotalSysUncer",		{}},
	{"Rap",			{}},	
	{"Xsec_AnAn_TotalSysUncer",			{}},	
	{"Xsec_0n0n_TotalSysUncer",			{}},
	{"Xsec_0nXnSum_TotalSysUncer",		{}},
	{"Xsec_XnXn_TotalSysUncer",			{}}
};

void Relative_Uncer( std::map<TString, std::vector<double>> &m1,	std::map<TString, std::vector<double>> m2, const TString param)
{
	std::vector<double> out;
	for (int i = 0; i < m1[param].size(); ++i)
	{
		out.push_back( 100 * abs(m1[param][i] - m2[param][i]) / m2[param][i] );
	}
	m1[param+"_Uncer"] = out;
}

void Cal_R_Uncer( std::map<TString, std::vector<double>> &map , const bool AddIA_Uncer = false)
{
	for (int i = 0; i < map.at("Xs").size(); ++i)
	{
		if (AddIA_Uncer)	map["R_Uncer"].push_back( 0.5 * 100 * map.at("Sigmas_IA_Err")[i] / map.at("Sigmas_IA")[i] );
		else 				map["R_Uncer"].push_back( 0.5 * map.at("Sigmas_Uncer")[i] 	);
	}
}

void CalculateUncer( std::map<TString, std::vector<double>> &Map, const std::map<TString, std::vector<double>> &Map_Default)
{
	//calculate the uncertainty if it has not been calculated
	if ( Map.find("Sigmas_Uncer") == Map.end() ){	Relative_Uncer(Map,	Map_Default,	"Sigmas");	}
	if ( Map.find("R_Uncer") 	  == Map.end() ){	Cal_R_Uncer(Map);					}

	if ( Map.find("Xsec_AnAn_Uncer") == Map.end() ){	Relative_Uncer(Map,	Map_Default,	"Xsec_AnAn");	}
	if ( Map.find("Xsec_0n0n_Uncer") == Map.end() ){	Relative_Uncer(Map,	Map_Default,	"Xsec_0n0n");	}
	if ( Map.find("Xsec_0nXnSum_Uncer") == Map.end() ){	Relative_Uncer(Map,	Map_Default,	"Xsec_0nXnSum");	}
	if ( Map.find("Xsec_XnXn_Uncer") == Map.end() ){	Relative_Uncer(Map,	Map_Default,	"Xsec_XnXn");	}
}

std::map<TString, std::vector<double>> getLargestUncer_Map(const std::vector< std::map<TString, std::vector<double>> > Maps, TString paramX, TString paramY)
{
	//pick the largest uncertainty for each point from the maps
	std::map<TString, std::vector<double>> LargestUncer_Map;
	LargestUncer_Map[paramX]	=	Maps[0].at(paramX);
	LargestUncer_Map[paramY+"_Uncer"]	=	std::vector<double>(LargestUncer_Map[paramX].size(),	0);

	for (int i = 0; i < LargestUncer_Map[paramX].size(); ++i)
	{
		for (int j = 0; j < Maps.size(); ++j)
		{
			if (LargestUncer_Map[paramY+"_Uncer"][i] < Maps[j].at(paramY+"_Uncer")[i])
			{
				LargestUncer_Map[paramY+"_Uncer"][i] = Maps[j].at(paramY+"_Uncer")[i];
			}
		}
	}
	return LargestUncer_Map;
}

std::map<TString, std::vector<double>> getCombinedUncer_Map(const std::vector< std::map<TString, std::vector<double>> > Maps, TString paramX, TString paramY)
{
	std::map<TString, std::vector<double>> CombinedUncer_Map;
	CombinedUncer_Map[paramX]	=	Maps[0].at(paramX);
	CombinedUncer_Map[paramY+"_Uncer"]	=	std::vector<double>(CombinedUncer_Map[paramX].size(),	0);

	for (int i = 0; i < CombinedUncer_Map[paramX].size(); ++i)
	{
		for (int j = 0; j < Maps.size(); ++j)
		{
			CombinedUncer_Map.at(paramY+"_Uncer")[i] = TMath::Hypot(CombinedUncer_Map.at(paramY+"_Uncer")[i],
																	Maps[j].at(paramY+"_Uncer")[i]				);
		}
	}
	return CombinedUncer_Map;
}


void AddToTotalSysUncer(	const std::map<TString, std::vector<double>> &Map, TString paramY)
{
	for (int i = 0; i < TotalSysUncer_Map.at(paramY+"_TotalSysUncer").size(); ++i)
	{
		TotalSysUncer_Map.at(paramY+"_TotalSysUncer")[i]		=	TMath::Hypot(TotalSysUncer_Map.at(paramY+"_TotalSysUncer")[i],		Map.at(paramY+"_Uncer")[i]);
	}
}

void plot_Uncer( std::map<TString, std::vector<double>> &Map, const std::map<TString, std::vector<double>> &Map_Default,
				TString paramX, TString paramY, int style, int color, TLegend *legend, TString legendName,	const bool AddIA_Uncer = false)
{
	if(AddIA_Uncer && paramY == "R")	Cal_R_Uncer(Map,	AddIA_Uncer);

	AddToTotalSysUncer(Map,	paramY);

	auto gr	= new TGraph(Map.at(paramX).size(),		&Map.at(paramX)[0],		&Map.at(paramY+"_Uncer")[0]);
	gr->SetMarkerSize(1.5);
	gr->SetMarkerColor(color);
	gr->SetMarkerStyle(style);
	gr->Draw("psame");

	legend->AddEntry(gr, legendName, "p");
}

void plot_ConstUncer(const std::map<TString, std::vector<double>> &Map,	TString paramX,	TString paramY,TLegend *legend)
{
	std::vector<double> v_Lumi_Uncer;
	std::vector<double> v_BR_Uncer;
	if (paramY == "R")
	{
		v_Lumi_Uncer	=	std::vector<double>	(Map.at(paramX).size(),	0.5*Lumi_Uncer*100);
		v_BR_Uncer		=	std::vector<double> (Map.at(paramX).size(),	0.5*Jpsi2MuMuBR_Uncer*100);
	}
	else{
		v_Lumi_Uncer	=	std::vector<double>	(Map.at(paramX).size(),	Lumi_Uncer*100);
		v_BR_Uncer		=	std::vector<double> (Map.at(paramX).size(),	Jpsi2MuMuBR_Uncer*100);
	}

	AddToTotalSysUncer({{paramY+"_Uncer",	v_Lumi_Uncer}},	paramY);
	AddToTotalSysUncer({{paramY+"_Uncer",	v_BR_Uncer}},	paramY);

	auto gr_Lumi_Uncer		= new TGraph(Map.at(paramX).size(),		&Map.at(paramX)[0],		&v_Lumi_Uncer[0]	);
	auto gr_BR_Uncer		= new TGraph(Map.at(paramX).size(),		&Map.at(paramX)[0],		&v_BR_Uncer[0]		);
	auto gr_TotalSys_Uncer	= new TGraph(TotalSysUncer_Map.at(paramX).size(),		&TotalSysUncer_Map.at(paramX)[0],		&TotalSysUncer_Map.at(paramY+"_TotalSysUncer")[0]		);

	gr_Lumi_Uncer->SetMarkerStyle(27);
	gr_Lumi_Uncer->SetMarkerSize(1.5);
	gr_Lumi_Uncer->SetMarkerColor(1);
	gr_Lumi_Uncer->Draw("psame");

	gr_BR_Uncer->SetMarkerStyle(27);
	gr_BR_Uncer->SetMarkerSize(1.5);
	gr_BR_Uncer->SetMarkerColor(2);
	gr_BR_Uncer->Draw("psame");

	gr_TotalSys_Uncer->SetMarkerStyle(47);
	gr_TotalSys_Uncer->SetMarkerSize(1.8);
	gr_TotalSys_Uncer->SetMarkerColor(1);
	gr_TotalSys_Uncer->Draw("psame");

	legend->AddEntry(gr_BR_Uncer,		"BR", 		"p");
	legend->AddEntry(gr_Lumi_Uncer,		"Lumi", 	"p");
	legend->AddEntry(gr_TotalSys_Uncer,	"Total", 	"p");
}

void plotSigmaVsW( std::map<TString, std::vector<double>> &Map, const double shift,
					 int style, int color, TString legendName, const bool LastCall = false)
{
	static bool FirstCall = true;
	static auto c1 = new TCanvas();
	static TH2D* SigmaVsW = new TH2D("SigmaVsW", ";W_{#gammaPb} (GeV);#sigma(#gamma A #rightarrow J/#psi A) (mb);", 10,0,420, 10, 0.005, 0.20);
	static TLegend  *legData =  new TLegend(0.25, 0.15, 0.4, 0.40);

	if ( FirstCall )
	{
		//gPad->SetPad(0.0,0.25,1.0,0.96);
		gPad->SetTopMargin(0.05);
		gPad->SetBottomMargin(0.12);
		//gPad->SetLeftMargin(0.05);
		gPad->SetRightMargin(0.05);
		c1->SetLogy();

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

		auto X_AXIS_ERR = std::vector<double>(Map.at("Ws").size(), 3);
		
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
		TGraphErrors* ge_IA        	= new TGraphErrors(Map.at("Ws_FitIA").size(),	&Map.at("Ws_FitIA")[0],	&Map.at("Sigmas_FitIA")[0],	0,	0);

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

		drawLatex(0.15, 0.87, "Pb+Pb UPC #sqrt{s_{NN}} = 5.02 TeV",  42,        0.05,      1 );
		
		TLegend  *legData =  new TLegend(0.13, 0.66, 0.35, 0.82);
		legData->SetFillStyle(0);
		legData->SetFillColor(0);
		legData->SetTextSize(0.040);
		legData->AddEntry(ge_ALICE_Run2_FwdRap, "ALICE (-4 < y < -3.5)", "p");
		legData->AddEntry(ge_ALICE_Run2_MidRap, "ALICE (|y| < 0.15)",  "p");
		legData->Draw("same");

		TLegend  *legTheory =  new TLegend(0.50, 0.20, 0.80, 0.30);
		legTheory->SetFillStyle(0);
		legTheory->SetFillColor(0);
		legTheory->SetTextSize(0.040);
		legTheory->AddEntry(ge_CGCnoFluct,        "CGC",                   "l");
		legTheory->AddEntry(ge_IA,                "Impulse Approximation", "l");
		legTheory->Draw("same");
	}

	std::vector<double> Ws_shift;
	for (int i = 0; i < Map["Ws"].size(); ++i){	Ws_shift.push_back(Map["Ws"][i]+shift); }

	TGraphErrors* ge_CMS	= new TGraphErrors(Map["Ws"].size(),	&Ws_shift[0],	&Map["Sigmas"][0],	0,	&Map["Sigmas_Err"][0]	);
	ge_CMS->SetMarkerStyle(style);
	ge_CMS->SetMarkerSize(1.1);
	ge_CMS->SetMarkerColor(color);
	ge_CMS->SetLineColor(color);
	ge_CMS->SetLineWidth(2);
	ge_CMS->Draw("pezsame");

	legData->AddEntry(ge_CMS,	legendName,	"p");
	legData->SetTextSize(0.030);

	FirstCall = false;
	if (LastCall) {
		legData->Draw("same");
		c1->SaveAs("outplots/ComparingSigmaVsW_logY.png");
		c1->SaveAs("outplots/ComparingSigmaVsW_logY.pdf");
		// delete c1;
	}
}

void plotRvsX( std::map<TString, std::vector<double>> &Map, const double shift,
					 int style, int color, TString legendName, const bool LastCall = false )
{	
	static bool FirstCall = true;
	static auto c = new TCanvas();
	static TH2D* htem2d = new TH2D("htem2d", "", 10,4.0e-5,5e-2, 10, 0.2, 1.1);
	static TLegend  *leg =  new TLegend(0.17, 0.60, 0.26, 0.85);

	if (FirstCall)
	{
		c->SetLogx();
		gPad->SetTopMargin(0.05);
		gPad->SetBottomMargin(0.12);
		gPad->SetRightMargin(0.05);

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

		drawLatex(0.15, 0.86, "Pb+Pb UPC #sqrt{s_{NN}} = 5.02 TeV",  42,        0.05,      1 );
		drawLatex(0.42, 0.22, "Nuclear suppression factor", 42, 0.06, 1);

		std::vector<double> X_AXIS_ERR = {5e-6,1.5e-4,4e-6,2e-4,3e-6,2.5e-4};
		auto R_TotalSysErr = TotalSysUncer_Map.at("R_TotalSysUncer");
		for (int i = 0; i < R_TotalSysErr.size(); ++i)
		{
			R_TotalSysErr[i] *= 0.01 * TotalSysUncer_Map.at("R")[i];
		}
		
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
		leg->AddEntry(ge_ALICE_Run2_FwdRap, "ALICE (-4 < y < -3.5)", "p");
		leg->AddEntry(ge_ALICE_Run2_MidRap, "ALICE (|y| < 0.15)",  "p");
	}
	
	std::vector<double> Xs_shift;
	for (int i = 0; i < Map["Xs"].size(); ++i)
	{	
		if 	(Map["Xs"][i] > 1e-3)	{ Xs_shift.push_back(Map["Xs"][i]+shift); 		}
		else 						{ Xs_shift.push_back(Map["Xs"][i]+shift/50); 	}
	}

	TGraphErrors* ge_CMS	= new TGraphErrors(Map["Xs"].size(), &Xs_shift[0], &Map["R"][0], &Map["Xs_Err"][0], &Map["R_Err"][0]);

	ge_CMS->SetMarkerStyle(style);
	ge_CMS->SetMarkerColor(color);
	ge_CMS->SetLineColor(color);
	ge_CMS->SetLineWidth(2);
	ge_CMS->Draw("pezsame");
	
	leg->AddEntry(ge_CMS,         legendName,        "lp");
	leg->SetTextSize(0.030);

	FirstCall = false;
	if (LastCall)
	{
		leg->Draw("same");
		c->SaveAs("outplots/ComparingRvsX.png");
		c->SaveAs("outplots/ComparingRvsX.pdf");
		// delete c;
	}
}


void getFitUncertainty()
{
	std::map<TString, std::vector<double>> CB_Poly3_Map 		= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP.root");
	std::map<TString, std::vector<double>> CB_Poly4_Map 		= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly4_PUShuai_6RapBins.appliedTnP.root");
	std::map<TString, std::vector<double>> CB_FixCBAN_Poly3_Map = getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_FixCBAN_Poly3_PUShuai_6RapBins.appliedTnP.root");
	std::map<TString, std::vector<double>> CBG_Poly3_Map 		= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CBG_Poly3_PUShuai_6RapBins.appliedTnP.root");
	std::map<TString, std::vector<double>> CB_Poly3_SdB_Map 	= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_SdB_PUShuai_6RapBins.appliedTnP.root");
	std::map<TString, std::vector<double>> CB_Poly3_NarrMass_Map= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_NarrowMass_6RapBins.appliedTnP.root");
	std::map<TString, std::vector<double>> CB_Poly3_WideMass_Map= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_WideMass_6RapBins.appliedTnP.root");
	std::map<TString, std::vector<double>> CB_Poly3_looseHF_Map = getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.looseHF.root");
	// std::map<TString, std::vector<double>> CB_Poly3_tightHF_Map = getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.tightHF.root");
	std::map<TString, std::vector<double>> CB_Poly3_fluxP_Map 	= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP.root", 1.05);
	std::map<TString, std::vector<double>> CB_Poly3_fluxM_Map 	= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP.root", 0.95);
	std::map<TString, std::vector<double>> CB_Poly3_TnP_Low_Map = getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP_Low.root");
	std::map<TString, std::vector<double>> CB_Poly3_TnP_Hig_Map	= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP_Hig.root");
	std::map<TString, std::vector<double>> CB_Poly3_PU_Map 		= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUATLAS_6RapBins.appliedTnP.root");

	const std::map<TString, std::vector<double>> JpsiXsec_Default		= CB_Poly3_Map;

	CalculateUncer(CB_Poly4_Map,			JpsiXsec_Default);		CalculateUncer(CB_FixCBAN_Poly3_Map,	JpsiXsec_Default);
	CalculateUncer(CBG_Poly3_Map,			JpsiXsec_Default);		CalculateUncer(CB_Poly3_SdB_Map,		JpsiXsec_Default);
	CalculateUncer(CB_Poly3_NarrMass_Map,	JpsiXsec_Default);		CalculateUncer(CB_Poly3_WideMass_Map,	JpsiXsec_Default);
	CalculateUncer(CB_Poly3_looseHF_Map,	JpsiXsec_Default);		CalculateUncer(CB_Poly3_fluxP_Map,		JpsiXsec_Default);
	CalculateUncer(CB_Poly3_fluxM_Map,		JpsiXsec_Default);		CalculateUncer(CB_Poly3_TnP_Low_Map,	JpsiXsec_Default);
	CalculateUncer(CB_Poly3_TnP_Hig_Map,	JpsiXsec_Default);		CalculateUncer(CB_Poly3_PU_Map,			JpsiXsec_Default);

	TotalSysUncer_Map["Ws"] = JpsiXsec_Default.at("Ws");	TotalSysUncer_Map["Xs"] = JpsiXsec_Default.at("Xs");
	TotalSysUncer_Map["Sigmas"] 	= JpsiXsec_Default.at("Sigmas");
	TotalSysUncer_Map["R"] 			= JpsiXsec_Default.at("R");
	TotalSysUncer_Map["Sigmas_TotalSysUncer"] 	= std::vector<double>(JpsiXsec_Default.at("Ws").size(),	0);
	TotalSysUncer_Map["R_TotalSysUncer"] 		= std::vector<double>(JpsiXsec_Default.at("Xs").size(),	0);	
	TotalSysUncer_Map["Rap"] 		= JpsiXsec_Default.at("Rap");
	TotalSysUncer_Map["Xsec_AnAn_TotalSysUncer"] 		= std::vector<double>(JpsiXsec_Default.at("Rap").size(),	0);	
	TotalSysUncer_Map["Xsec_0n0n_TotalSysUncer"] 		= std::vector<double>(JpsiXsec_Default.at("Rap").size(),	0);	
	TotalSysUncer_Map["Xsec_0nXnSum_TotalSysUncer"] 	= std::vector<double>(JpsiXsec_Default.at("Rap").size(),	0);	
	TotalSysUncer_Map["Xsec_XnXn_TotalSysUncer"] 		= std::vector<double>(JpsiXsec_Default.at("Rap").size(),	0);	

	//-----------------------------------------------------------------------------------
	//-----------------------------------Comparing SigmaVsW------------------------------
	plotSigmaVsW(CB_Poly3_Map,			0,	47,	1,	"Default"				);
	plotSigmaVsW(CB_Poly4_Map,			0,	24,	1,	"CB_Poly4"				);
	plotSigmaVsW(CB_FixCBAN_Poly3_Map,	0,	24,	2,	"CB_FixCBAN_Poly3"		);
	plotSigmaVsW(CBG_Poly3_Map,			0,	24,	6,	"CBG_Poly3"				);
	plotSigmaVsW(CB_Poly3_SdB_Map,		-5,	24,	28,	"QEDPtShape"			);
	plotSigmaVsW(CB_Poly3_NarrMass_Map,	-5,	24,	29,	"NarrMass"				);
	plotSigmaVsW(CB_Poly3_WideMass_Map,	-5,	24,	30,	"WideMass"				);

	plotSigmaVsW(CB_Poly3_looseHF_Map,	5,	26,	1,	"HFveto"		);
	// plotSigmaVsW(CB_Poly3_fluxP_Map,	10,	26,	1,	"CB_Poly3_fluxP"		);
	plotSigmaVsW(CB_Poly3_fluxM_Map,	5,	26,	2,	"PhotonFlux" 			);
	plotSigmaVsW(CB_Poly3_PU_Map,		5,	26,	4,	"n-PileUp"				);
	plotSigmaVsW(CB_Poly3_TnP_Low_Map,	5,	26,	6,	"TnP", 				true);
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//-----------------------------------Comparing RVsX----------------------------------
	plotRvsX(CB_Poly3_Map,			0,	47,	1,	"Default"				);
	plotRvsX(CB_Poly4_Map,			0,	24,	1,	"CB_Poly4"				);
	plotRvsX(CB_FixCBAN_Poly3_Map,	0,	24,	2,	"CB_FixCBAN_Poly3"		);
	plotRvsX(CBG_Poly3_Map,			0,	24,	6,	"CBG_Poly3"				);
	plotRvsX(CB_Poly3_SdB_Map,		-4e-4,	24,	28,	"QEDPtShape"			);
	plotRvsX(CB_Poly3_NarrMass_Map,	-4e-4,	24,	29,	"NarrMass"				);
	plotRvsX(CB_Poly3_WideMass_Map,	-4e-4,	24,	30,	"WideMass"				);

	plotRvsX(CB_Poly3_looseHF_Map,	4e-4,	26,	1,	"HFveto"			);
	// plotRvsX(CB_Poly3_fluxP_Map,	0,	26,	1,	"CB_Poly3_fluxP"		);
	plotRvsX(CB_Poly3_fluxM_Map,	4e-4,	26,	2,	"PhotonFlux"		);
	plotRvsX(CB_Poly3_PU_Map,		4e-4,	26,	4,	"n-PileUp"			);
	plotRvsX(CB_Poly3_TnP_Low_Map,	4e-4,	26,	6,	"TnP", 			true);
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//-----------------------------------W Unicertainty----------------------------------
	auto c2 = new TCanvas();
	c2->SetLogx();
	auto Sigma_FitUncer	= new TH2D("Sigma_FitUncertainty", ";W_{#gammaPb} (GeV);#sigma Uncer.(%);", 10,30,520, 10,0,15);
	Sigma_FitUncer->Draw();

	TLegend  *legendW =  new TLegend(0.35, 0.60, 0.50, 0.90);
	legendW->SetFillStyle(0);
	legendW->SetTextSize(0.035);

	auto SigmasMassFitRangeUncer_Map	= getLargestUncer_Map( { CB_Poly3_NarrMass_Map,	CB_Poly3_WideMass_Map},	"Ws",	"Sigmas");
	auto SigmasTotalFitUncer_Map 		= getCombinedUncer_Map( { CB_Poly4_Map,		CB_FixCBAN_Poly3_Map,	CBG_Poly3_Map,
																CB_Poly3_SdB_Map,	SigmasMassFitRangeUncer_Map},"Ws",	"Sigmas");
	auto SigmasFluxUncer_Map			= getLargestUncer_Map( { CB_Poly3_fluxM_Map,	CB_Poly3_fluxP_Map},	"Ws",	"Sigmas" );
	auto SigmasTnPUncer_Map				= getLargestUncer_Map( { CB_Poly3_TnP_Low_Map,	CB_Poly3_TnP_Hig_Map},	"Ws",	"Sigmas" );

	plot_Uncer(SigmasTotalFitUncer_Map,	JpsiXsec_Default,	"Ws",	"Sigmas",	24,	1,	legendW,	"Signal Ext");
	plot_Uncer(SigmasFluxUncer_Map,		JpsiXsec_Default,	"Ws",	"Sigmas",	25,	1,	legendW,	"PhotonFlux");
	plot_Uncer(CB_Poly3_looseHF_Map,	JpsiXsec_Default,	"Ws",	"Sigmas",	26,	1,	legendW,	"HFveto");
	plot_Uncer(CB_Poly3_PU_Map,			JpsiXsec_Default,	"Ws",	"Sigmas",	28,	1,	legendW,	"n-PileUp");
	plot_Uncer(SigmasTnPUncer_Map,		JpsiXsec_Default,	"Ws",	"Sigmas",	32,	1,	legendW,	"TnP");

	plot_ConstUncer(JpsiXsec_Default,	"Ws",	"Sigmas",	legendW);
	legendW->Draw("same");

	c2->SaveAs("outplots/Uncertainty_SigmaVsW.png");
	c2->SaveAs("outplots/Uncertainty_SigmaVsW.pdf");
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//-----------------------------------R Unicertainty----------------------------------
	auto c3 = new TCanvas();
	c3->SetLogx();
	auto R_FitUncer	= new TH2D("R_FitUncertainty", ";x;R^{Pb}_{g} Uncer.(%);", 10,4.0e-5,1.1e-2, 10, 0, 15);
	R_FitUncer->Draw();

	TLegend  *legendR =  new TLegend(0.35, 0.60, 0.50, 0.90);
	legendR->SetFillStyle(0);
	legendR->SetTextSize(0.035);

	auto RMassFitRangeUncer_Map		= getLargestUncer_Map( { CB_Poly3_NarrMass_Map,	CB_Poly3_WideMass_Map},	"Xs",	"R");
	auto RTotalFitUncer_Map 		= getCombinedUncer_Map( { CB_Poly4_Map,		CB_FixCBAN_Poly3_Map,	CBG_Poly3_Map,
																CB_Poly3_SdB_Map,	RMassFitRangeUncer_Map},"Xs",	"R");
	auto RFluxUncer_Map				= getLargestUncer_Map( { CB_Poly3_fluxM_Map,	CB_Poly3_fluxP_Map},	"Xs",	"R" );
	auto RTnPUncer_Map				= getLargestUncer_Map( { CB_Poly3_TnP_Low_Map,	CB_Poly3_TnP_Hig_Map},	"Xs",	"R" );

	plot_Uncer(RTotalFitUncer_Map,	JpsiXsec_Default,	"Xs",	"R",	24,	1,	legendR,	"Signal Ext");
	plot_Uncer(RFluxUncer_Map,		JpsiXsec_Default,	"Xs",	"R",	25,	1,	legendR,	"PhotonFlux");
	plot_Uncer(CB_Poly3_looseHF_Map,JpsiXsec_Default,	"Xs",	"R",	26,	1,	legendR,	"HFveto");
	plot_Uncer(CB_Poly3_PU_Map,		JpsiXsec_Default,	"Xs",	"R",	28,	1,	legendR,	"n-PileUp");
	plot_Uncer(RTnPUncer_Map,		JpsiXsec_Default,	"Xs",	"R",	32,	1,	legendR,	"TnP");

	plot_Uncer(CB_Poly3_Map,		JpsiXsec_Default,	"Xs",	"R",	27,	4,	legendR,	"IA", 	true);
	plot_ConstUncer(JpsiXsec_Default,	"Xs",	"R",	legendR);
	legendR->Draw("same");

	c3->SaveAs("outplots/Uncertainty_RVsX.png");
	c3->SaveAs("outplots/Uncertainty_RVsX.pdf");
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//-----------------------------------Xsec_AnAn Unicertainty--------------------------
	auto c4 = new TCanvas();
	auto hXsec_AnAn_Uncer	= new TH2D("hXsec_AnAn_Uncer", "hXsec_AnAn_Uncer;y;d#sigma_{J/#psi}/dy (mb);", 10, -3, -1, 10, 0, 20);
	hXsec_AnAn_Uncer->Draw();

	TLegend  *legend_AnAn =  new TLegend(0.35, 0.60, 0.50, 0.90);
	legend_AnAn->SetFillStyle(0);
	legend_AnAn->SetTextSize(0.035);

	auto Xsec_AnAnMassFitRangeUncer_Map		= getLargestUncer_Map( { CB_Poly3_NarrMass_Map,	CB_Poly3_WideMass_Map},	"Rap",	"Xsec_AnAn");
	auto Xsec_AnAnTotalFitUncer_Map 		= getCombinedUncer_Map( { CB_Poly4_Map,		CB_FixCBAN_Poly3_Map,	CBG_Poly3_Map,
																CB_Poly3_SdB_Map,	Xsec_AnAnMassFitRangeUncer_Map},"Rap",	"Xsec_AnAn");
	auto Xsec_AnAnTnPUncer_Map				= getLargestUncer_Map( { CB_Poly3_TnP_Low_Map,	CB_Poly3_TnP_Hig_Map},	"Rap",	"Xsec_AnAn" );

	plot_Uncer(Xsec_AnAnTotalFitUncer_Map,	JpsiXsec_Default,	"Rap",	"Xsec_AnAn",	24,	1,	legend_AnAn,	"Signal Ext");
	plot_Uncer(CB_Poly3_looseHF_Map,		JpsiXsec_Default,	"Rap",	"Xsec_AnAn",	26,	1,	legend_AnAn,	"HFveto");
	// plot_Uncer(CB_Poly3_PU_Map,				JpsiXsec_Default,	"Rap",	"Xsec_AnAn",	28,	1,	legend_AnAn,	"n-PileUp");
	plot_Uncer(Xsec_AnAnTnPUncer_Map,		JpsiXsec_Default,	"Rap",	"Xsec_AnAn",	32,	1,	legend_AnAn,	"TnP");

	plot_ConstUncer(JpsiXsec_Default,	"Rap",	"Xsec_AnAn",	legend_AnAn);
	legend_AnAn->Draw("same");

	c4->SaveAs("outplots/Uncertainty_Xsec_AnAn.png");
	c4->SaveAs("outplots/Uncertainty_Xsec_AnAn.pdf");
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//-----------------------------------Xsec_0n0n Unicertainty--------------------------
	auto c5 = new TCanvas();
	auto hXsec_0n0n_Uncer	= new TH2D("hXsec_0n0n_Uncer", "hXsec_0n0n_Uncer;y;d#sigma_{J/#psi}/dy (mb);", 10, -3, -1, 10, 0, 20);
	hXsec_0n0n_Uncer->Draw();

	TLegend  *legend_0n0n =  new TLegend(0.35, 0.60, 0.50, 0.90);
	legend_0n0n->SetFillStyle(0);
	legend_0n0n->SetTextSize(0.035);

	auto Xsec_0n0nMassFitRangeUncer_Map		= getLargestUncer_Map( { CB_Poly3_NarrMass_Map,	CB_Poly3_WideMass_Map},	"Rap",	"Xsec_0n0n");
	auto Xsec_0n0nTotalFitUncer_Map 		= getCombinedUncer_Map( { CB_Poly4_Map,		CB_FixCBAN_Poly3_Map,	CBG_Poly3_Map,
																CB_Poly3_SdB_Map,	Xsec_0n0nMassFitRangeUncer_Map},"Rap",	"Xsec_0n0n");
	auto Xsec_0n0nTnPUncer_Map				= getLargestUncer_Map( { CB_Poly3_TnP_Low_Map,	CB_Poly3_TnP_Hig_Map},	"Rap",	"Xsec_0n0n" );

	plot_Uncer(Xsec_0n0nTotalFitUncer_Map,	JpsiXsec_Default,	"Rap",	"Xsec_0n0n",	24,	1,	legend_0n0n,	"Signal Ext");
	plot_Uncer(CB_Poly3_looseHF_Map,		JpsiXsec_Default,	"Rap",	"Xsec_0n0n",	26,	1,	legend_0n0n,	"HFveto");
	plot_Uncer(CB_Poly3_PU_Map,				JpsiXsec_Default,	"Rap",	"Xsec_0n0n",	28,	1,	legend_0n0n,	"n-PileUp");
	plot_Uncer(Xsec_0n0nTnPUncer_Map,		JpsiXsec_Default,	"Rap",	"Xsec_0n0n",	32,	1,	legend_0n0n,	"TnP");

	plot_ConstUncer(JpsiXsec_Default,	"Rap",	"Xsec_0n0n",	legend_0n0n);
	legend_0n0n->Draw("same");

	c5->SaveAs("outplots/Uncertainty_Xsec_0n0n.png");
	c5->SaveAs("outplots/Uncertainty_Xsec_0n0n.pdf");
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//-----------------------------------Xsec_0nXnSum Unicertainty-----------------------
	auto c6 = new TCanvas();
	auto hXsec_0nXnSum_Uncer	= new TH2D("hXsec_0nXnSum_Uncer", "hXsec_0nXnSum_Uncer;y;d#sigma_{J/#psi}/dy (mb);", 10, -3, -1, 10, 0, 20);
	hXsec_0nXnSum_Uncer->Draw();

	TLegend  *legend_0nXnSum =  new TLegend(0.35, 0.60, 0.50, 0.90);
	legend_0nXnSum->SetFillStyle(0);
	legend_0nXnSum->SetTextSize(0.035);

	auto Xsec_0nXnSumMassFitRangeUncer_Map		= getLargestUncer_Map( { CB_Poly3_NarrMass_Map,	CB_Poly3_WideMass_Map},	"Rap",	"Xsec_0nXnSum");
	auto Xsec_0nXnSumTotalFitUncer_Map 		= getCombinedUncer_Map( { CB_Poly4_Map,		CB_FixCBAN_Poly3_Map,	CBG_Poly3_Map,
																CB_Poly3_SdB_Map,	Xsec_0nXnSumMassFitRangeUncer_Map},"Rap",	"Xsec_0nXnSum");
	auto Xsec_0nXnSumTnPUncer_Map				= getLargestUncer_Map( { CB_Poly3_TnP_Low_Map,	CB_Poly3_TnP_Hig_Map},	"Rap",	"Xsec_0nXnSum" );

	plot_Uncer(Xsec_0nXnSumTotalFitUncer_Map,	JpsiXsec_Default,	"Rap",	"Xsec_0nXnSum",	24,	1,	legend_0nXnSum,	"Signal Ext");
	plot_Uncer(CB_Poly3_looseHF_Map,			JpsiXsec_Default,	"Rap",	"Xsec_0nXnSum",	26,	1,	legend_0nXnSum,	"HFveto");
	plot_Uncer(CB_Poly3_PU_Map,					JpsiXsec_Default,	"Rap",	"Xsec_0nXnSum",	28,	1,	legend_0nXnSum,	"n-PileUp");
	plot_Uncer(Xsec_0nXnSumTnPUncer_Map,		JpsiXsec_Default,	"Rap",	"Xsec_0nXnSum",	32,	1,	legend_0nXnSum,	"TnP");

	plot_ConstUncer(JpsiXsec_Default,	"Rap",	"Xsec_0nXnSum",	legend_0nXnSum);
	legend_0nXnSum->Draw("same");

	c6->SaveAs("outplots/Uncertainty_Xsec_0nXnSum.png");
	c6->SaveAs("outplots/Uncertainty_Xsec_0nXnSum.pdf");
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//-----------------------------------Xsec_XnXn Unicertainty-----------------------
	auto c7 = new TCanvas();
	auto hXsec_XnXn_Uncer	= new TH2D("hXsec_XnXn_Uncer", "hXsec_XnXn_Uncer;y;d#sigma_{J/#psi}/dy (mb);", 10, -3, -1, 10, 0, 20);
	hXsec_XnXn_Uncer->Draw();

	TLegend  *legend_XnXn =  new TLegend(0.35, 0.60, 0.50, 0.90);
	legend_XnXn->SetFillStyle(0);
	legend_XnXn->SetTextSize(0.035);

	auto Xsec_XnXnMassFitRangeUncer_Map		= getLargestUncer_Map( { CB_Poly3_NarrMass_Map,	CB_Poly3_WideMass_Map},	"Rap",	"Xsec_XnXn");
	auto Xsec_XnXnTotalFitUncer_Map 		= getCombinedUncer_Map( { CB_Poly4_Map,		CB_FixCBAN_Poly3_Map,	CBG_Poly3_Map,
																CB_Poly3_SdB_Map,	Xsec_XnXnMassFitRangeUncer_Map},"Rap",	"Xsec_XnXn");
	auto Xsec_XnXnTnPUncer_Map				= getLargestUncer_Map( { CB_Poly3_TnP_Low_Map,	CB_Poly3_TnP_Hig_Map},	"Rap",	"Xsec_XnXn" );

	plot_Uncer(Xsec_XnXnTotalFitUncer_Map,	JpsiXsec_Default,	"Rap",	"Xsec_XnXn",	24,	1,	legend_XnXn,	"Signal Ext");
	plot_Uncer(CB_Poly3_looseHF_Map,		JpsiXsec_Default,	"Rap",	"Xsec_XnXn",	26,	1,	legend_XnXn,	"HFveto");
	plot_Uncer(CB_Poly3_PU_Map,				JpsiXsec_Default,	"Rap",	"Xsec_XnXn",	28,	1,	legend_XnXn,	"n-PileUp");
	plot_Uncer(Xsec_XnXnTnPUncer_Map,		JpsiXsec_Default,	"Rap",	"Xsec_XnXn",	32,	1,	legend_XnXn,	"TnP");

	plot_ConstUncer(JpsiXsec_Default,	"Rap",	"Xsec_XnXn",	legend_XnXn);
	legend_XnXn->Draw("same");

	c7->SaveAs("outplots/Uncertainty_Xsec_XnXn.png");
	c7->SaveAs("outplots/Uncertainty_Xsec_XnXn.pdf");
	//-----------------------------------------------------------------------------------


	for (int i = 0; i < TotalSysUncer_Map.at("Ws").size(); ++i)
	{
		cout<<Form("W:%f,	dSigma:%f,	X:%f,	dR:%f", TotalSysUncer_Map.at("Ws")[i],	TotalSysUncer_Map.at("Sigmas_TotalSysUncer")[i],
													TotalSysUncer_Map.at("Xs")[i],	TotalSysUncer_Map.at("R_TotalSysUncer")[i])<<endl;
	}
	for (int i = 0; i < TotalSysUncer_Map.at("Rap").size(); ++i)
	{
		cout<<Form("Rap:%f,	AnAn dSigmady:%f,	dSigmady_Uncer:%f", TotalSysUncer_Map.at("Rap")[i],	JpsiXsec_Default.at("Xsec_AnAn")[i],	TotalSysUncer_Map.at("Xsec_AnAn_TotalSysUncer")[i])<<endl;
	}
	for (int i = 0; i < TotalSysUncer_Map.at("Rap").size(); ++i)
	{
		cout<<Form("Rap:%f,	0n0n dSigmady:%f,	dSigmady_Uncer:%f", TotalSysUncer_Map.at("Rap")[i],	JpsiXsec_Default.at("Xsec_0n0n")[i],	TotalSysUncer_Map.at("Xsec_0n0n_TotalSysUncer")[i])<<endl;
	}
	for (int i = 0; i < TotalSysUncer_Map.at("Rap").size(); ++i)
	{
		cout<<Form("Rap:%f,	0nXnSum dSigmady:%f,	dSigmady_Uncer:%f", TotalSysUncer_Map.at("Rap")[i],	JpsiXsec_Default.at("Xsec_0nXnSum")[i],	TotalSysUncer_Map.at("Xsec_0nXnSum_TotalSysUncer")[i])<<endl;
	}

	saveMap(TotalSysUncer_Map,	"rootfiles/TotalSysUncer_Map.root");
}