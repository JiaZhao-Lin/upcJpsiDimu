#include "plotShadowingRatio.C"


std::map<TString, std::vector<double>> TotalSysUncer_Map = {
	{"Xs",			{}},	{"Ws",					{}},
	{"Sigmas",		{}},	{"R",					{}},
	{"R_TotalSysUncer",	{}},	{"Sigmas_TotalSysUncer",		{}}
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
	for (int i = 0; i < map["Xs"].size(); ++i)
	{
		if (AddIA_Uncer)	map["R_Uncer"].push_back( 0.5 * 100 * map["Sigmas_IA_Err"][i] / map["Sigmas_IA"][i] );
		else 				map["R_Uncer"].push_back( 0.5 * map["Sigmas_Uncer"][i] 	);
	}
}

void CalculateUncer( std::map<TString, std::vector<double>> &Map, const std::map<TString, std::vector<double>> &Map_Default)
{
	//calculate the uncertainty if it has not been calculated
	if ( Map.find("Sigmas_Uncer") == Map.end() ){	Relative_Uncer(Map,	Map_Default,	"Sigmas");	}
	if ( Map.find("R_Uncer") 	  == Map.end() ){	Cal_R_Uncer(Map);					}
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
	for (int i = 0; i < Map.at(paramY+"_Uncer").size(); ++i)
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
	std::vector<double> v_Lumi_Uncer	(Map.at(paramX).size(),	Lumi_Uncer*100);
	std::vector<double> v_BR_Uncer		(Map.at(paramX).size(),	Jpsi2MuMuBR_Uncer*100);

	AddToTotalSysUncer({{"Sigmas_Uncer",	v_Lumi_Uncer}},	"Sigmas");	AddToTotalSysUncer({{"R_Uncer",	v_Lumi_Uncer}},	"R");
	AddToTotalSysUncer({{"Sigmas_Uncer",	v_BR_Uncer}},	"Sigmas");	AddToTotalSysUncer({{"R_Uncer",	v_BR_Uncer}},	"R");

	auto gr_Lumi_Uncer		= new TGraph(Map.at(paramX).size(),		&Map.at(paramX)[0],		&v_Lumi_Uncer[0]	);
	auto gr_BR_Uncer		= new TGraph(Map.at(paramX).size(),		&Map.at(paramX)[0],		&v_BR_Uncer[0]		);
	auto gr_TotalSys_Uncer	= new TGraph(TotalSysUncer_Map.at(paramX).size(),		&TotalSysUncer_Map.at(paramX)[0],		&TotalSysUncer_Map.at(paramY+"_TotalSysUncer")[0]		);

	gr_Lumi_Uncer->SetMarkerStyle(56);
	gr_Lumi_Uncer->SetMarkerSize(1.5);
	gr_Lumi_Uncer->SetMarkerColor(1);
	gr_Lumi_Uncer->Draw("psame");

	gr_BR_Uncer->SetMarkerStyle(56);
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
	static TLegend  *legData =  new TLegend(0.30, 0.17, 0.55, 0.40);

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
			ALICE_Run2_FwdRap_W        .push_back( y2W(ALICE_Run2_FwdRap_y[i]) );
			ALICE_Run2_FwdRap_Sigma    .push_back( ALICE_Run2_FwdRap_dSigmady[i] / ALICE_Run2_FwdRap_Flux[i] );
			ALICE_Run2_FwdRap_Sigma_Err.push_back( ALICE_Run2_FwdRap_Sigma[i] * ALICE_Run2_FwdRap_dSigmady_Err[i] / ALICE_Run2_FwdRap_dSigmady[i] );
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

		TGraphAsymmErrors *gae_Xsec_ALICE_Run1 = new TGraphAsymmErrors(ALICE_W_Run1.size(), &ALICE_W_Run1[0],&ALICE_Xsec_Run1[0],0,0, &ALICE_XsecErr2_Run1[0], &ALICE_XsecErr1_Run1[0]);
		
		TGraphErrors* ge_ALICE_Run2_MidRap	= new TGraphErrors(ALICE_Run2_MidRap_W.size(),	&ALICE_Run2_MidRap_W[0],	&ALICE_Run2_MidRap_Sigma[0],	0,	&ALICE_Run2_MidRap_Sigma_Err[0]	);
		TGraphErrors* ge_ALICE_Run2_FwdRap	= new TGraphErrors(ALICE_Run2_FwdRap_W.size(),	&ALICE_Run2_FwdRap_W[0],	&ALICE_Run2_FwdRap_Sigma[0],	0,	&ALICE_Run2_FwdRap_Sigma_Err[0]	);

		TGraphErrors* ge_CGCnoFluct = new TGraphErrors(CGC_JpsiNoFluct_W.size(),	&CGC_JpsiNoFluct_W[0],	&CGC_JpsiNoFluct_CohXsec[0],	0,	0);
		TGraphErrors* ge_IA        	= new TGraphErrors(Map["Ws_FitIA"].size(),	&Map["Ws_FitIA"][0],	&Map["Sigmas_FitIA"][0],	0,	0);
		
		ge_ALICE_Run2_MidRap->SetMarkerStyle(24);
		ge_ALICE_Run2_MidRap->SetMarkerColor(4);
		ge_ALICE_Run2_MidRap->SetLineColor(4);
		ge_ALICE_Run2_MidRap->SetLineWidth(2);
		ge_ALICE_Run2_MidRap->Draw("pezsame");

		ge_ALICE_Run2_FwdRap->SetMarkerStyle(25);
		ge_ALICE_Run2_FwdRap->SetMarkerColor(4);
		ge_ALICE_Run2_FwdRap->SetLineColor(4);
		ge_ALICE_Run2_FwdRap->SetLineWidth(2);
		ge_ALICE_Run2_FwdRap->Draw("pezsame");
		
		// gae_Xsec_ALICE_Run1->SetMarkerStyle(24);
		// gae_Xsec_ALICE_Run1->SetMarkerColor(1);
		// gae_Xsec_ALICE_Run1->SetLineColor(1);
		// gae_Xsec_ALICE_Run1->SetLineWidth(2);
		// gae_Xsec_ALICE_Run1->Draw("Pesame");

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

		drawLatex(0.15, 0.87, "Pb+Pb UPC #sqrt{s_{NN}} = 5.02 TeV",  42,        0.06,      1 );

		TLegend  *legTheory =  new TLegend(0.75, 0.20, 0.95, 0.30);
		legTheory->SetFillStyle(0);
		legTheory->SetFillColor(0);
		legTheory->SetTextSize(0.04);
		legTheory->AddEntry(ge_CGCnoFluct,        "CGC",                   "l");
		legTheory->AddEntry(ge_IA,                "IA", "l");
		legTheory->Draw("same");

		legData->SetFillStyle(0);
		legData->SetFillColor(0);
		legData->SetTextSize(0.035);
		legData->AddEntry(ge_ALICE_Run2_FwdRap, "ALICE (-4 < y < -3.5)", "p");
		legData->AddEntry(ge_ALICE_Run2_MidRap, "ALICE (|y| < 0.15)",    "p");
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
	
	FirstCall = false;
	if (LastCall) {
		legData->Draw("same");
		// c1->SaveAs("outplots/SigmaVsW_logY.png");
		// c1->SaveAs("outplots/SigmaVsW_logY.pdf");
		// delete c1;
	}
}

void plotRvsX( std::map<TString, std::vector<double>> &Map, const double shift,
					 int style, int color, TString legendName, const bool LastCall = false )
{	
	static bool FirstCall = true;
	static auto c = new TCanvas();
	static TH2D* htem2d = new TH2D("htem2d", "", 10,4.0e-5,1.1e-2, 10, 0.2, 1.1);
	static TLegend  *leg =  new TLegend(0.17, 0.60, 0.55, 0.85);

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

		TGraphErrors* ge_ALICE_Run1 	= new TGraphErrors(ALICE_x.size(),      &ALICE_x[0],      &ALICE_R[0],      0,          &ALICE_R_Err[0] );
		TGraphErrors* ge_ALICE_Run2 	= new TGraphErrors(ALICE_Run2_x.size(), &ALICE_Run2_x[0], &ALICE_Run2_R[0], 0,          &ALICE_Run2_R_Err[0] );
		TGraphErrors* ge_ALICE_Run2_Psi = new TGraphErrors(ALICE_Run2_Psi_x.size(), &ALICE_Run2_Psi_x[0], &ALICE_Run2_Psi_R[0], 0,          &ALICE_Run2_R_Err[0] );

		ge_ALICE_Run1->SetMarkerStyle(24);
		ge_ALICE_Run1->SetMarkerColor(4);
		ge_ALICE_Run1->SetLineColor(4);
		ge_ALICE_Run1->SetLineWidth(2);
		//ge_ALICE_Run1->Draw("pesame");

		ge_ALICE_Run2->SetMarkerStyle(24);
		ge_ALICE_Run2->SetMarkerColor(4);
		ge_ALICE_Run2->SetLineColor(4);
		ge_ALICE_Run2->SetLineWidth(2);
		ge_ALICE_Run2->Draw("pezsame");

		ge_ALICE_Run2_Psi->SetMarkerStyle(24);
		ge_ALICE_Run2_Psi->SetMarkerColor(2);
		ge_ALICE_Run2_Psi->SetLineColor(2);
		ge_ALICE_Run2_Psi->SetLineWidth(2);
		// ge_ALICE_Run2_Psi->Draw("pesame");

		leg->SetFillStyle(0);
		leg->SetFillColor(0);
		leg->SetTextSize(0.035);

		//leg->AddEntry(ge_ALICE_Run1,  "ALICE Run1", "lp"); //Run1
		leg->AddEntry(ge_ALICE_Run2,  "ALICE", "lp"); //Run2
		// leg->AddEntry(ge_ALICE_Run2_Psi,  "ALICE Run2 Psi", "lp");
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
	
	FirstCall = false;
	if (LastCall)
	{
		leg->Draw("same");
		// c->SaveAs("outplots/ShadowingRatiovsX.png");
		// c->SaveAs("outplots/ShadowingRatiovsX.pdf");
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

	//-----------------------------------------------------------------------------------
	//-----------------------------------Comparing SigmaVsW------------------------------
	plotSigmaVsW(CB_Poly3_Map,			0,	47,	1,	"Default"				);
	plotSigmaVsW(CB_Poly4_Map,			0,	24,	1,	"CB_Poly4"				);
	plotSigmaVsW(CB_FixCBAN_Poly3_Map,	0,	24,	2,	"CB_FixCBAN_Poly3"		);
	plotSigmaVsW(CBG_Poly3_Map,			0,	24,	6,	"CBG_Poly3"				);
	plotSigmaVsW(CB_Poly3_SdB_Map,		0,	24,	28,	"QEDPtShape"			);
	plotSigmaVsW(CB_Poly3_NarrMass_Map,	0,	24,	29,	"NarrMass"				);
	plotSigmaVsW(CB_Poly3_WideMass_Map,	0,	24,	30,	"WideMass"				);

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
	plotRvsX(CB_Poly3_SdB_Map,		0,	24,	28,	"QEDPtShape"			);
	plotRvsX(CB_Poly3_NarrMass_Map,	0,	24,	29,	"NarrMass"				);
	plotRvsX(CB_Poly3_WideMass_Map,	0,	24,	30,	"WideMass"				);

	plotRvsX(CB_Poly3_looseHF_Map,	4e-4,	26,	1,	"HFveto"			);
	// plotRvsX(CB_Poly3_fluxP_Map,	0,	55,	1,	"CB_Poly3_fluxP"		);
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

	plot_Uncer(SigmasTotalFitUncer_Map,	JpsiXsec_Default,	"Ws",	"Sigmas",	53,	1,	legendW,	"Signal Ext");
	plot_Uncer(SigmasFluxUncer_Map,		JpsiXsec_Default,	"Ws",	"Sigmas",	54,	1,	legendW,	"PhotonFlux");
	plot_Uncer(CB_Poly3_looseHF_Map,	JpsiXsec_Default,	"Ws",	"Sigmas",	55,	1,	legendW,	"HFveto");
	plot_Uncer(CB_Poly3_PU_Map,			JpsiXsec_Default,	"Ws",	"Sigmas",	57,	1,	legendW,	"n-PileUp");
	plot_Uncer(SigmasTnPUncer_Map,		JpsiXsec_Default,	"Ws",	"Sigmas",	59,	1,	legendW,	"TnP");

	plot_ConstUncer(JpsiXsec_Default,	"Ws",	"Sigmas",	legendW);
	legendW->Draw("same");
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

	plot_Uncer(RTotalFitUncer_Map,	JpsiXsec_Default,	"Xs",	"R",	53,	1,	legendR,	"Signal Ext");
	plot_Uncer(RFluxUncer_Map,		JpsiXsec_Default,	"Xs",	"R",	54,	1,	legendR,	"PhotonFlux");
	plot_Uncer(CB_Poly3_looseHF_Map,JpsiXsec_Default,	"Xs",	"R",	55,	1,	legendR,	"HFveto");
	plot_Uncer(CB_Poly3_PU_Map,		JpsiXsec_Default,	"Xs",	"R",	57,	1,	legendR,	"n-PileUp");
	plot_Uncer(RTnPUncer_Map,		JpsiXsec_Default,	"Xs",	"R",	59,	1,	legendR,	"TnP");

	plot_Uncer(CB_Poly3_Map,		JpsiXsec_Default,	"Xs",	"R",	56,	4,	legendR,	"IA", 	true);
	plot_ConstUncer(JpsiXsec_Default,	"Xs",	"R",	legendR);
	legendR->Draw("same");
	//-----------------------------------------------------------------------------------

	for (int i = 0; i < TotalSysUncer_Map.at("Ws").size(); ++i)
	{
		cout<<Form("W:%f,	dSigma:%f,	X:%f,	dR:%f", TotalSysUncer_Map.at("Ws")[i],	TotalSysUncer_Map.at("Sigmas_TotalSysUncer")[i],
													TotalSysUncer_Map.at("Xs")[i],	TotalSysUncer_Map.at("R_TotalSysUncer")[i])<<endl;
	}

	saveMap(TotalSysUncer_Map,	"rootfiles/TotalSysUncer_Map.root");
}