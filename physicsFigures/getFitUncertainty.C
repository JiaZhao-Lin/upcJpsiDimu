#include "plotShadowingRatio.C"


void Relative_Uncer( std::map<TString, std::vector<double>> &m1,	std::map<TString, std::vector<double>> m2, const TString param)
{
	std::vector<double> out;
	for (int i = 0; i < m1[param].size(); ++i)
	{
		out.push_back( 100 * abs(m1[param][i] - m2[param][i]) / m2[param][i] );
	}
	m1[param+"_Uncer"] = out;
}

void plot_Uncer(std::map<TString, std::vector<double>> &Map, const std::map<TString, std::vector<double>> &Map_Default,
				TString paramX, TString paramY, int style, int color, TLegend *legend, TString legendName)
{
	Relative_Uncer(Map,	Map_Default,	paramY);

	auto gr	= new TGraph(Map[paramX].size(),		&Map[paramX][0],		&Map[paramY+"_Uncer"][0]);
	gr->SetMarkerSize(1.1);
	gr->SetMarkerColor(color);
	gr->SetMarkerStyle(style);
	gr->Draw("psame");

	legend->AddEntry(gr, legendName, "p");
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
		leg->SetTextSize(0.040);

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
	std::map<TString, std::vector<double>> CB_Poly3_looseHF_Map = getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.looseHF.root");
	// std::map<TString, std::vector<double>> CB_Poly3_tightHF_Map = getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.tightHF.root");
	std::map<TString, std::vector<double>> CB_Poly3_fluxP_Map 	= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP.root", 1.05);
	std::map<TString, std::vector<double>> CB_Poly3_fluxM_Map 	= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP.root", 0.95);
	std::map<TString, std::vector<double>> CB_Poly3_TnP_Low_Map 	= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP_Low.root");
	std::map<TString, std::vector<double>> CB_Poly3_TnP_Hig_Map 	= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP_Hig.root");

	const std::map<TString, std::vector<double>> JpsiXsec_Default		= CB_Poly3_Map;


	//-----------------------------------Comparing SigmaVsW------------------------------
	plotSigmaVsW(CB_Poly3_Map,			0,	24,	1,	"CB_Poly3"				);
	plotSigmaVsW(CB_Poly4_Map,			0,	24,	2,	"CB_Poly4"				);
	plotSigmaVsW(CB_FixCBAN_Poly3_Map,	0,	24,	3,	"CB_FixCBAN_Poly3"		);
	plotSigmaVsW(CBG_Poly3_Map,			0,	24,	28,	"CBG_Poly3"				);
	plotSigmaVsW(CB_Poly3_SdB_Map,		0,	24,	6,	"CB_Poly3_SdB"			);

	plotSigmaVsW(CB_Poly3_looseHF_Map,	5,	25,	1,	"CB_Poly3_looseHF"		);
	// plotSigmaVsW(CB_Poly3_fluxP_Map,	10,	26,	1,	"CB_Poly3_fluxP"		);
	plotSigmaVsW(CB_Poly3_fluxM_Map,	10,	26,	1,	"CB_Poly3_fluxM", 	true);
	//-----------------------------------------------------------------------------------


	//-----------------------------------Comparing RVsW----------------------------------
	plotRvsX(CB_Poly3_Map,			0,	24,	1,	"CB_Poly3"				);
	plotRvsX(CB_Poly4_Map,			0,	24,	2,	"CB_Poly4"				);
	plotRvsX(CB_FixCBAN_Poly3_Map,	0,	24,	3,	"CB_FixCBAN_Poly3"		);
	plotRvsX(CBG_Poly3_Map,			0,	24,	28,	"CBG_Poly3"				);
	plotRvsX(CB_Poly3_SdB_Map,		0,	24,	6,	"CB_Poly3_SdB"			);

	plotRvsX(CB_Poly3_looseHF_Map,	4e-4,	25,	1,	"CB_Poly3_looseHF"		);
	// plotRvsX(CB_Poly3_fluxP_Map,	0,	55,	1,	"CB_Poly3_fluxP"		);
	plotRvsX(CB_Poly3_fluxM_Map,	8e-4,	26,	1,	"CB_Poly3_fluxM", 	true);
	//-----------------------------------------------------------------------------------


	//-----------------------------------W Unicertainty----------------------------------
	auto c2 = new TCanvas();
	auto Sigma_FitUncer	= new TH2D("Sigma_FitUncertainty", ";W_{#gammaPb} (GeV);#sigma Uncer.(%);", 10,14,520, 10,0,20);
	Sigma_FitUncer->Draw();

	TLegend  *legendW =  new TLegend(0.15, 0.70, 0.30, 0.90);
	legendW->SetFillStyle(0);
	legendW->SetTextSize(0.040);
	plot_Uncer(CB_Poly4_Map,			JpsiXsec_Default,	"Ws",	"Sigmas",	53,	1,	legendW,	"CB_Poly4");
	plot_Uncer(CB_FixCBAN_Poly3_Map,	JpsiXsec_Default,	"Ws",	"Sigmas",	53,	2,	legendW,	"CB_FixCBAN_Poly3");
	plot_Uncer(CBG_Poly3_Map,			JpsiXsec_Default,	"Ws",	"Sigmas",	53,	3,	legendW,	"CBG_Poly3");
	plot_Uncer(CB_Poly3_SdB_Map,		JpsiXsec_Default,	"Ws",	"Sigmas",	53,	4,	legendW,	"CB_Poly3_SdB");

	plot_Uncer(CB_Poly3_looseHF_Map,	JpsiXsec_Default,	"Ws",	"Sigmas",	54,	1,	legendW,	"CB_Poly3_looseHF");
	// plot_Uncer(CB_Poly3_tightHF_Map,	JpsiXsec_Default,	"Ws",	"Sigmas",	54,	2,	legendW,	"CB_Poly3_tightHF");
	// plot_Uncer(CB_Poly3_fluxP_Map,		JpsiXsec_Default,	"Ws",	"Sigmas",	55,	1,	legendW,	"CB_Poly3_fluxP");
	plot_Uncer(CB_Poly3_fluxM_Map,		JpsiXsec_Default,	"Ws",	"Sigmas",	55,	2,	legendW,	"CB_Poly3_fluxM");

	// plot_Uncer(CB_Poly3_TnP_Low_Map,	JpsiXsec_Default,	"Ws",	"Sigmas",	57,	1,	legendW,	"CB_Poly3_TnP_Low_Map");
	// plot_Uncer(CB_Poly3_TnP_Hig_Map,	JpsiXsec_Default,	"Ws",	"Sigmas",	57,	2,	legendW,	"CB_Poly3_TnP_Hig_Map");
	legendW->Draw("same");
	//-----------------------------------------------------------------------------------


	//-----------------------------------R Unicertainty----------------------------------
	auto c3 = new TCanvas();
	c3->SetLogx();
	auto R_FitUncer	= new TH2D("R_FitUncertainty", ";x;R^{Pb}_{g} Uncer.(%);", 10,4.0e-5,1.1e-2, 10, 0, 10);
	R_FitUncer->Draw();

	TLegend  *legendR =  new TLegend(0.15, 0.70, 0.30, 0.90);
	legendR->SetFillStyle(0);
	legendR->SetTextSize(0.040);
	plot_Uncer(CB_Poly4_Map,			JpsiXsec_Default,	"Xs",	"R",	53,	1,	legendR,	"CB_Poly4");
	plot_Uncer(CB_FixCBAN_Poly3_Map,	JpsiXsec_Default,	"Xs",	"R",	53,	2,	legendR,	"CB_FixCBAN_Poly3");
	plot_Uncer(CBG_Poly3_Map,			JpsiXsec_Default,	"Xs",	"R",	53,	3,	legendR,	"CBG_Poly3");
	plot_Uncer(CB_Poly3_SdB_Map,		JpsiXsec_Default,	"Xs",	"R",	53,	4,	legendR,	"CB_Poly3_SdB");

	plot_Uncer(CB_Poly3_looseHF_Map,	JpsiXsec_Default,	"Xs",	"R",	54,	1,	legendR,	"CB_Poly3_looseHF");
	// plot_Uncer(CB_Poly3_tightHF_Map,	JpsiXsec_Default,	"Xs",	"R",	54,	2,	legendR,	"CB_Poly3_tightHF");
	// plot_Uncer(CB_Poly3_fluxP_Map,		JpsiXsec_Default,	"Xs",	"R",	55,	1,	legendR,	"CB_Poly3_fluxP");
	plot_Uncer(CB_Poly3_fluxM_Map,		JpsiXsec_Default,	"Xs",	"R",	55,	2,	legendR,	"CB_Poly3_fluxM");

	// plot_Uncer(CB_Poly3_TnP_Low_Map,	JpsiXsec_Default,	"Xs",	"R",	57,	1,	legendR,	"CB_Poly3_TnP_Low_Map");
	// plot_Uncer(CB_Poly3_TnP_Hig_Map,	JpsiXsec_Default,	"Xs",	"R",	57,	2,	legendR,	"CB_Poly3_TnP_Hig_Map");
	legendR->Draw("same");
	//-----------------------------------------------------------------------------------
}