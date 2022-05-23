#include "./ImpulseApproximation/runUPC_AAModel.C"
#include "fit2D.C"
#include "../common/LoadSignal.C"

const double JpsiMass   = 3.096916;
const double Sqrt_s     = 5020;
const double Gamma_beam = 2672.9;
const double Mass_N     = (0.93827+0.93957) / 2;

void Cal_R_Error(std::map<TString, std::vector<double>> &ShadowRatio_ParamsMap)
{
	for (int i = 0; i < ShadowRatio_ParamsMap["Raps"].size(); ++i)
	{
		double tem = ShadowRatio_ParamsMap["R"][i]*ShadowRatio_ParamsMap["R"][i] * sqrt(pow(ShadowRatio_ParamsMap["Sigmas_Err"][i]/ShadowRatio_ParamsMap["Sigmas"][i],2) + pow(ShadowRatio_ParamsMap["Sigmas_IA_Err"][i]/ShadowRatio_ParamsMap["Sigmas_IA"][i],2) );
		ShadowRatio_ParamsMap["R_Err"].push_back( tem/ShadowRatio_ParamsMap["R"][i]/2);
	}
}

double y2x(const double y)
{
	return (JpsiMass / Sqrt_s) * exp(-y);
}

double y2W(const double y)
{
	return sqrt( (2 * Gamma_beam * Mass_N * JpsiMass) * exp(y) );
}

void plotSigmaVsW( std::map<TString, std::vector<double>> &ShadowRatio_ParamsMap, const int flag4Axis = 0 ) //0: logY only, 1: logX and LogY
{
	for (int i = 0; i < ALICE_Run2_FwdRap_y.size(); ++i)
	{
		ALICE_Run2_FwdRap_W        .push_back( y2W(ALICE_Run2_FwdRap_y[i]) );
		ALICE_Run2_FwdRap_Sigma    .push_back( ALICE_Run2_FwdRap_dSigmady[i] / ALICE_Run2_FwdRap_Flux[i] );
		ALICE_Run2_FwdRap_Sigma_Err.push_back( ALICE_Run2_FwdRap_Sigma[i] * ALICE_Run2_FwdRap_dSigmady_Err[i] / ALICE_Run2_FwdRap_dSigmady[i] );
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

	TGraphErrors* ge_CMS        		= new TGraphErrors(ShadowRatio_ParamsMap["Ws"].size(),	&ShadowRatio_ParamsMap["Ws"][0],	&ShadowRatio_ParamsMap["Sigmas"][0],	0,	&ShadowRatio_ParamsMap["Sigmas_Err"][0]	);
	
	TGraphAsymmErrors *gae_Xsec_ALICE_Run1 = new TGraphAsymmErrors(ALICE_W_Run1.size(), &ALICE_W_Run1[0],&ALICE_Xsec_Run1[0],0,0, &ALICE_XsecErr2_Run1[0], &ALICE_XsecErr1_Run1[0]);
	
	TGraphErrors* ge_ALICE_Run2_MidRap	= new TGraphErrors(ALICE_Run2_MidRap_W.size(),	&ALICE_Run2_MidRap_W[0],	&ALICE_Run2_MidRap_Sigma[0],	0,	&ALICE_Run2_MidRap_Sigma_Err[0]	);
	TGraphErrors* ge_ALICE_Run2_FwdRap	= new TGraphErrors(ALICE_Run2_FwdRap_W.size(),	&ALICE_Run2_FwdRap_W[0],	&ALICE_Run2_FwdRap_Sigma[0],	0,	&ALICE_Run2_FwdRap_Sigma_Err[0]	);

	TGraphErrors* ge_CGCnoFluct = new TGraphErrors(CGC_JpsiNoFluct_W.size(),	&CGC_JpsiNoFluct_W[0],	&CGC_JpsiNoFluct_CohXsec[0],	0,	0);
	TGraphErrors* ge_IA        	= new TGraphErrors(ShadowRatio_ParamsMap["Ws_FitIA"].size(),	&ShadowRatio_ParamsMap["Ws_FitIA"][0],	&ShadowRatio_ParamsMap["Sigmas_FitIA"][0],	0,	0);

	ge_CMS->SetMarkerStyle(20);
	ge_CMS->SetMarkerSize(1.1);
	ge_CMS->SetMarkerColor(1);
	ge_CMS->SetLineColor(1);
	ge_CMS->SetLineWidth(2);
	ge_CMS->Draw("pezsame");

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

		c->SaveAs("outplots/SigmaVsW_logY.png");
		c->SaveAs("outplots/SigmaVsW_logY.pdf");
	}

	delete c;
}

void plotRvsX( std::map<TString, std::vector<double>> &ShadowRatio_ParamsMap )
{	
	auto c = new TCanvas();

	gPad->SetTopMargin(0.05);
	gPad->SetBottomMargin(0.12);
	//gPad->SetLeftMargin(0.05);
	gPad->SetRightMargin(0.05);


	c->SetLogx();

	//TH2D* htem2d = new TH2D("htem2d", "", 10,3.0e-5,5.0e-2, 10, 0, 1);
	TH2D* htem2d = new TH2D("htem2d", "", 10,4.0e-5,1.1e-2, 10, 0.2, 1.1);
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

	TGraphErrors* ge_CMS        	= new TGraphErrors(ShadowRatio_ParamsMap["Xs"].size(),    &ShadowRatio_ParamsMap["Xs"][0], 
					&ShadowRatio_ParamsMap["R"][0],	&ShadowRatio_ParamsMap["Xs_Err"][0],	&ShadowRatio_ParamsMap["R_Err"][0]);
	TGraphErrors* ge_ALICE_Run1 	= new TGraphErrors(ALICE_x.size(),      &ALICE_x[0],      &ALICE_R[0],      0,          &ALICE_R_Err[0] );
	TGraphErrors* ge_ALICE_Run2 	= new TGraphErrors(ALICE_Run2_x.size(), &ALICE_Run2_x[0], &ALICE_Run2_R[0], 0,          &ALICE_Run2_R_Err[0] );
	TGraphErrors* ge_ALICE_Run2_Psi = new TGraphErrors(ALICE_Run2_Psi_x.size(), &ALICE_Run2_Psi_x[0], &ALICE_Run2_Psi_R[0], 0,          &ALICE_Run2_R_Err[0] );

	ge_CMS->SetMarkerStyle(20);
	ge_CMS->SetMarkerColor(1);
	ge_CMS->SetLineColor(1);
	ge_CMS->SetLineWidth(2);
	ge_CMS->Draw("pezsame");

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
	
	
	// TGraph *grshade = new TGraph(2*Theory_x_max.size());
	// for (int i = 0; i < Theory_x_max.size(); i++) {
	// 	grshade->SetPoint(i, Theory_x_max[i], Theory_R_max[i]);
	// 	grshade->SetPoint(Theory_x_max.size()+i,Theory_x_min[i],Theory_R_min[i]);
	// }
	// grshade->SetFillStyle(3013);
	// grshade->SetFillColor(16);
	// grshade->Draw("f");
	// grmin->Draw("l");
	// grmax->Draw("l");

	TLegend  *leg =  new TLegend(0.17, 0.60, 0.55, 0.80);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.050);
	leg->AddEntry(ge_CMS,         "CMS",        "lp");
	//leg->AddEntry(ge_ALICE_Run1,  "ALICE Run1", "lp"); //Run1
	leg->AddEntry(ge_ALICE_Run2,  "ALICE", "lp"); //Run2
	// leg->AddEntry(ge_ALICE_Run2_Psi,  "ALICE Run2 Psi", "lp");
	leg->Draw("same");

	drawLatex(0.15, 0.86, "Pb+Pb UPC #sqrt{s_{NN}} = 5.02 TeV",  42,        0.05,      1 );
	drawLatex(0.42, 0.22, "Nuclear suppression factor", 42, 0.06, 1);

	c->SaveAs("outplots/ShadowingRatiovsX.png");
	c->SaveAs("outplots/ShadowingRatiovsX.pdf");
	delete c;

	return ShadowRatio_ParamsMap;
}

std::map<TString, std::vector<double>> getParamsMap(const TString infile = "../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP.root", const double flux_uncer = 1.0)
{
	std::map<TString, std::vector<double>> ShadowRatio_ParamsMap = {
		{"Raps",		{}},	{"Raps_Err",		{}},
		{"Sigmas",		{}},	{"Sigmas_Err",		{}},
		{"Sigmas_IA",	{}},	{"Sigmas_IA_Err",	{}},
		{"R",			{}},	{"R_Err",			{}},
		{"Xs",			{}},	{"Xs_Err",			{}},
		{"Ws",			{}},	{"Ws_Err",			{}},

		{"Ws_FitIA",		{}},
		{"Sigmas_FitIA",	{}},	{"Sigmas_FitIA_Err",	{}}
	};

	LoadJpsiXsec loadJpsiXsec(infile);

	{auto Temp = fit2D(loadJpsiXsec.GetMap(), flux_uncer); ShadowRatio_ParamsMap["Raps"] = Temp[0]; ShadowRatio_ParamsMap["Raps_Err"] = Temp[1]; ShadowRatio_ParamsMap["Sigmas"] = Temp[2]; ShadowRatio_ParamsMap["Sigmas_Err"] = Temp[3];}
	for (int i = 0; i < ShadowRatio_ParamsMap["Raps"].size(); ++i)
	{
		ShadowRatio_ParamsMap["Xs"] 		.push_back( y2x(ShadowRatio_ParamsMap["Raps"][i]) 	);
		ShadowRatio_ParamsMap["Xs_Err"]		.push_back( 0 );
		ShadowRatio_ParamsMap["Ws"] 		.push_back( y2W(ShadowRatio_ParamsMap["Raps"][i]) 	);
		ShadowRatio_ParamsMap["Ws_Err"] 	.push_back( 0 );

		cout<<Form("y: %f,	x: %f,	W: %f", ShadowRatio_ParamsMap["Raps"][i], ShadowRatio_ParamsMap["Xs"][i], ShadowRatio_ParamsMap["Ws"][i])<<endl;
		cout<<"x from W: "<<pow(JpsiMass,2)/pow(ShadowRatio_ParamsMap["Ws"][i],2)<<endl;
	}
	
	{runUPC_AAModel(ShadowRatio_ParamsMap);}

	for (int i = 0; i < ShadowRatio_ParamsMap["Raps"].size(); ++i)
	{
		ShadowRatio_ParamsMap["R"].push_back( sqrt( ShadowRatio_ParamsMap["Sigmas"][i]/ShadowRatio_ParamsMap["Sigmas_IA"][i] ) );
	}
	
	Cal_R_Error(ShadowRatio_ParamsMap);

	return ShadowRatio_ParamsMap;
}

void plotShadowingRatio()
{
	auto ShadowRatio_ParamsMap = getParamsMap();
	plotSigmaVsW(ShadowRatio_ParamsMap);
	plotRvsX(ShadowRatio_ParamsMap);
}