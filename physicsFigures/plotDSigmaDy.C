#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"
#include "readTheory.C"
#include "Map_IO.C"

TGraphErrors*       ge_XsecVsY_Alice2019;
TGraphAsymmErrors* gae_XsecVsY_Alice2019;
TGraphErrors*       ge_XsecVsY_Alice2021;
TGraphAsymmErrors* gae_XsecVsY_Alice2021;
TGraphErrors*       ge_XsecVsY_LHCb2021;
TGraphAsymmErrors* gae_XsecVsY_LHCb2021;
TGraphErrors*       ge_XsecVsY_LHCb2022;
TGraphAsymmErrors* gae_XsecVsY_LHCb2022;
TGraphErrors*       ge_XsecVsY_CMS2022;
TGraphAsymmErrors* gae_XsecVsY_CMS2022;
TGraphErrors*       ge_XsecVsY_0n0n_CMS2022;
TGraphAsymmErrors* gae_XsecVsY_0n0n_CMS2022;
TGraphErrors*       ge_XsecVsY_0nXnSum_CMS2022;
TGraphAsymmErrors* gae_XsecVsY_0nXnSum_CMS2022;
TGraphErrors*       ge_XsecVsY_XnXn_CMS2022;
TGraphAsymmErrors* gae_XsecVsY_XnXn_CMS2022;

void getCMSData();
void getAliceData();
void drawCMSvsAlice();
void drawCMS_NeuConfig();
void drawDDP(TString Name,	TLegend* legend);

void plotDSigmaDy( )
{
	getAliceData();
	getCMSData();
	drawCMSvsAlice();
	drawCMS_NeuConfig();
}

void getAliceData()
{
	TFile* infile_Alice2019 = new TFile("./inputfiles/UPC_JpsiXsec_HEPData_Alice_2019PLB.root", "read");

	TDirectory* dir_Alice2019 = gFile->GetDirectory("Table 1");
	//
	//	TH1F* h_XsecVsY, h_XsecStatErrVsY, h_XsecSystHigErrVsY, h_XsecSystLowErrVsY;
	//
	//	dir->GetObject("Hist1D_y1",         h_XsecVsY           );
	//	dir->GetObject("Hist1D_y1_e1",      h_XsecStatErrVsY    );
	//	dir->GetObject("Hist1D_y1_e2plus",  h_XsecSystHigErrVsY );
	//	dir->GetObject("Hist1D_y1_e2minus", h_XsecSystLowErrVsY );
	//
	//TGraphAsymmErrors* gae;
	//dir_Alice2019->GetObject("Graph1D_y1", gae);

	//for(int ib=0; ib<7; ib++)
	//{
	//	cout<<"ib: "<<ib<<" "<<gae->GetErrorYhigh(ib)<<endl;
	//}

	const int     nPots_Alice2019  = 6;
	const double  Rap_Alice2019[nPots_Alice2019]        = {-3.875, -3.625, -3.375, -3.125, -2.875, -2.625};
	const double  Xsec_Alice2019[nPots_Alice2019]       = {1.615,   1.938,  2.377,  2.831,  3.018,  3.531};
	const double  RapErrLow_Alice2019[nPots_Alice2019]  = {0.125,   0.125,  0.125,  0.125,  0.125,  0.125};
	const double  RapErrHig_Alice2019[nPots_Alice2019]  = {0.125,   0.125,  0.125,  0.125,  0.125,  0.125};
	const double  XsecErrLow_Alice2019[nPots_Alice2019] = {0.147,   0.190,  0.229,  0.280,  0.294,  0.362};
	const double  XsecErrHig_Alice2019[nPots_Alice2019] = {0.135,   0.166,  0.212,  0.253,  0.259,  0.294};
	const double  XsecErr_Alice2019[nPots_Alice2019]    = {0.060,   0.042,  0.040,  0.047,  0.061,  0.139};

	ge_XsecVsY_Alice2019  = new TGraphErrors(
				nPots_Alice2019, 
				Rap_Alice2019,       Xsec_Alice2019, 
				RapErrLow_Alice2019, XsecErr_Alice2019);
	gae_XsecVsY_Alice2019 = new TGraphAsymmErrors(
				nPots_Alice2019,
				Rap_Alice2019,        Xsec_Alice2019,
				RapErrLow_Alice2019,  RapErrHig_Alice2019,
				XsecErrLow_Alice2019, XsecErrHig_Alice2019);


	const int     nPots_Alice2021  = 6;
	const double  Rap_Alice2021[nPots_Alice2021]        = {-0.575, -0.250, -0.075, 0.075, 0.250, 0.575 };
	const double  Xsec_Alice2021[nPots_Alice2021]       = { 3.870,  4.210,  4.070, 4.070, 4.210, 3.870 };
	const double  RapErrLow_Alice2021[nPots_Alice2021]  = { 0.225,  0.100,  0.075, 0.075, 0.100, 0.225 };
	const double  RapErrHig_Alice2021[nPots_Alice2021]  = { 0.225,  0.100,  0.075, 0.075, 0.100, 0.225 };
	const double  XsecErrLow_Alice2021[nPots_Alice2021] = { 0.220,  0.240,  0.230, 0.230, 0.240, 0.220 };
	const double  XsecErrHig_Alice2021[nPots_Alice2021] = { 0.220,  0.240,  0.230, 0.230, 0.240, 0.220 };
	const double  XsecErr_Alice2021[nPots_Alice2021]    = { 0.110,  0.110,  0.110, 0.110, 0.110, 0.110 };

	ge_XsecVsY_Alice2021  = new TGraphErrors(
				nPots_Alice2021, 
				Rap_Alice2021,       Xsec_Alice2021, 
				RapErrLow_Alice2021, XsecErr_Alice2021);
	gae_XsecVsY_Alice2021 = new TGraphAsymmErrors(
				nPots_Alice2021,
				Rap_Alice2021,        Xsec_Alice2021,
				RapErrLow_Alice2021,  RapErrHig_Alice2021,
				XsecErrLow_Alice2021, XsecErrHig_Alice2021);

	const int     nPots_LHCb2021  = 5;
	const double  dy_LHCb = 0.5;
	const double  Rap_LHCb2021[nPots_LHCb2021]        = {-4.25, -3.75, -3.25, -2.75, -2.25 };
	const double  Xsec_LHCb2021[nPots_LHCb2021]       = { 0.48/dy_LHCb,  0.65/dy_LHCb,  0.89/dy_LHCb,  1.09/dy_LHCb,  1.35/dy_LHCb };
	const double  RapErrLow_LHCb2021[nPots_LHCb2021]  = { 0.25,  0.25,  0.25,  0.25,  0.25 };
	const double  RapErrHig_LHCb2021[nPots_LHCb2021]  = { 0.25,  0.25,  0.25,  0.25,  0.25 };
	const double  XsecErrLow_LHCb2021[nPots_LHCb2021] = { 0.02/dy_LHCb,  0.03/dy_LHCb,  0.04/dy_LHCb,  0.05/dy_LHCb,  0.06/dy_LHCb };
	const double  XsecErrHig_LHCb2021[nPots_LHCb2021] = { 0.02/dy_LHCb,  0.03/dy_LHCb,  0.04/dy_LHCb,  0.05/dy_LHCb,  0.06/dy_LHCb };
	const double  XsecErr_LHCb2021[nPots_LHCb2021]    = { 0.09/dy_LHCb,  0.06/dy_LHCb,  0.07/dy_LHCb,  0.09/dy_LHCb,  0.19/dy_LHCb };

	ge_XsecVsY_LHCb2021  = new TGraphErrors(
			nPots_LHCb2021,
			Rap_LHCb2021,       Xsec_LHCb2021,
			RapErrLow_LHCb2021, XsecErr_LHCb2021);
	gae_XsecVsY_LHCb2021 = new TGraphAsymmErrors(
			nPots_LHCb2021,
			Rap_LHCb2021,        Xsec_LHCb2021,
			RapErrLow_LHCb2021,  RapErrHig_LHCb2021,
			XsecErrLow_LHCb2021, XsecErrHig_LHCb2021);

	const int     nPots_LHCb2022  = 5;
	const double  Rap_LHCb2022[nPots_LHCb2022]        = {-4.25, -3.75, -3.25, -2.75, -2.25 };
	const double  Xsec_LHCb2022[nPots_LHCb2022]       = { 1.181,  1.922,  2.395,  2.896,  3.392 };
	const double  RapErrLow_LHCb2022[nPots_LHCb2022]  = { 0.25,  0.25,  0.25,  0.25,  0.25 };
	const double  RapErrHig_LHCb2022[nPots_LHCb2022]  = { 0.25,  0.25,  0.25,  0.25,  0.25 };
	const double  XsecErrLow_LHCb2022[nPots_LHCb2022] = { TMath::Hypot(0.049, 0.052),	TMath::Hypot(0.072, 0.084),  TMath::Hypot(0.089, 0.105),  TMath::Hypot(0.117, 0.127),  TMath::Hypot(0.165, 0.147) };
	const double  XsecErrHig_LHCb2022[nPots_LHCb2022] = { TMath::Hypot(0.049, 0.052),	TMath::Hypot(0.072, 0.084),  TMath::Hypot(0.089, 0.105),  TMath::Hypot(0.117, 0.127),  TMath::Hypot(0.165, 0.147) };
	const double  XsecErr_LHCb2022[nPots_LHCb2022]    = { 0.054,  0.039,  0.040,  0.053,  0.108 };

	ge_XsecVsY_LHCb2022  = new TGraphErrors(
			nPots_LHCb2022,
			Rap_LHCb2022,       Xsec_LHCb2022,
			RapErrLow_LHCb2022, XsecErr_LHCb2022);
	gae_XsecVsY_LHCb2022 = new TGraphAsymmErrors(
			nPots_LHCb2022,
			Rap_LHCb2022,        Xsec_LHCb2022,
			RapErrLow_LHCb2022,  RapErrHig_LHCb2022,
			XsecErrLow_LHCb2022, XsecErrHig_LHCb2022);
}

void drawCMSvsAlice()
{
	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
	setPad(0.12, 0.08, 0.07, 0.13);

	//TH2D* htem2d = new TH2D("htem2d", "", 10, -4.1, 1.0, 10, 0, 14.0);
	TH2D* htem2d = new TH2D("AnAn", "AnAn;y;d#sigma_{J/#psi}/dy (mb);", 10, -4.5, 0, 10, 0, 8.0);
	htem2d ->GetYaxis()->CenterTitle();
	//htem2d ->GetYaxis()->SetNdivisions(6);
	htem2d ->GetYaxis()->SetTitleSize(0.07);
	htem2d ->GetYaxis()->SetTitleOffset(0.63);
	htem2d ->GetYaxis()->SetLabelSize(0.06);
	htem2d ->GetXaxis()->CenterTitle();
	htem2d ->GetXaxis()->SetTitleSize(0.07);
	htem2d ->GetXaxis()->SetTitleOffset(0.75);
	htem2d ->GetXaxis()->SetLabelSize(0.06);
	//htem2d ->SetTickLength(0.08);

	htem2d->Draw();

	gae_XsecVsY_Alice2019  ->SetMarkerStyle(21);
	gae_XsecVsY_Alice2019 ->SetFillColor(8);
	gae_XsecVsY_Alice2019 ->SetFillStyle(3001);
	gae_XsecVsY_Alice2019 ->Draw("2same");
	ge_XsecVsY_Alice2019  ->SetMarkerStyle(21);
	ge_XsecVsY_Alice2019  ->SetMarkerSize(1.);
	ge_XsecVsY_Alice2019  ->Draw("pezsame");
	
	gae_XsecVsY_Alice2021  ->SetMarkerStyle(21);
	gae_XsecVsY_Alice2021 ->SetFillColor(kGreen);
	gae_XsecVsY_Alice2021 ->SetFillStyle(3001);
	gae_XsecVsY_Alice2021 ->Draw("2same");
	ge_XsecVsY_Alice2021  ->SetMarkerStyle(21);
	ge_XsecVsY_Alice2021  ->SetMarkerSize(1.);
	ge_XsecVsY_Alice2021  ->Draw("pezsame");

	gae_XsecVsY_LHCb2022  ->SetMarkerStyle(24);
	gae_XsecVsY_LHCb2022 ->SetFillColor(kGray);
	gae_XsecVsY_LHCb2022 ->SetFillStyle(3001);
	gae_XsecVsY_LHCb2022 ->Draw("2same");
	ge_XsecVsY_LHCb2022  ->SetMarkerStyle(26);
	ge_XsecVsY_LHCb2022  ->SetMarkerSize(1.);
	ge_XsecVsY_LHCb2022  ->Draw("pezsame");

	gae_XsecVsY_CMS2022 ->SetMarkerStyle(20);
	gae_XsecVsY_CMS2022 ->SetMarkerColor(2);
	gae_XsecVsY_CMS2022 ->SetLineColor(2);
	gae_XsecVsY_CMS2022 ->SetFillColorAlpha(16, 0.5);
	gae_XsecVsY_CMS2022 ->SetFillStyle(1001);
	gae_XsecVsY_CMS2022 ->Draw("2same");
    ge_XsecVsY_CMS2022  ->SetMarkerStyle(20);
	ge_XsecVsY_CMS2022  ->SetMarkerSize(1.);
    ge_XsecVsY_CMS2022  ->SetMarkerColor(2);
    ge_XsecVsY_CMS2022  ->SetLineColor(2);
    ge_XsecVsY_CMS2022  ->SetLineWidth(2);
	ge_XsecVsY_CMS2022  ->Draw("pezsame");

	drawLatex(0.15, 0.85, "Pb+Pb #rightarrow Pb+Pb+J/#psi #sqrt{s_{NN}} = 5.02 TeV",      42,       0.06,      1);

	TLegend  *leg =  new TLegend(0.14, 0.60, 0.55, 0.82);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.050);
	leg->AddEntry(gae_XsecVsY_CMS2022,              "CMS",               "lpf" );
	leg->AddEntry(gae_XsecVsY_Alice2019,            "ALICE 2019",        "lpf");
	leg->AddEntry(gae_XsecVsY_Alice2021,            "ALICE 2021",        "lpf");
	leg->AddEntry(gae_XsecVsY_LHCb2022,             "LHCb  2022",        "lpf");
	// drawDDP("DDP",	leg);
	leg->Draw("same");

	c1->SaveAs("./outplots/xsec_jpsi_upc_AnAn_CMSvsALICE.png");
	c1->SaveAs("./outplots/xsec_jpsi_upc_AnAn_CMSvsALICE.pdf");

}

void drawCMS_NeuConfig()
{
	TCanvas* c2 = new TCanvas("c2", "c2", 0, 0, 800, 600);
	c2->Divide(2,2);
	setPad(0.12, 0.08, 0.07, 0.13);

	c2->cd(1);
    c2->GetPad(1)->SetRightMargin(0.01);
	TH2D* htem2d_0n0n = new TH2D("htem2d_0n0n", "htem2d_0n0n;y;d#sigma_{J/#psi}/dy (mb)", 10, -4.1, 0, 10, 0.0, 7.0);
	htem2d_0n0n ->GetYaxis()->CenterTitle();
	htem2d_0n0n ->GetYaxis()->SetNdivisions(6);
	htem2d_0n0n ->GetXaxis()->SetNdivisions(6);
	htem2d_0n0n ->GetYaxis()->SetTitleSize(0.065);
	htem2d_0n0n ->GetYaxis()->SetTitleOffset(0.8);
	htem2d_0n0n ->GetYaxis()->SetLabelSize(0.05);
	htem2d_0n0n ->GetXaxis()->CenterTitle();
	htem2d_0n0n ->GetXaxis()->SetTitleSize(0.065);
	htem2d_0n0n ->GetXaxis()->SetTitleOffset(0.75);
	htem2d_0n0n ->GetXaxis()->SetLabelSize(0.05);
	//htem2d_0n0n ->SetTickLength(0.08);
	htem2d_0n0n->Draw();
	drawLatex(0.15, 0.85, "Pb+Pb #rightarrow Pb+Pb+J/#psi",      42,       0.053,      1);
    drawLatex(0.8, 0.8, "#bf{0n0n}",      42,       0.075,      1);
    drawLatex(0.58,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
    drawLatex(0.15,0.94,"#bf{CMS}",42, 0.05, 1);

	gae_XsecVsY_0n0n_CMS2022 ->SetMarkerStyle(20);
    gae_XsecVsY_0n0n_CMS2022 ->SetMarkerColor(2);
    gae_XsecVsY_0n0n_CMS2022 ->SetLineColor(2);
	gae_XsecVsY_0n0n_CMS2022 ->SetFillColorAlpha(16, 0.5);
	gae_XsecVsY_0n0n_CMS2022 ->SetFillStyle(1001);
	gae_XsecVsY_0n0n_CMS2022 ->Draw("2same");
	ge_XsecVsY_0n0n_CMS2022  ->SetMarkerStyle(20);
	ge_XsecVsY_0n0n_CMS2022  ->SetMarkerSize(1.5);
    ge_XsecVsY_0n0n_CMS2022  ->SetMarkerColor(2);
    ge_XsecVsY_0n0n_CMS2022  ->SetLineColor(2);
	ge_XsecVsY_0n0n_CMS2022  ->Draw("pezsame");

	TLegend  *leg_0n0n =  new TLegend(0.14, 0.60, 0.55, 0.82);
	leg_0n0n->SetFillStyle(0);
	leg_0n0n->SetFillColor(0);
	leg_0n0n->SetTextSize(0.050);
	leg_0n0n->AddEntry(gae_XsecVsY_0n0n_CMS2022,		"CMS",		"lpf");
	drawLTA("0n0n",leg_0n0n);
	leg_0n0n->Draw("same");
	
	c2->cd(2);
    c2->GetPad(2)->SetRightMargin(0.01);
	TH2D* htem2d_0nXnSum = new TH2D("htem2d_0nXnSum", "htem2d_0nXnSum;y;d#sigma_{J/#psi}/dy (mb)", 10, -4.1, 0, 10, 0.0, 2.0);
	htem2d_0nXnSum ->GetYaxis()->CenterTitle();
	htem2d_0nXnSum ->GetYaxis()->SetNdivisions(6);
    htem2d_0nXnSum ->GetXaxis()->SetNdivisions(6);
	htem2d_0nXnSum ->GetYaxis()->SetTitleSize(0.065);
	htem2d_0nXnSum ->GetYaxis()->SetTitleOffset(0.8);
	htem2d_0nXnSum ->GetYaxis()->SetLabelSize(0.05);
	htem2d_0nXnSum ->GetXaxis()->CenterTitle();
	htem2d_0nXnSum ->GetXaxis()->SetTitleSize(0.065);
	htem2d_0nXnSum ->GetXaxis()->SetTitleOffset(0.75);
	htem2d_0nXnSum ->GetXaxis()->SetLabelSize(0.05);
	//htem2d_0nXnSum ->SetTickLength(0.08);
	htem2d_0nXnSum->Draw();
    drawLatex(0.15, 0.85, "Pb+Pb #rightarrow Pb+Pb'+J/#psi",      42,       0.053,      1);
    drawLatex(0.8, 0.8, "#bf{0nXn}",      42,       0.075,      1);
    drawLatex(0.58,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
    drawLatex(0.15,0.94,"#bf{CMS}",42, 0.05, 1);

	gae_XsecVsY_0nXnSum_CMS2022 ->SetMarkerStyle(20);
    gae_XsecVsY_0nXnSum_CMS2022 ->SetMarkerColor(2);
    gae_XsecVsY_0nXnSum_CMS2022 ->SetLineColor(2);
	gae_XsecVsY_0nXnSum_CMS2022 ->SetFillColorAlpha(16, 0.5);
	gae_XsecVsY_0nXnSum_CMS2022 ->SetFillStyle(1001);
	gae_XsecVsY_0nXnSum_CMS2022 ->Draw("2same");
	ge_XsecVsY_0nXnSum_CMS2022  ->SetMarkerStyle(20);
	ge_XsecVsY_0nXnSum_CMS2022  ->SetMarkerSize(1.5);
    ge_XsecVsY_0nXnSum_CMS2022  ->SetMarkerColor(2);
    ge_XsecVsY_0nXnSum_CMS2022  ->SetLineColor(2);
	ge_XsecVsY_0nXnSum_CMS2022  ->Draw("pezsame");

	TLegend  *leg_0nXnSum =  new TLegend(0.14, 0.60, 0.55, 0.82);
	leg_0nXnSum->SetFillStyle(0);
	leg_0nXnSum->SetFillColor(0);
	leg_0nXnSum->SetTextSize(0.050);
	leg_0nXnSum->AddEntry(gae_XsecVsY_0nXnSum_CMS2022,		"CMS",		"lpf");
	drawLTA("0nXnSum",leg_0nXnSum);
	leg_0nXnSum->Draw("same");

	c2->cd(3);
    c2->GetPad(3)->SetRightMargin(0.01);
	TH2D* htem2d_XnXn = new TH2D("htem2d_XnXn", "htem2d_XnXn;y;d#sigma_{J/#psi}/dy (mb)", 10, -4.1, 0, 10, 0.0, 0.5);
	htem2d_XnXn ->GetYaxis()->CenterTitle();
	htem2d_XnXn ->GetYaxis()->SetNdivisions(6);
    htem2d_XnXn ->GetXaxis()->SetNdivisions(6);
	htem2d_XnXn ->GetYaxis()->SetTitleSize(0.065);
	htem2d_XnXn ->GetYaxis()->SetTitleOffset(0.9);
	htem2d_XnXn ->GetYaxis()->SetLabelSize(0.05);
	htem2d_XnXn ->GetXaxis()->CenterTitle();
	htem2d_XnXn ->GetXaxis()->SetTitleSize(0.065);
	htem2d_XnXn ->GetXaxis()->SetTitleOffset(0.75);
	htem2d_XnXn ->GetXaxis()->SetLabelSize(0.05);
	//htem2d_XnXn ->SetTickLength(0.08);
	htem2d_XnXn->Draw();
    drawLatex(0.15, 0.85, "Pb+Pb #rightarrow Pb'+Pb'+J/#psi",      42,       0.053,      1);
    drawLatex(0.8,  0.8,  "#bf{XnXn}",      42,       0.075,      1);
    drawLatex(0.58, 0.94, "PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
    drawLatex(0.15, 0.94, "#bf{CMS}",42, 0.05, 1);

	gae_XsecVsY_XnXn_CMS2022 ->SetMarkerStyle(20);
    gae_XsecVsY_XnXn_CMS2022 ->SetMarkerColor(2);
    gae_XsecVsY_XnXn_CMS2022 ->SetLineColor(2);
	gae_XsecVsY_XnXn_CMS2022 ->SetFillColorAlpha(16, 0.5);
	gae_XsecVsY_XnXn_CMS2022 ->SetFillStyle(1001);
	gae_XsecVsY_XnXn_CMS2022 ->Draw("2same");
	ge_XsecVsY_XnXn_CMS2022  ->SetMarkerStyle(20);
	ge_XsecVsY_XnXn_CMS2022  ->SetMarkerSize(1.5);
    ge_XsecVsY_XnXn_CMS2022  ->SetMarkerColor(2);
    ge_XsecVsY_XnXn_CMS2022  ->SetLineColor(2);
	ge_XsecVsY_XnXn_CMS2022  ->Draw("pezsame");

	TLegend  *leg_XnXn =  new TLegend(0.14, 0.60, 0.55, 0.82);
	leg_XnXn->SetFillStyle(0);
	leg_XnXn->SetFillColor(0);
	leg_XnXn->SetTextSize(0.050);
	leg_XnXn->AddEntry(gae_XsecVsY_XnXn_CMS2022,		"CMS",		"lpf");
	drawLTA("XnXn",leg_XnXn);
	leg_XnXn->Draw("same");

	c2->cd(4);
    c2->GetPad(4)->SetRightMargin(0.01);
	TH2D* htem2d_AnAn = new TH2D("htem2d_AnAn", "htem2d_AnAn;y;d#sigma_{J/#psi}/dy (mb);", 10, -4.1, 0, 10, 0, 8.0);
	htem2d_AnAn ->GetYaxis()->CenterTitle();
	htem2d_AnAn ->GetYaxis()->SetNdivisions(6);
    htem2d_AnAn ->GetXaxis()->SetNdivisions(6);
	htem2d_AnAn ->GetYaxis()->SetTitleSize(0.065);
	htem2d_AnAn ->GetYaxis()->SetTitleOffset(0.63);
	htem2d_AnAn ->GetYaxis()->SetLabelSize(0.05);
	htem2d_AnAn ->GetXaxis()->CenterTitle();
	htem2d_AnAn ->GetXaxis()->SetTitleSize(0.065);
	htem2d_AnAn ->GetXaxis()->SetTitleOffset(0.75);
	htem2d_AnAn ->GetXaxis()->SetLabelSize(0.05);
	//htem2d_AnAn ->SetTickLength(0.08);
	htem2d_AnAn->Draw();
    drawLatex(0.15, 0.85, "Pb+Pb #rightarrow Pb'+Pb'+J/#psi",      42,       0.053,      1);
    drawLatex(0.75, 0.8,  "#bf{AnAn}",                             42,       0.075,      1);
    //drawLatex(0.75, 0.8, "#bf{Inclusive}",      42,       0.065,      1);
    drawLatex(0.58,0.94,  "PbPb 1.52 nb^{-1} (5.02 TeV)",          42,       0.050,      1);
    drawLatex(0.15,0.94,  "#bf{CMS}",                              42,       0.050,      1);

	gae_XsecVsY_Alice2019 ->SetMarkerStyle(25);
	gae_XsecVsY_Alice2019 ->SetMarkerColor(4);
	gae_XsecVsY_Alice2019 ->SetLineColor(4);
	gae_XsecVsY_Alice2019 ->SetFillColor(kGray);
	gae_XsecVsY_Alice2019 ->SetFillStyle(3001);
	ge_XsecVsY_Alice2019  ->SetMarkerStyle(25);
	ge_XsecVsY_Alice2019  ->SetMarkerSize(1.5);
    ge_XsecVsY_Alice2019  ->SetMarkerColor(4);
    ge_XsecVsY_Alice2019  ->SetLineColor(4);
	
	gae_XsecVsY_Alice2021 ->SetMarkerStyle(24);
	gae_XsecVsY_Alice2021 ->SetMarkerColor(4);
	gae_XsecVsY_Alice2021 ->SetLineColor(4);
	gae_XsecVsY_Alice2021 ->SetFillColor(kGray);
	gae_XsecVsY_Alice2021 ->SetFillStyle(3001);
	ge_XsecVsY_Alice2021  ->SetMarkerStyle(24);
	ge_XsecVsY_Alice2021  ->SetMarkerSize(1.6);
    ge_XsecVsY_Alice2021  ->SetMarkerColor(4);
    ge_XsecVsY_Alice2021  ->SetLineColor(4);

	gae_XsecVsY_LHCb2022  ->SetMarkerStyle(26);
	gae_XsecVsY_LHCb2022 ->SetFillColor(kGray);
	gae_XsecVsY_LHCb2022 ->SetFillStyle(3001);
	ge_XsecVsY_LHCb2022  ->SetMarkerStyle(26);
	ge_XsecVsY_LHCb2022  ->SetMarkerSize(1.5);

	gae_XsecVsY_CMS2022 ->SetMarkerStyle(20);
	gae_XsecVsY_CMS2022 ->SetMarkerColor(2);
	gae_XsecVsY_CMS2022 ->SetFillColorAlpha(16, 0.5);
	gae_XsecVsY_CMS2022 ->SetFillStyle(1001);
    ge_XsecVsY_CMS2022  ->SetMarkerStyle(20);
	ge_XsecVsY_CMS2022  ->SetMarkerSize(1.5);
    ge_XsecVsY_CMS2022  ->SetMarkerColor(2);
    ge_XsecVsY_CMS2022  ->SetLineColor(2);
    ge_XsecVsY_CMS2022  ->SetLineWidth(2);


	TLegend  *leg_AnAn   =  new TLegend(0.14, 0.55, 0.55, 0.82);
	leg_AnAn->SetFillStyle(0);
	leg_AnAn->SetFillColor(0);
	leg_AnAn->SetTextSize(0.050);
	leg_AnAn->AddEntry(gae_XsecVsY_CMS2022,			"CMS",               "lpf" );
	leg_AnAn->AddEntry(gae_XsecVsY_Alice2019,		"ALICE 2019",        "lpf");
	leg_AnAn->AddEntry(gae_XsecVsY_Alice2021,		"ALICE 2021",        "lpf");
	leg_AnAn->AddEntry(gae_XsecVsY_LHCb2022,		"LHCb  2022",        "lpf");
	leg_AnAn->Draw("same");
	
	TLegend  *leg_LTA =  new TLegend(0.48, 0.18, 0.70, 0.34);
	leg_LTA->SetTextSize(0.048);
	drawLTA("AnAn", leg_LTA);
	// drawbBK("Xsec",leg_LTA);
	//drawCD("Xsec",  leg_LTA);
	leg_LTA->Draw("same");
	
	TLegend  *leg_CD =  new TLegend(0.71, 0.15, 0.97, 0.35);
	leg_CD->SetTextSize(0.050);
	drawCD("Xsec",  leg_CD);
	leg_CD->Draw("same");


	gae_XsecVsY_Alice2019 ->Draw("2same");
	ge_XsecVsY_Alice2019  ->Draw("pezsame");
	gae_XsecVsY_Alice2021 ->Draw("2same");
	ge_XsecVsY_Alice2021  ->Draw("pezsame");
	gae_XsecVsY_LHCb2022 ->Draw("2same");
	ge_XsecVsY_LHCb2022  ->Draw("pezsame");
	gae_XsecVsY_CMS2022 ->Draw("2same");
	ge_XsecVsY_CMS2022  ->Draw("pezsame");

	c2->SaveAs("outplots/Xsec_NeuConfig.png");
	c2->SaveAs("outplots/Xsec_NeuConfig.pdf");
}

void drawDDP(TString Name,	TLegend* legend)
{
	auto Xsec_AnAn_DataDrivenPrediction	= readMap("rootfiles/DataDrivenPrediction.root");

	auto gr_AnAn_DataDrivenPrediction = new TGraph(Xsec_AnAn_DataDrivenPrediction.at("Xsec_AnAn_Prediction").size(),
												Xsec_AnAn_DataDrivenPrediction.at("Raps").data(),
												Xsec_AnAn_DataDrivenPrediction.at("Xsec_AnAn_Prediction").data());

	gr_AnAn_DataDrivenPrediction->Draw("same");
	gr_AnAn_DataDrivenPrediction->SetLineColor(1);
	gr_AnAn_DataDrivenPrediction->SetLineWidth(2);
	legend->AddEntry(gr_AnAn_DataDrivenPrediction,	Name.Data(),		"l" );
}

void getCMSData()
{
	const int template_option  = 1; //0:OldCohJpsi; 1:NewCohJpsi (w R+1fm);
	const TString template_Name[2] = {"", "_NewCohJpsi"};

	auto ParamsMap 	= readMap(Form("rootfiles/Results%s_Map.root", template_Name[template_option].Data()));
	auto TotalSysUncer_Map 		= readMap(Form("rootfiles/TotalSysUncer%s_Map.root", template_Name[template_option].Data()));

	std::vector<double> X_AXIS_ERR = ParamsMap.at("RapErr");
	auto Xsec_AnAn_TotalSysErr = ParamsMap.at("Xsec_AnAn");
	auto Xsec_0n0n_TotalSysErr = ParamsMap.at("Xsec_0n0n");
	auto Xsec_0nXnSum_TotalSysErr = ParamsMap.at("Xsec_0nXnSum");
	auto Xsec_XnXn_TotalSysErr = ParamsMap.at("Xsec_XnXn");
	for (int i = 0; i < Xsec_AnAn_TotalSysErr.size(); ++i)
	{
		Xsec_AnAn_TotalSysErr[i] 	*= 0.01 * TotalSysUncer_Map.at("Xsec_AnAn_TotalSysUncer")[i];
		Xsec_0n0n_TotalSysErr[i] 	*= 0.01 * TotalSysUncer_Map.at("Xsec_0n0n_TotalSysUncer")[i];
		Xsec_0nXnSum_TotalSysErr[i] *= 0.01 * TotalSysUncer_Map.at("Xsec_0nXnSum_TotalSysUncer")[i];
		Xsec_XnXn_TotalSysErr[i] 	*= 0.01 * TotalSysUncer_Map.at("Xsec_XnXn_TotalSysUncer")[i];
		cout<<Form("Rap:%f,	AnAn_TotalSysErr:%f,	0n0n_TotalSysErr:%f,	0nXnSum_TotalSysErr:%f,	XnXn_TotalSysErr:%f",
				ParamsMap.at("Rap")[i],Xsec_AnAn_TotalSysErr[i],Xsec_0n0n_TotalSysErr[i],Xsec_0nXnSum_TotalSysErr[i],Xsec_XnXn_TotalSysErr[i])<<endl;
	}

	ge_XsecVsY_CMS2022  = new TGraphErrors(
				ParamsMap.at("Rap").size(), 
				ParamsMap.at("Rap").data(),       ParamsMap.at("Xsec_AnAn").data(), 
				ParamsMap.at("RapErr").data(),    ParamsMap.at("XsecErr_AnAn").data());
	gae_XsecVsY_CMS2022 = new TGraphAsymmErrors(
				ParamsMap.at("Rap").size(),
				&ParamsMap.at("Rap")[0],        &ParamsMap.at("Xsec_AnAn")[0],
				&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
				&Xsec_AnAn_TotalSysErr[0], &Xsec_AnAn_TotalSysErr[0]);

	ge_XsecVsY_0n0n_CMS2022  = new TGraphErrors(
				ParamsMap.at("Rap").size(), 
				ParamsMap.at("Rap").data(),       ParamsMap.at("Xsec_0n0n").data(), 
				ParamsMap.at("RapErr").data(),    ParamsMap.at("XsecErr_0n0n").data());
	gae_XsecVsY_0n0n_CMS2022 = new TGraphAsymmErrors(
				ParamsMap.at("Rap").size(),
				&ParamsMap.at("Rap")[0],        &ParamsMap.at("Xsec_0n0n")[0],
				&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
				&Xsec_0n0n_TotalSysErr[0], &Xsec_0n0n_TotalSysErr[0]);

	ge_XsecVsY_0nXnSum_CMS2022  = new TGraphErrors(
				ParamsMap.at("Rap").size(), 
				ParamsMap.at("Rap").data(),       ParamsMap.at("Xsec_0nXnSum").data(), 
				ParamsMap.at("RapErr").data(),    ParamsMap.at("XsecErr_0nXnSum").data());
	gae_XsecVsY_0nXnSum_CMS2022 = new TGraphAsymmErrors(
				ParamsMap.at("Rap").size(),
				&ParamsMap.at("Rap")[0],        &ParamsMap.at("Xsec_0nXnSum")[0],
				&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
				&Xsec_0nXnSum_TotalSysErr[0], &Xsec_0nXnSum_TotalSysErr[0]);

	ge_XsecVsY_XnXn_CMS2022  = new TGraphErrors(
				ParamsMap.at("Rap").size(), 
				ParamsMap.at("Rap").data(),       ParamsMap.at("Xsec_XnXn").data(), 
				ParamsMap.at("RapErr").data(),    ParamsMap.at("XsecErr_XnXn").data());
	gae_XsecVsY_XnXn_CMS2022 = new TGraphAsymmErrors(
				ParamsMap.at("Rap").size(),
				&ParamsMap.at("Rap")[0],        &ParamsMap.at("Xsec_XnXn")[0],
				&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
				&Xsec_XnXn_TotalSysErr[0], &Xsec_XnXn_TotalSysErr[0]);
}
