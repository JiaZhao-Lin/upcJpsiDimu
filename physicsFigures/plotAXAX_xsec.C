#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"

TGraphErrors*       ge_XsecVsY_Alice2019;
TGraphAsymmErrors* gae_XsecVsY_Alice2019;
TGraphErrors*       ge_XsecVsY_Alice2021;
TGraphAsymmErrors* gae_XsecVsY_Alice2021;
TGraphErrors*       ge_XsecVsY_LHCb2021;
TGraphAsymmErrors* gae_XsecVsY_LHCb2021;
TGraphErrors*       ge_XsecVsY_CMS2022;

void getCMSData();
void getAliceData();
void drawCMSvsAlice();

void plotAXAX_xsec( )
{
	getAliceData();
	getCMSData();
	drawCMSvsAlice();
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
}

void drawCMSvsAlice()
{
	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
	setPad(0.12, 0.08, 0.07, 0.13);

	//TH2D* htem2d = new TH2D("htem2d", "", 10, -4.1, 1.0, 10, 0, 14.0);
	TH2D* htem2d = new TH2D("htem2d", "", 10, -4.1, 0, 10, 0, 8.0);
	htem2d ->SetYTitle("d#sigma_{J/#psi}/dy (mb)");
	htem2d ->SetXTitle("y");
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

	gae_XsecVsY_Alice2019  ->SetMarkerStyle(20);
	gae_XsecVsY_Alice2019 ->SetFillColor(8);
	gae_XsecVsY_Alice2019 ->SetFillStyle(3001);
	gae_XsecVsY_Alice2019 ->Draw("2same");
	ge_XsecVsY_Alice2019  ->SetMarkerStyle(20);
	ge_XsecVsY_Alice2019  ->SetMarkerSize(0.5);
	ge_XsecVsY_Alice2019  ->Draw("psame");
	
	gae_XsecVsY_Alice2021  ->SetMarkerStyle(20);
	gae_XsecVsY_Alice2021 ->SetFillColor(kGreen);
	gae_XsecVsY_Alice2021 ->SetFillStyle(3001);
	gae_XsecVsY_Alice2021 ->Draw("2same");
	ge_XsecVsY_Alice2021  ->SetMarkerStyle(20);
	ge_XsecVsY_Alice2021  ->SetMarkerSize(0.5);
	ge_XsecVsY_Alice2021  ->Draw("psame");

	gae_XsecVsY_LHCb2021  ->SetMarkerStyle(24);
	gae_XsecVsY_LHCb2021 ->SetFillColor(kGray);
	gae_XsecVsY_LHCb2021 ->SetFillStyle(3001);
	gae_XsecVsY_LHCb2021 ->Draw("2same");
	ge_XsecVsY_LHCb2021  ->SetMarkerStyle(24);
	ge_XsecVsY_LHCb2021  ->SetMarkerSize(0.5);
	ge_XsecVsY_LHCb2021  ->Draw("psame");


    ge_XsecVsY_CMS2022  ->SetMarkerStyle(21);
	ge_XsecVsY_CMS2022  ->SetMarkerSize(0.5);
    ge_XsecVsY_CMS2022  ->SetMarkerColor(4);
    ge_XsecVsY_CMS2022  ->SetLineColor(4);
    ge_XsecVsY_CMS2022  ->SetLineWidth(2);
	ge_XsecVsY_CMS2022  ->Draw("psame");

	drawLatex(0.15, 0.85, "Pb+Pb #rightarrow Pb+Pb+J/#psi #sqrt{s_{NN}} = 5.02 TeV",      42,       0.06,      1);

	TLegend  *leg =  new TLegend(0.14, 0.60, 0.55, 0.82);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.050);
	leg->AddEntry(ge_XsecVsY_CMS2022,               "CMS",               "lp" );
	leg->AddEntry(gae_XsecVsY_Alice2019,            "ALICE 2019",        "lpf");
	leg->AddEntry(gae_XsecVsY_Alice2021,            "ALICE 2021",        "lpf");
	leg->AddEntry(gae_XsecVsY_LHCb2021,             "LHCb  2021",        "lpf");
	leg->Draw("same");

	c1->SaveAs("./outplots/xsec_jpsi_upc_AnAn_CMSvsALICE.png");
	c1->SaveAs("./outplots/xsec_jpsi_upc_AnAn_CMSvsALICE.pdf");

}

void getCMSData()
{
	const int     nPots_CMS2022 = 3;
	const double  Rap_CMS2022[nPots_CMS2022]        = { -2.25, -2.0,   -1.75 };
	const double  Xsec_CMS2022[nPots_CMS2022]       = { 3.690,  4.289,  4.582 };
	const double  RapErrLow_CMS2022[nPots_CMS2022]  = { 0.15,  0.1,	   0.15 };
	const double  RapErrHig_CMS2022[nPots_CMS2022]  = { 0.15,  0.1,	   0.15 };
	const double  XsecErr_CMS2022[nPots_CMS2022]    = { 0.041,  0.037,  0.068};

	ge_XsecVsY_CMS2022  = new TGraphErrors(
				nPots_CMS2022, 
				Rap_CMS2022,       Xsec_CMS2022, 
				RapErrLow_CMS2022, XsecErr_CMS2022);

}
