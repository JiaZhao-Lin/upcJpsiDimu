#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"

void plotRvsN()
{
	const double x[3]    = {0.5, 1.5, 2.5};
	const double xerr[3] = {0, 0, 0};
	const double y_negRap[3]    = {0.146,  0.196,  0.319 };
	const double yerr_negRap[3] = {0.0205, 0.0407, 0.0579};

	TGraphErrors* ge_negRap = new TGraphErrors(3, x, y_negRap, xerr, yerr_negRap);
	
	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
	setPad(0.12, 0.08, 0.07, 0.13);

	TH2D* htem2d = new TH2D("htem2d", "", 3, 0.0, 3.0, 10, 0.09, 0.42);
	//TH2D* htem2d = new TH2D("htem2d", "", 10, -4.1, 1.0, 10, 0, 8.0);
	htem2d ->SetYTitle("#sigma_{#psi(2S)}/#sigma_{J/#psi}");
	htem2d ->SetXTitle("");
	htem2d ->GetYaxis()->CenterTitle();
	htem2d ->GetYaxis()->SetTitleSize(0.07);
	htem2d ->GetYaxis()->SetTitleOffset(0.80);
	htem2d ->GetYaxis()->SetLabelSize(0.04);
	htem2d ->GetXaxis()->SetNdivisions(3);
	htem2d ->GetXaxis()->CenterTitle();
	htem2d ->GetXaxis()->SetTitleSize(0.07);
	htem2d ->GetXaxis()->SetTitleOffset(0.75);
	htem2d ->GetXaxis()->SetLabelSize(0.09);
	htem2d ->SetTickLength(0.04);
	
	const TString labelNames[3] = {"0n0n", "0nXnSum", "XnXn"};
	for(int i=0;i<3;i++) htem2d->GetXaxis()->SetBinLabel( i+1, labelNames[i] );
	htem2d->Draw();
	
	ge_negRap->SetMarkerStyle(21);
	ge_negRap->SetMarkerSize(2);
	ge_negRap->SetLineWidth(2);
	ge_negRap->Draw("pesame");

	drawLatex(0.15, 0.85, "UPC Pb+Pb #sqrt{s_{NN}} = 5.02 TeV",      42,       0.06,      1);
	drawLatex(0.18, 0.78, "-2.0 < y < -1.6",      42,       0.05,      1);
	
	c1->SaveAs("./outplots/RvsN_UPC.png");
	c1->SaveAs("./outplots/RvsN_UPC.pdf");

}
