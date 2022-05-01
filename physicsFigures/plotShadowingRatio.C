#include "./ImpulseApproximation/runUPC_AAModel.C"
#include "fit2D.C"

const double JpsiMass   = 3.096916;
const double Sqrt_s     = 5020;
const double Gamma_beam = 2672.9;
const double Mass_N     = (0.93827+0.93957) / 2;

std::vector<double> Raps,      Raps_Err;
std::vector<double> Sigmas,    Sigmas_Err;
std::vector<double> Sigmas_IA, Sigmas_IA_Err;
std::vector<double> R,         R_Err;

std::vector<double> 	Xs = {},  Xs_Err = {};
std::vector<double> 	Ws = {},  Ws_Err = {};

void Cal_R_Error()
{
	for (int i = 0; i < Raps.size(); ++i)
	{
		double tem = R[i]*R[i] * sqrt(pow(Sigmas_Err[i]/Sigmas[i],2) + pow(Sigmas_IA_Err[i]/Sigmas_IA[i],2) );
		R_Err.push_back( tem/R[i]/2);
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

void plotSigmaVsW()
{
	auto c = new TCanvas();
	TH2D* SigmaVsW = new TH2D("SigmaVsW", ";W_{#gamma p} (GeV);#sigma_{#gamma A #rightarrow J/#Psi A} (mb)", 10,0,450, 10, 0, 0.15);

	SigmaVsW->GetYaxis()->SetTitleSize(0.065);
	SigmaVsW->GetYaxis()->SetTitleOffset(0.85);
	SigmaVsW->GetYaxis()->SetLabelSize(0.04);
	SigmaVsW->GetXaxis()->SetTitleSize(0.05);
	SigmaVsW->GetXaxis()->SetTitleOffset(0.89);
	SigmaVsW->GetXaxis()->SetLabelSize(0.04);
	SigmaVsW->SetTickLength(0.04);
	SigmaVsW->Draw();

	TGraphErrors* ge_CMS        		= new TGraphErrors(Ws.size(),	&Ws[0],	&Sigmas[0],	0,	&Sigmas_Err[0]	);
	TGraphErrors* ge_ALICE_Run2_MidRap	= new TGraphErrors(ALICE_Run2_MidRap_W.size(),	&ALICE_Run2_MidRap_W[0],	&ALICE_Run2_MidRap_Sigma[0],	0,	&ALICE_Run2_MidRap_Sigma_Err[0]	);

	TGraphErrors* ge_CGCnoFluct = new TGraphErrors(CGC_JpsiNoFluct_W.size(),	&CGC_JpsiNoFluct_W[0],	&CGC_JpsiNoFluct_CohXsec[0],	0,	0);
	TGraphErrors* ge_IA        	= new TGraphErrors(W_IA.size(),	&W_IA[0],	&Sigma_IA[0],	0,	0);

	ge_CMS->SetMarkerStyle(20);
	ge_CMS->SetMarkerColor(1);
	ge_CMS->SetLineColor(1);
	ge_CMS->SetLineWidth(2);
	ge_CMS->Draw("Pesame");

	ge_ALICE_Run2_MidRap->SetMarkerStyle(24);
	ge_ALICE_Run2_MidRap->SetMarkerColor(4);
	ge_ALICE_Run2_MidRap->SetLineColor(4);
	ge_ALICE_Run2_MidRap->SetLineWidth(2);
	ge_ALICE_Run2_MidRap->Draw("Pesame");

	ge_CGCnoFluct->SetMarkerColor(4);
	ge_CGCnoFluct->SetLineColor(4);
	ge_CGCnoFluct->SetLineWidth(2);
	ge_CGCnoFluct->Draw("csame");

	ge_IA->SetMarkerColor(2);
	ge_IA->SetLineColor(2);
	ge_IA->SetLineWidth(2);
	ge_IA->Draw("csame");

	TLegend  *leg =  new TLegend(0.15, 0.65, 0.48, 0.90);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.055);
	leg->AddEntry(ge_CMS,	"CMS",			"lp");
	leg->AddEntry(ge_ALICE_Run2_MidRap,	"ALICE MidRap",			"lp");
	leg->AddEntry(ge_CGCnoFluct,	"CGC NoFluct",	"l");
	leg->AddEntry(ge_IA,	"IA",	        "l");
	leg->Draw("same");

	c->SaveAs("outplots/SigmaVsW.png");
	c->SaveAs("outplots/SigmaVsW.pdf");
}

void plotShadowingRatio( )
{
	{auto Temp = fit2D(); Raps = Temp[0]; Raps_Err = Temp[1]; Sigmas = Temp[2]; Sigmas_Err = Temp[3];}
	for (int i = 0; i < Raps.size(); ++i)
	{
		Xs 		.push_back( y2x(Raps[i]) 	);
		Xs_Err	.push_back( 0 );
		Ws 		.push_back( y2W(Raps[i]) 	);
		Ws_Err 	.push_back( 0 );

		cout<<Form("y: %f,	x: %f,	W: %f", Raps[i], Xs[i], Ws[i])<<endl;
		cout<<"x from W: "<<pow(JpsiMass,2)/pow(Ws[i],2)<<endl;
	}
	
	{auto Temp = runUPC_AAModel(Ws); Sigmas_IA = Temp[0]; Sigmas_IA_Err = Temp[1];}

	for (int i = 0; i < Raps.size(); ++i)
	{
		R.push_back( sqrt( Sigmas[i]/Sigmas_IA[i] ) );
	}
	
	Cal_R_Error();
	plotSigmaVsW();

	auto c = new TCanvas();
	c->SetLogx();

	TH2D* htem2d = new TH2D("htem2d", "", 10,3.0e-5,5.0e-2, 10, 0, 1);
	htem2d->SetYTitle("R^{Pb}_{g}(x, #mu^{2}=0.24 GeV^{2})");
	htem2d->SetXTitle("x");

	htem2d->GetYaxis()->SetTitleSize(0.06);
	htem2d->GetYaxis()->SetTitleOffset(0.85);
	htem2d->GetYaxis()->SetLabelSize(0.04);
	htem2d->GetXaxis()->SetTitleSize(0.07);
	htem2d->GetXaxis()->SetTitleOffset(0.60);
	htem2d->GetXaxis()->SetLabelSize(0.04);
	htem2d->SetTickLength(0.04);
	htem2d->Draw();

	TGraphErrors* ge_CMS        	= new TGraphErrors(Xs.size(),           &Xs[0],           &R[0],            &Xs_Err[0], &R_Err[0]       );
	TGraphErrors* ge_ALICE_Run1 	= new TGraphErrors(ALICE_x.size(),      &ALICE_x[0],      &ALICE_R[0],      0,          &ALICE_R_Err[0] );
	TGraphErrors* ge_ALICE_Run2 	= new TGraphErrors(ALICE_Run2_x.size(), &ALICE_Run2_x[0], &ALICE_Run2_R[0], 0,          &ALICE_Run2_R_Err[0] );
	TGraphErrors* ge_ALICE_Run2_Psi = new TGraphErrors(ALICE_Run2_Psi_x.size(), &ALICE_Run2_Psi_x[0], &ALICE_Run2_Psi_R[0], 0,          &ALICE_Run2_R_Err[0] );

	ge_CMS->SetMarkerStyle(20);
	ge_CMS->SetMarkerColor(1);
	ge_CMS->SetLineColor(1);
	ge_CMS->SetLineWidth(2);
	ge_CMS->Draw("Pesame");

	ge_ALICE_Run1->SetMarkerStyle(24);
	ge_ALICE_Run1->SetMarkerColor(4);
	ge_ALICE_Run1->SetLineColor(4);
	ge_ALICE_Run1->SetLineWidth(2);
	//ge_ALICE_Run1->Draw("pesame");
	ge_ALICE_Run2->SetMarkerStyle(24);
	ge_ALICE_Run2->SetMarkerColor(4);
	ge_ALICE_Run2->SetLineColor(4);
	ge_ALICE_Run2->SetLineWidth(2);
	ge_ALICE_Run2->Draw("pesame");

	ge_ALICE_Run2_Psi->SetMarkerStyle(24);
	ge_ALICE_Run2_Psi->SetMarkerColor(2);
	ge_ALICE_Run2_Psi->SetLineColor(2);
	ge_ALICE_Run2_Psi->SetLineWidth(2);
	ge_ALICE_Run2_Psi->Draw("pesame");
	
	drawLatex(0.15, 0.84, "Nuclear suppresion factor", 42, 0.06, 1);
	
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

	TLegend  *leg =  new TLegend(0.55, 0.25, 0.88, 0.50);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.055);
	leg->AddEntry(ge_CMS,         "CMS",        "lp");
	//leg->AddEntry(ge_ALICE_Run1,  "ALICE Run1", "lp");
	leg->AddEntry(ge_ALICE_Run2,  "ALICE Run2", "lp");
	leg->AddEntry(ge_ALICE_Run2_Psi,  "ALICE Run2 Psi", "lp");
	leg->Draw("same");

	c->SaveAs("outplots/ShadowingRatiovsX.png");
	c->SaveAs("outplots/ShadowingRatiovsX.pdf");
}
