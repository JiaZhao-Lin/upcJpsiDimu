#include "../common/LoadSignal.C"


void getEffRatio(){
	LoadEfficiency Eff 		("out4effAndTemp/Efficiency_AllSpecs_6RapBins.root", 			false);
	LoadEfficiency Eff_TnP	("out4effAndTemp/Efficiency_AllSpecs_6RapBins.appliedTnP.root", false);
	
	TH1D* hEffVsY = (TH1D*) Eff.H_EffVsY_CohJpsi->Clone();
	TH1D* hEffVsY_TnP = (TH1D*) Eff_TnP.H_EffVsY_CohJpsi->Clone();
	TH1D* hEffVsY_Ratio = (TH1D*) Eff_TnP.H_EffVsY_CohJpsi->Clone();
	hEffVsY_Ratio->Divide( hEffVsY);


	TH2D* EffVsY = new TH2D("EffVsY", ";Rapidity;Efficiency;", 10,-2.4,2.4, 10, 0, 1.1);
	EffVsY->Draw();
	

	hEffVsY->Draw("Same");

	hEffVsY_TnP->SetMarkerStyle(24);
	hEffVsY_TnP->SetMarkerColor(2);
	hEffVsY_TnP->SetLineColor(2);
	hEffVsY_TnP->Draw("Same");

	hEffVsY_Ratio->SetMarkerStyle(3);
	hEffVsY_Ratio->Draw("Same");

	TLegend  *leg =  new TLegend(0.4, 0.4, 0.6, 0.6);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.04);
	leg->AddEntry(hEffVsY,               "Eff_MC",                  "p");
	leg->AddEntry(hEffVsY_TnP,           "Eff_TnP",                 "p");
	leg->AddEntry(hEffVsY_Ratio,         "TnP/MC",                   "p");
	leg->Draw("same");
}