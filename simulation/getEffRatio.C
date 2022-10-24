#include "../common/LoadSignal.C"


void getEffRatio(){
	const int template_option  = 1; //0:OldCohJpsi; 1:NewCohJpsi (w R+1fm);
	const TString template_Name[2] = {"", "_NewCohJpsi"};
	const TString template_Dir[2] = {"out4effAndTemp", "out4effAndTemp_NewCohJpsi"};

	LoadEfficiency Eff 		(Form("%s/Efficiency_AllSpecs_6RapBins.root", template_Dir[template_option].Data()), 			false);
	LoadEfficiency Eff_TnP	(Form("%s/Efficiency_AllSpecs_6RapBins.appliedTnP.root", template_Dir[template_option].Data()), false);
	
	TH1D* hEffVsY = (TH1D*) Eff.H_EffVsY_CohJpsi->Clone();
	TH1D* hEffVsY_TnP = (TH1D*) Eff_TnP.H_EffVsY_CohJpsi->Clone();
	TH1D* hEffVsY_Ratio = (TH1D*) Eff_TnP.H_EffVsY_CohJpsi->Clone();
	hEffVsY_Ratio->Divide( hEffVsY);

	TCanvas* c = new TCanvas();

	TH2D* EffVsY = new TH2D("EffVsY", ";Rapidity;Efficiency;", 10,-2.4,2.4, 10, 0, 1.2);
	EffVsY->Draw();
	
	for (int i = 0; i < hEffVsY->GetNbinsX(); ++i)
	{
		auto a = hEffVsY->GetBinContent(i+1);
		auto b = hEffVsY_TnP->GetBinContent(i+1);
		cout<<a<<" "<<b<<" "<<b/a<<endl;
	}

	hEffVsY->Draw("Same");

	hEffVsY_TnP->SetMarkerStyle(24);
	hEffVsY_TnP->SetMarkerColor(2);
	hEffVsY_TnP->SetLineColor(2);
	hEffVsY_TnP->Draw("Same");

	hEffVsY_Ratio->SetMarkerStyle(3);
	hEffVsY_Ratio->SetMarkerColor(1);
	hEffVsY_Ratio->SetLineColor(1);
	hEffVsY_Ratio->Draw("Same");

	TLegend  *leg =  new TLegend(0.4, 0.4, 0.6, 0.6);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.04);
	leg->AddEntry(hEffVsY,               "Eff_MC",                  "p");
	leg->AddEntry(hEffVsY_TnP,           "Eff_TnP",                 "p");
	leg->AddEntry(hEffVsY_Ratio,         "TnP/MC",                   "p");
	leg->Draw("same");

	c->SaveAs("./out4effAndTemp/Eff_Ratio.pdf");
	c->SaveAs("./out4effAndTemp/Eff_Ratio.png");
}