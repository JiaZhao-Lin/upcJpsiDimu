#include "/Users/zaochenrice/myFunction.h"

void CompCohJpsiPt_DataVsSim( )
{
	TFile *inf_Temps = new TFile("../../../simulation/effAndTemp/MassPtTemp_AllSpecs.root", "read");
	TH1D  *hPt_CohJpsi_Sim = (TH1D *)inf_Temps->Get("hCohJpsiPt"          );

	//read in the data file from the mass fit at low pt region
	TFile* infile_fromFit = new TFile("../outplots/outdata_pt_QED_CohJpsi.root", "read");
	TH1D* hPt_Jpsi_fromFit = (TH1D*)infile_fromFit->Get("hPt_CohJpsi");//hPt_CohJpsi within Jpsi mass window, only pt<0.15 GeV is cohJpsi
	
	TH1D* hPt_CohJpsi_fromFit = new TH1D("hPt_CohJpsi_fromFit", "", 15, 0., 0.15);

	for( int ib=0; ib<hPt_CohJpsi_fromFit->GetNbinsX(); ib++ )
	{
		cout<<"ib: "<<ib+1<<" "<<hPt_CohJpsi_fromFit->GetBinCenter(ib+1)<<endl;

		hPt_CohJpsi_fromFit ->SetBinContent( ib+1, hPt_Jpsi_fromFit->GetBinContent(ib+1) );
		hPt_CohJpsi_fromFit ->SetBinError(   ib+1, hPt_Jpsi_fromFit->GetBinError(ib+1)   );
	}
	
	//Normalize for a fair comparison
	hPt_CohJpsi_fromFit -> Scale(1./hPt_CohJpsi_fromFit->Integral()); //Normalized within 0-0.15 GeV/c 
	
	const int temLBin = hPt_CohJpsi_Sim->FindBin(0.00);
	const int temHBin = hPt_CohJpsi_Sim->FindBin(0.15);

	cout<<"bins to normalize the CohJpsi in Sim"<<endl;
	cout<<"temLBin: "<<temLBin<<endl;
	cout<<"temHBin: "<<temHBin<<endl;
	
	hPt_CohJpsi_Sim   -> Scale(1./hPt_CohJpsi_Sim->Integral(temLBin, temHBin));
	
	//overlay
	gStyle->SetFrameLineWidth(3);
	TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
	c1 -> cd();
	setPad(0.13, 0.05, 0.05, 0.13);

	gPad->SetLogy(0);

	double x1= -0.01, x2= 0.15;

	hPt_CohJpsi_Sim ->GetYaxis()->SetTitleSize(0.05);
	hPt_CohJpsi_Sim ->GetYaxis()->SetTitleOffset(1.15);
	hPt_CohJpsi_Sim ->GetXaxis()->SetTitleSize(0.05);
	hPt_CohJpsi_Sim ->GetXaxis()->SetTitleOffset(1.10);
	hPt_CohJpsi_Sim ->SetMarkerStyle(20);
	hPt_CohJpsi_Sim ->SetMarkerSize(0.9);
	hPt_CohJpsi_Sim ->SetMarkerColor(2);
	hPt_CohJpsi_Sim ->SetLineColor(2);
	hPt_CohJpsi_Sim ->SetLineWidth(2);
	hPt_CohJpsi_Sim ->SetAxisRange(x1,    x2+0.4,    "x");
	hPt_CohJpsi_Sim ->SetAxisRange(3.e-6, 0.135,     "y");
	hPt_CohJpsi_Sim ->SetYTitle("");
	hPt_CohJpsi_Sim ->Draw("pe");
	hPt_CohJpsi_Sim ->Draw("histsame");

	hPt_CohJpsi_fromFit ->SetMarkerStyle(20);
	hPt_CohJpsi_fromFit ->SetMarkerSize(0.9);
	hPt_CohJpsi_fromFit ->SetMarkerColor(4);
	hPt_CohJpsi_fromFit ->SetLineColor(4);
	hPt_CohJpsi_fromFit ->SetLineWidth(2);
	hPt_CohJpsi_fromFit ->Draw("pesame");
	hPt_CohJpsi_fromFit ->Draw("histsame");

	drawLatex(0.30, 0.88, "p_{T} of Coherent J/#psi #rightarrow #mu#mu",  42,  0.06,  1);
	
	TLegend  *leg2 =  new TLegend(0.44, 0.55, 0.88, 0.80);
	leg2->SetFillStyle(0);
	leg2->SetFillColor(0);
	leg2->SetTextFont(42);
	leg2->SetTextSize(0.052);
	leg2->AddEntry( hPt_CohJpsi_fromFit,    "From Data Fit",     "lp");
	leg2->AddEntry( hPt_CohJpsi_Sim,        "From Simulation",   "lp");
	leg2->Draw("same");

	c1->SaveAs("./outplots/hPt_CohJpsi_FromSimVsFromDataFit_comp.png");
	c1->SaveAs("./outplots/hPt_CohJpsi_FromSimVsFromDataFit_comp.pdf");
	
	//-----------------------------------------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------------------------------
	//Normalize simulated template to have the same height as the pt shape from Data Fit
	TH1D* hPt_CohJpsi_Sim_SameHeight = (TH1D*) hPt_CohJpsi_Sim -> Clone("hPt_CohJpsi_Sim_SameHeight");
	hPt_CohJpsi_Sim_SameHeight -> Scale( hPt_CohJpsi_fromFit->GetBinContent(5) / hPt_CohJpsi_Sim->GetBinContent(6) );
	
	hPt_CohJpsi_Sim_SameHeight ->SetAxisRange(x1,    x2+0.4,    "x");
	hPt_CohJpsi_Sim_SameHeight ->SetAxisRange(3.e-6, 0.135,     "y");
	hPt_CohJpsi_Sim_SameHeight ->Draw("pe");
	hPt_CohJpsi_Sim_SameHeight ->Draw("histsame");
	
	hPt_CohJpsi_fromFit ->Draw("pesame");
	hPt_CohJpsi_fromFit ->Draw("histsame");

	drawLatex(0.30, 0.88, "p_{T} of Coherent J/#psi #rightarrow #mu#mu",  42,  0.06,  1);
	
	TLegend  *leg3 =  new TLegend(0.44, 0.55, 0.88, 0.80);
	leg3->SetFillStyle(0);
	leg3->SetFillColor(0);
	leg3->SetTextFont(42);
	leg3->SetTextSize(0.052);
	leg3->AddEntry( hPt_CohJpsi_fromFit,        "From Data Fit",     "lp");
	leg3->AddEntry( hPt_CohJpsi_Sim_SameHeight, "From Simulation",   "lp");
	leg3->Draw("same");

	c1->SaveAs("./outplots/hPt_CohJpsi_FromSimVsFromDataFit_SameHeight_comp.png");
	c1->SaveAs("./outplots/hPt_CohJpsi_FromSimVsFromDataFit_SameHeight_comp.pdf");
	
	//-----------------------------------------------------------------------------------------------------------------
	// Shift the simulated data to left by 0.01 GeV/c and compare
	//-----------------------------------------------------------------------------------------------------------------
	c1->SetLogy(0);

	TH1D* hPt_CohJpsi_Sim_Shifted = (TH1D*) hPt_CohJpsi_Sim_SameHeight->Clone("hPt_CohJpsi_Sim_Shifted");

	for(int ib=1; ib<hPt_CohJpsi_Sim_Shifted->GetNbinsX(); ib++)
	{
		double ibContent = hPt_CohJpsi_Sim_SameHeight->GetBinContent( ib+1 );
		double ibError   = hPt_CohJpsi_Sim_SameHeight->GetBinError(   ib+1 );

		hPt_CohJpsi_Sim_Shifted ->SetBinContent( ib, ibContent );
		hPt_CohJpsi_Sim_Shifted ->SetBinError(   ib, ibError   );
	}

	hPt_CohJpsi_Sim_Shifted ->SetAxisRange(x1,    x2+0.4,   "x");
	hPt_CohJpsi_Sim_Shifted ->SetAxisRange(3.e-6, 0.135,     "y");
	hPt_CohJpsi_Sim_Shifted ->Draw("pe");
	hPt_CohJpsi_Sim_Shifted ->Draw("histsame");
	
	hPt_CohJpsi_fromFit ->Draw("pesame");
	hPt_CohJpsi_fromFit ->Draw("histsame");

	drawLatex(0.30, 0.88, "p_{T} of Coherent J/#psi #rightarrow #mu#mu",  42,  0.06,  1);
	
	TLegend  *leg4 =  new TLegend(0.40, 0.55, 0.88, 0.80);
	leg4->SetFillStyle(0);
	leg4->SetFillColor(0);
	leg4->SetTextFont(42);
	leg4->SetTextSize(0.052);
	leg4->AddEntry( hPt_CohJpsi_fromFit,       "From Data Fit",        "lp" );
	leg4->AddEntry( hPt_CohJpsi_Sim_Shifted,   "From Simu. Shifted",   "lp" );
	leg4->Draw("same");

	c1->SaveAs("./outplots/hPt_CohJpsi_FromSimVsFromDataFit_SameHeight_withShift.png");
	c1->SaveAs("./outplots/hPt_CohJpsi_FromSimVsFromDataFit_SameHeight_withShift.pdf");

	//save the shifted template
	TFile* outfile = new TFile("./outplots/outTemp_CohJpsi_shifted.root", "recreate");
	cout<<"shifted coherent template is save into: "<<outfile->GetName()<<endl;
	outfile ->cd();
	hPt_CohJpsi_Sim_Shifted -> Write();

	outfile -> Close();

	//Now let's find the stretch factor by matching the streched tempalte to what in real data
}
