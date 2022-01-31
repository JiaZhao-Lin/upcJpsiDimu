#include "/Users/zaochenrice/myFunction.h"
#include "../../../common/constants.h"

void compSideBand()
{
	TFile* infile = new TFile("../../outdata/outfile_sideband_pt.root", "read");
	
	TH1D *hPt_LowMassBand_Rap[nRapBins], *hPt_HiMassBand_Rap[nRapBins];
	
	TH1D *hPt_LeftBand, *hPt_RightBand;
	TH1D *hPt_SideBand; 

	for( Int_t irap=0; irap<nRapBins; irap++ )
	{
		hPt_LowMassBand_Rap[irap] = (TH1D*) infile->Get( Form("hPt_LowMassBand_RapBin%d", irap) );
		hPt_HiMassBand_Rap[irap]  = (TH1D*) infile->Get( Form("hPt_HiMassBand_RapBin%d",  irap) );

		if( irap==0 )
		{
			hPt_LeftBand  = (TH1D*) hPt_LowMassBand_Rap[irap] ->Clone("hPt_LeftBand" );
			hPt_RightBand = (TH1D*) hPt_HiMassBand_Rap[irap]  ->Clone("hPt_RightBand");
		}
		else
		{
			hPt_LeftBand  -> Add( hPt_LowMassBand_Rap[irap] );
			hPt_RightBand -> Add( hPt_HiMassBand_Rap[irap]  );
		}

	}//irap

	hPt_SideBand  = (TH1D*) hPt_LeftBand -> Clone("hPt_SideBand");
	hPt_SideBand  -> Add( hPt_RightBand );

	hPt_LeftBand  -> Scale( 1./hPt_LeftBand ->Integral() );
	hPt_RightBand -> Scale( 1./hPt_RightBand->Integral() );
	hPt_SideBand  -> Scale( 1./hPt_SideBand ->Integral() );

	double x1= -0.01, x2= 0.50;
	//overlay
	gStyle->SetFrameLineWidth(3);
	TCanvas* c2 = new TCanvas("c2", "c2", 900, 900);
	c2 -> Divide(1,2);
	c2 -> cd(1);
	gPad->SetPad(0.0,0.25,1.0,0.96);
	gPad->SetBottomMargin(0);
	gPad->SetRightMargin(0.05);
	gPad->SetTopMargin(0);

	hPt_LeftBand ->GetYaxis()->SetTitleSize(0.05);
	hPt_LeftBand ->GetYaxis()->SetTitleOffset(1.25);
	hPt_LeftBand ->GetXaxis()->SetTitleSize(0.05);
	hPt_LeftBand ->GetXaxis()->SetTitleOffset(1.30);
	hPt_LeftBand ->SetMarkerStyle(20);
	hPt_LeftBand ->SetMarkerSize(0.9);
	hPt_LeftBand ->SetMarkerColor(1);
	hPt_LeftBand ->SetLineColor(1);
	hPt_LeftBand ->SetLineWidth(2);
	hPt_LeftBand ->SetAxisRange(x1,    x2,   "x");
	hPt_LeftBand ->SetAxisRange(3.e-4, 0.145, "y");
	hPt_LeftBand ->Draw("pe");
	hPt_LeftBand ->Draw("histsame");

	hPt_RightBand ->SetMarkerStyle(20);
	hPt_RightBand ->SetMarkerSize(0.9);
	hPt_RightBand ->SetMarkerColor(4);
	hPt_RightBand ->SetLineColor(4);
	hPt_RightBand ->SetLineWidth(2);
	hPt_RightBand ->Draw("pesame");
	hPt_RightBand ->Draw("histsame");
	
	hPt_SideBand ->SetMarkerStyle(20);
	hPt_SideBand ->SetMarkerSize(0.9);
	hPt_SideBand ->SetMarkerColor(2);
	hPt_SideBand ->SetLineColor(2);
	hPt_SideBand ->SetLineWidth(2);
	hPt_SideBand ->Draw("histsame");

	TLegend  *leg =  new TLegend(0.44, 0.55, 0.88, 0.80);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.052);
	leg->AddEntry( hPt_LeftBand,    "2.75<M_{#mu#mu}<2.90 GeV/c^{2}",        "lp");
	leg->AddEntry( hPt_RightBand,   "3.30<M_{#mu#mu}<3.45 GeV/c^{2}",        "lp");
	leg->AddEntry( hPt_SideBand,    "Left+Right SideBand",                   "lp");
	leg->Draw("same");

	drawLatex(0.45, 0.92, "p_{T} of #gamma#gamma #rightarrow #mu#mu in Data",       42,  0.07,  1);
	
	c2->cd(2);
	gPad->SetLogy(0);
	gPad->SetPad(0.0,0.0,1.0,0.25);
	gPad->SetTopMargin(0);
	gPad->SetRightMargin(0.05);
	gPad->SetBottomMargin(0.45);
	
	TH1D* HRatio_Right2Left = (TH1D*) hPt_RightBand->Clone("HRatio_Right2Left");
	HRatio_Right2Left -> Divide(hPt_LeftBand);
	
	HRatio_Right2Left ->SetMarkerStyle(24);
	HRatio_Right2Left ->SetMarkerSize(1.2);
	HRatio_Right2Left ->SetMarkerColor(1);
	HRatio_Right2Left ->SetLineColor(1);
	HRatio_Right2Left ->SetLineWidth(2);

	HRatio_Right2Left ->SetAxisRange(x1, x2, "x");
	HRatio_Right2Left ->GetYaxis()->SetRangeUser(0.01, 2.0);
	HRatio_Right2Left ->SetYTitle("#frac{RightBand}{LeftBand}");
	HRatio_Right2Left ->SetXTitle("p_{T} (GeV/c)");
	HRatio_Right2Left ->GetYaxis()->CenterTitle();
	HRatio_Right2Left ->GetYaxis()->SetNdivisions(6);
	HRatio_Right2Left ->GetYaxis()->SetTitleSize(0.12);
	HRatio_Right2Left ->GetYaxis()->SetTitleOffset(0.43);
	HRatio_Right2Left ->GetYaxis()->SetLabelSize(0.10);
	//HRatio_Right2Left ->GetYaxis()->SetLabelFont(20);
	HRatio_Right2Left ->GetXaxis()->SetTitleSize(0.20);
	HRatio_Right2Left ->GetXaxis()->SetTitleOffset(0.95);
	HRatio_Right2Left ->GetXaxis()->SetLabelSize(0.16);
	HRatio_Right2Left ->GetXaxis()->SetLabelFont(40);
	HRatio_Right2Left ->SetTickLength(0.08);
	HRatio_Right2Left ->Draw("histe");
	
	drawLine(x1,      1,  x2, 1,  2, 2, 4);
	
	c2->SaveAs("./outplots/hPt_SideBands_comp.png");
	c2->SaveAs("./outplots/hPt_SideBands_comp.pdf");

	//--------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------
	//read in the data file from the mass fit at low pt region
	TFile* infile_fromFit = new TFile("../outplots/outdata_pt_QED_CohJpsi.root", "read");
	TH1D* hPt_QED_fromFit = (TH1D*)infile_fromFit->Get("hPt_QED");
	
	hPt_QED_fromFit -> Scale(1./hPt_QED_fromFit->Integral());
	
	//overlay
	gStyle->SetFrameLineWidth(3);
	TCanvas* c1 = new TCanvas("c1", "c1", 900, 900);
	c1 -> Divide(1,2);
	c1 -> cd(1);
	gPad->SetPad(0.0,0.25,1.0,0.96);
	gPad->SetBottomMargin(0);
	gPad->SetRightMargin(0.05);
	gPad->SetTopMargin(0);

	x1= -0.01, x2= 0.50;

	hPt_SideBand ->GetYaxis()->SetTitleSize(0.05);
	hPt_SideBand ->GetYaxis()->SetTitleOffset(1.25);
	hPt_SideBand ->GetXaxis()->SetTitleSize(0.05);
	hPt_SideBand ->GetXaxis()->SetTitleOffset(1.30);
	hPt_SideBand ->SetMarkerStyle(20);
	hPt_SideBand ->SetMarkerSize(0.9);
	hPt_SideBand ->SetMarkerColor(2);
	hPt_SideBand ->SetLineColor(2);
	hPt_SideBand ->SetLineWidth(2);
	hPt_SideBand ->SetAxisRange(x1,    x2,    "x");
	hPt_SideBand ->SetAxisRange(3.e-4, 0.145, "y");
	hPt_SideBand ->Draw("pe");
	hPt_SideBand ->Draw("histsame");

	hPt_QED_fromFit ->SetMarkerStyle(20);
	hPt_QED_fromFit ->SetMarkerSize(0.9);
	hPt_QED_fromFit ->SetMarkerColor(4);
	hPt_QED_fromFit ->SetLineColor(4);
	hPt_QED_fromFit ->SetLineWidth(2);
	hPt_QED_fromFit ->Draw("pesame");
	hPt_QED_fromFit ->Draw("histsame");

	drawLatex(0.40, 0.90, "p_{T} of #gamma#gamma #rightarrow #mu#mu in Data",       42,  0.08,  1);
	
	TLegend  *leg2 =  new TLegend(0.44, 0.55, 0.88, 0.80);
	leg2->SetFillStyle(0);
	leg2->SetFillColor(0);
	leg2->SetTextFont(42);
	leg2->SetTextSize(0.052);
	leg2->AddEntry( hPt_QED_fromFit, "From Data Fit",      "lp");
	leg2->AddEntry( hPt_SideBand,    "From SideBand",      "lp");
	leg2->Draw("same");
	
	c1->cd(2);
	gPad->SetLogy(0);
	gPad->SetPad(0.0,0.0,1.0,0.25);
	gPad->SetTopMargin(0);
	gPad->SetRightMargin(0.05);
	gPad->SetBottomMargin(0.45);
	
	TH1D* HRatio_SideBand2DataFit = (TH1D*) hPt_QED_fromFit->Clone("HRatio_SideBand2DataFit");
	HRatio_SideBand2DataFit ->Reset();
	cout<<"make sure all avaliable bins of hPt_QED_fromFit are exactly same as hPt_SideBand !!!!!!!!! "<<endl;

	for(int ib=0; ib<HRatio_SideBand2DataFit->GetNbinsX(); ib++)
	{
		double ibContent = hPt_SideBand->GetBinContent(ib+1);
		double ibError   = hPt_SideBand->GetBinError(  ib+1);

		HRatio_SideBand2DataFit -> SetBinContent(ib+1, ibContent);
		HRatio_SideBand2DataFit -> SetBinError(  ib+1, ibError  );
	}

	HRatio_SideBand2DataFit ->Divide(hPt_QED_fromFit);

	HRatio_SideBand2DataFit ->SetMarkerStyle(24);
	HRatio_SideBand2DataFit ->SetMarkerSize(1.2);
	HRatio_SideBand2DataFit ->SetMarkerColor(1);
	HRatio_SideBand2DataFit ->SetLineColor(1);
	HRatio_SideBand2DataFit ->SetLineWidth(2);

	HRatio_SideBand2DataFit ->SetAxisRange(x1,    x2, "x");
	HRatio_SideBand2DataFit ->GetYaxis()->SetRangeUser(0.40, 1.8);
	HRatio_SideBand2DataFit ->SetYTitle("#frac{SideBand}{DataFit}");
	HRatio_SideBand2DataFit ->SetXTitle("p_{T} (GeV/c)");
	HRatio_SideBand2DataFit ->GetYaxis()->CenterTitle();
	HRatio_SideBand2DataFit ->GetYaxis()->SetNdivisions(6);
	HRatio_SideBand2DataFit ->GetYaxis()->SetTitleSize(0.12);
	HRatio_SideBand2DataFit ->GetYaxis()->SetTitleOffset(0.43);
	HRatio_SideBand2DataFit ->GetYaxis()->SetLabelSize(0.10);
	//HRatio_SideBand2DataFit ->GetYaxis()->SetLabelFont(20);
	HRatio_SideBand2DataFit ->GetXaxis()->SetTitleSize(0.20);
	HRatio_SideBand2DataFit ->GetXaxis()->SetTitleOffset(0.95);
	HRatio_SideBand2DataFit ->GetXaxis()->SetLabelSize(0.16);
	HRatio_SideBand2DataFit ->GetXaxis()->SetLabelFont(40);
	HRatio_SideBand2DataFit ->SetTickLength(0.08);
	HRatio_SideBand2DataFit ->Draw("histe");
	
	drawLine(x1,      1,  x2, 1,  2, 2, 4);
	
	c1->SaveAs("./outplots/hPt_QED_SideBandVsFromFit_comp.png");
	c1->SaveAs("./outplots/hPt_QED_SideBandVsFromFit_comp.pdf");

}
