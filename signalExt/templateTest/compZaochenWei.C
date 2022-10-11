const TString datanames[2] = {"~/Downloads/pttemplate_new_R77_Pt05_comb.root", "pttemplate_new_R77_Pt05_zy1.root"};

const TString histnames[3] = {"hPt_1619_Combo", "hPt_1921_Combo", "hPt_2124_Combo"};

TH1D* h1d_pt[2][3];


void compZaochenWei()
{
	for(int id=0; id<2; id++)
	{
		TFile* infile = new TFile(datanames[id], "read");

		for(int ih=0; ih<3; ih++)
		{
			h1d_pt[id][ih] = (TH1D*) infile->Get( histnames[ih] );
		}
	}

	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);

	const int colors[2] = {1,4};
	const TString titles[2] = {"WeiDrawfromTree", "ZaochenDrawfromTree"};
	for(int ih=0; ih<3; ih++)
	{
		c1->cd();
		c1->SetLogy(1);
		for(int id=0; id<2; id++)
		{
			h1d_pt[id][ih] -> SetLineColor(colors[id]);
			h1d_pt[id][ih] ->Scale(1./h1d_pt[id][ih]->Integral());
				
			if(id==0)
			{
				h1d_pt[id][ih] -> SetAxisRange(0.,1.0, "x");
				h1d_pt[id][ih] -> Draw("hist");
			}
			else
			{
				h1d_pt[id][ih] -> Draw("histsame");
			}
		}

		TLegend  *leg =  new TLegend(0.52, 0.65, 0.88, 0.88);
		leg->SetFillStyle(0);
		leg->SetFillColor(0);
		leg->SetTextSize(0.035);

		for(int id=0; id<2; id++)
		{
			leg->AddEntry(h1d_pt[id][ih],    titles[id],           "l");
		}
		leg->Draw("same");

		c1->SaveAs("./plots/comp_"+histnames[ih]+".png");
	}

}

