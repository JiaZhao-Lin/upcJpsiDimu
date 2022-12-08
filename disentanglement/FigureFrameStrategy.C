#include "../common/function.C"

enum class FigureFrameStrategyList
{
	Sigma_Log,
	Sigma_LogLog,
	R,
	DSigmaDy,
	DSigmaDy_NeuConfig
};


struct FigureFrameStrategy
{
	std::unique_ptr< TCanvas > c = std::make_unique< TCanvas >	();
	std::unique_ptr< TH1 > htemp;

	virtual void Apply() = 0;
	virtual ~FigureFrameStrategy() = default;

	void SaveAs(TString fileName){	c->SaveAs(fileName);	};
};

struct Sigma_Log_Strategy : FigureFrameStrategy
{

	void Apply()
	{
		c->SetLogy();
		gPad->SetTopMargin(0.08);
		gPad->SetBottomMargin(0.12);
		gPad->SetRightMargin(0.05);

		htemp = std::make_unique< TH2D >	("SigmaVsW", ";W_{#gammaN}^{Pb} (GeV);#sigma(#gamma Pb #rightarrow J/#psi Pb) (mb);", 10,0,420, 10, 0.004, 0.19);
		htemp->GetYaxis()->CenterTitle();
		htemp->GetXaxis()->CenterTitle();
		htemp->GetYaxis()->SetTitleSize(0.05);
		htemp->GetYaxis()->SetTitleSize(0.05);
		htemp->GetYaxis()->SetTitleOffset(0.99);
		htemp->GetYaxis()->SetLabelSize(0.04);
		htemp->GetXaxis()->SetTitleSize(0.05);
		htemp->GetXaxis()->SetTitleOffset(0.98);
		htemp->GetXaxis()->SetLabelSize(0.04);
		htemp->SetTickLength(0.04);
		htemp->Draw("");

		drawLatex(0.18,0.86,"x~3.8#times10^{-3}",42, 0.025, 1);
		drawLatex(0.29,0.86,"x~9.6#times10^{-4}",42, 0.025, 1);
		drawLatex(0.48,0.86,"x~2.4#times10^{-4}",42, 0.025, 1);
		drawLatex(0.67,0.86,"x~1.1#times10^{-4}",42, 0.025, 1);
		drawLatex(0.87,0.86,"x~6.0#times10^{-5}",42, 0.025, 1);

		drawLatex(0.15, 0.80, "Pb+Pb #rightarrow Pb'+Pb'+J/#psi",  42,        0.05,      1 );
		drawLatex(0.56,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
		drawLatex(0.15,0.94,"#bf{CMS} #it{Preliminary}",42, 0.05, 1);

	}
};

struct Sigma_LogLog_Strategy : FigureFrameStrategy
{
	void Apply()
	{
		c->SetLogx();
		c->SetLogy();
		gPad->SetTopMargin(0.08);
		gPad->SetBottomMargin(0.12);
		gPad->SetRightMargin(0.05);

		htemp = std::make_unique< TH2D >	("SigmaVsW", ";W_{#gammaN}^{Pb} (GeV);#sigma(#gamma Pb #rightarrow J/#psi Pb) (mb);", 10,14,520, 10, 0.004, 0.13);
		htemp->GetYaxis()->CenterTitle();
		htemp->GetXaxis()->CenterTitle();
		htemp->GetYaxis()->SetTitleSize(0.05);
		htemp->GetYaxis()->SetTitleSize(0.05);
		htemp->GetYaxis()->SetTitleOffset(0.99);
		htemp->GetYaxis()->SetLabelSize(0.04);
		htemp->GetXaxis()->SetTitleSize(0.05);
		htemp->GetXaxis()->SetTitleOffset(0.98);
		htemp->GetXaxis()->SetLabelSize(0.04);
		htemp->SetTickLength(0.04);
		htemp->Draw("");

		drawLatex(0.18,0.86,"x~3.8#times10^{-3}",42, 0.025, 1);
		drawLatex(0.29,0.86,"x~9.6#times10^{-4}",42, 0.025, 1);
		drawLatex(0.48,0.86,"x~2.4#times10^{-4}",42, 0.025, 1);
		drawLatex(0.67,0.86,"x~1.1#times10^{-4}",42, 0.025, 1);
		drawLatex(0.87,0.86,"x~6.0#times10^{-5}",42, 0.025, 1);

		drawLatex(0.15, 0.80, "Pb+Pb #rightarrow Pb'+Pb'+J/#psi",  42,        0.05,      1 );
		drawLatex(0.56,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
		drawLatex(0.15,0.94,"#bf{CMS} #it{Preliminary}",42, 0.05, 1);

	}
};

struct R_Strategy : FigureFrameStrategy
{
	void Apply()
	{
		c->SetLogx();
		gPad->SetTopMargin(0.08);
		gPad->SetBottomMargin(0.12);
		//gPad->SetLeftMargin(0.05);
		gPad->SetRightMargin(0.05);

		//htemp = std::make_unique< TH2D >	("htemp", "", 10,3.0e-5,5.0e-2, 10, 0, 1);
		htemp = std::make_unique< TH2D >	("RvsX", ";R^{Pb}_{g}(x, #mu^{2}=2.4 GeV^{2});x;", 10,1.0e-5,5e-2, 10, 0.2, 1.05);
		htemp->GetYaxis()->SetTitleSize(0.06);
		htemp->GetYaxis()->SetTitleOffset(0.85);
		htemp->GetYaxis()->SetLabelSize(0.04);
		htemp->GetYaxis()->CenterTitle();
		htemp->GetXaxis()->CenterTitle();
		htemp->GetXaxis()->SetTitleSize(0.07);
		htemp->GetXaxis()->SetTitleOffset(0.69);
		htemp->GetXaxis()->SetLabelSize(0.04);
		htemp->SetTickLength(0.04);
		htemp->Draw();

		drawLatex(0.15, 0.84, "Pb+Pb #rightarrow Pb'+Pb'+J/#psi",  42,        0.05,      1 );
		drawLatex(0.56,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
		drawLatex(0.15,0.94,"#bf{CMS} #it{Preliminary}",42, 0.05, 1);

	}
};

struct DSigmaDy_Strategy : FigureFrameStrategy
{
	void Apply()
	{
		setPad(0.12, 0.08, 0.07, 0.13);

		//htemp = std::make_unique< TH2D >	("htemp", "", 10, -4.1, 1.0, 10, 0, 14.0);
		htemp = std::make_unique< TH2D >	("AnAn", "AnAn;y;d#sigma_{J/#psi}/dy (mb);", 10, -4.5, 0, 10, 0, 8.0);
		htemp ->GetYaxis()->CenterTitle();
		//htemp ->GetYaxis()->SetNdivisions(6);
		htemp ->GetYaxis()->SetTitleSize(0.07);
		htemp ->GetYaxis()->SetTitleOffset(0.63);
		htemp ->GetYaxis()->SetLabelSize(0.06);
		htemp ->GetXaxis()->CenterTitle();
		htemp ->GetXaxis()->SetTitleSize(0.07);
		htemp ->GetXaxis()->SetTitleOffset(0.75);
		htemp ->GetXaxis()->SetLabelSize(0.06);
		//htemp ->SetTickLength(0.08);
		htemp->Draw();

		drawLatex(0.15, 0.85, "Pb+Pb #rightarrow Pb+Pb+J/#psi #sqrt{s_{NN}} = 5.02 TeV",      42,       0.06,      1);

	}
};

struct DSigmaDy_mini_Strategy : FigureFrameStrategy
{
	TString frameName;
	DSigmaDy_mini_Strategy(TString frameName_):	frameName{frameName_}	{};

	void Apply()
	{
		TString n = "htemp_" + frameName;
		htemp = std::make_unique< TH2D >	(n, n + ";y;d#sigma_{J/#psi}/dy (mb)", 10, -4.1, 0, 10, 0.0, 7.0);
		htemp ->GetYaxis()->CenterTitle();
		htemp ->GetYaxis()->SetNdivisions(6);
		htemp ->GetXaxis()->SetNdivisions(6);
		htemp ->GetYaxis()->SetTitleSize(0.065);
		htemp ->GetYaxis()->SetTitleOffset(0.8);
		htemp ->GetYaxis()->SetLabelSize(0.05);
		htemp ->GetXaxis()->CenterTitle();
		htemp ->GetXaxis()->SetTitleSize(0.065);
		htemp ->GetXaxis()->SetTitleOffset(0.75);
		htemp ->GetXaxis()->SetLabelSize(0.05);
		//htemp ->SetTickLength(0.08);
		htemp->Draw();

		drawLatex(0.15, 0.85, "Pb+Pb #rightarrow Pb+Pb+J/#psi",      42,       0.053,      1);
		drawLatex(0.8, 0.8, Form("#bf{%s}", frameName.Data()),      42,       0.075,      1);
		drawLatex(0.58,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
		drawLatex(0.15,0.94,"#bf{CMS} #it{Preliminary}",42, 0.05, 1);
	}
};

struct DSigmaDy_NeuConfig_Strategy : FigureFrameStrategy
{
	std::vector< std::unique_ptr< DSigmaDy_mini_Strategy > > miniFrames;
	std::vector< TString > frameNames	=	{"0n0n",	"0nXn",	"XnXn",	"AnAn"};
	DSigmaDy_NeuConfig_Strategy()
	{
		for (int i = 0; i < frameNames.size(); ++i)
		{
			miniFrames.push_back( std::make_unique<DSigmaDy_mini_Strategy> (	frameNames[i]	) );
		}
	}

	void Apply()
	{
		c->Divide(2,2);
		setPad(0.12, 0.08, 0.07, 0.13);

		for (int i = 1; i < 5; ++i)
		{
			c->cd(i);
			c->GetPad(i)->SetRightMargin(0.01);

			miniFrames[i]->Apply();
		}
	}
};

// void FigureFrameStrategy()
// {
//     auto c = new TCanvas();
//     Sigma_Log_Strategy s;
//  void Apply();
// }