#include "../common/function.h"

enum class FigureFrameStrategyList
{
	Sigma_Log,
	Sigma_LogLog,
	R,
	DSigmaDy
};


struct FigureFrameStrategy
{
	virtual TH1* Apply() = 0;
	virtual ~FigureFrameStrategy() = default;
};

struct Sigma_Log_Strategy : FigureFrameStrategy
{
	Sigma_Log_Strategy(){};
	TH1* Apply()
	{
		gPad->SetTopMargin(0.08);
		gPad->SetBottomMargin(0.12);
		gPad->SetRightMargin(0.05);

		TH2D* hist = new TH2D("SigmaVsW", ";W_{#gammaN}^{Pb} (GeV);#sigma(#gamma Pb #rightarrow J/#psi Pb) (mb);", 10,0,420, 10, 0.004, 0.19);
		hist->GetYaxis()->CenterTitle();
		hist->GetXaxis()->CenterTitle();
		hist->GetYaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->SetTitleOffset(0.99);
		hist->GetYaxis()->SetLabelSize(0.04);
		hist->GetXaxis()->SetTitleSize(0.05);
		hist->GetXaxis()->SetTitleOffset(0.98);
		hist->GetXaxis()->SetLabelSize(0.04);
		hist->SetTickLength(0.04);
		hist->Draw("");

		drawLatex(0.18,0.86,"x~3.8#times10^{-3}",42, 0.025, 1);
		drawLatex(0.29,0.86,"x~9.6#times10^{-4}",42, 0.025, 1);
		drawLatex(0.48,0.86,"x~2.4#times10^{-4}",42, 0.025, 1);
		drawLatex(0.67,0.86,"x~1.1#times10^{-4}",42, 0.025, 1);
		drawLatex(0.87,0.86,"x~6.0#times10^{-5}",42, 0.025, 1);

		drawLatex(0.15, 0.80, "Pb+Pb #rightarrow Pb'+Pb'+J/#psi",  42,        0.05,      1 );
		drawLatex(0.56,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
		drawLatex(0.15,0.94,"#bf{CMS} #it{Preliminary}",42, 0.05, 1);

		return hist;
	}
};

struct Sigma_LogLog_Strategy : FigureFrameStrategy
{
	TH1* Apply()
	{
		gPad->SetTopMargin(0.08);
		gPad->SetBottomMargin(0.12);
		gPad->SetRightMargin(0.05);

		TH2D* hist = new TH2D("SigmaVsW", ";W_{#gammaN}^{Pb} (GeV);#sigma(#gamma Pb #rightarrow J/#psi Pb) (mb);", 10,14,520, 10, 0.004, 0.13);
		hist->GetYaxis()->CenterTitle();
		hist->GetXaxis()->CenterTitle();
		hist->GetYaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->SetTitleOffset(0.99);
		hist->GetYaxis()->SetLabelSize(0.04);
		hist->GetXaxis()->SetTitleSize(0.05);
		hist->GetXaxis()->SetTitleOffset(0.98);
		hist->GetXaxis()->SetLabelSize(0.04);
		hist->SetTickLength(0.04);
		hist->Draw("");

		drawLatex(0.18,0.86,"x~3.8#times10^{-3}",42, 0.025, 1);
		drawLatex(0.29,0.86,"x~9.6#times10^{-4}",42, 0.025, 1);
		drawLatex(0.48,0.86,"x~2.4#times10^{-4}",42, 0.025, 1);
		drawLatex(0.67,0.86,"x~1.1#times10^{-4}",42, 0.025, 1);
		drawLatex(0.87,0.86,"x~6.0#times10^{-5}",42, 0.025, 1);

		drawLatex(0.15, 0.80, "Pb+Pb #rightarrow Pb'+Pb'+J/#psi",  42,        0.05,      1 );
		drawLatex(0.56,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
		drawLatex(0.15,0.94,"#bf{CMS} #it{Preliminary}",42, 0.05, 1);

		return hist;
	}
};

// void FigureFrameStrategy()
// {
//     auto c = new TCanvas();
//     Sigma_Log_Strategy s;
//     s.Apply();
// }