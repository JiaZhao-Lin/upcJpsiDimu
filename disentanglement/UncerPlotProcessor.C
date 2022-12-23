#include "ResultFrameStrategy.C"
#include "UncerAnalyzer.C"
#include "../common/PlotStrategy.h"

//shift the x axis by a certain amount when plotting the AnalysisData. Used to compare the results
struct ShiftPlot : PlotStrategy
{
	static void Plot(const AnalysisData& data_, const TString paramX, const TString paramY, const double shift,
				const int style, const int color, TLegend* leg, const TString legendName)
	{
		auto x 		= data_.Get(paramX);
		auto y		= data_.Get(paramY);
		auto y_Err	= data_.Get(paramY+"_Err");

		// the shift is used to shift the x-axis. it is different for log scale and linear scale
		if(paramX == "X")
		{
			for (int i = 0; i < x.size(); ++i)
			{	
				if 	(x[i] > 1e-3)	{ x[i] += shift; 		}
				else 				{ x[i] += shift/50; 	}
			}
		}
		else
		{
			for (int i = 0; i < x.size(); ++i){	x[i] += shift; }
		}

		TGraphErrors* ge	= new TGraphErrors(x.size(),	x.data(),	y.data(),	0,	y_Err.data()	);
		ge->SetMarkerStyle(style);
		ge->SetMarkerSize(1.1);
		ge->SetMarkerColor(color);
		ge->SetLineColor(color);
		ge->SetLineWidth(2);
		ge->Draw("pezsame");

		leg->AddEntry(ge,	legendName,	"p");
	}
};


struct UncerPlotProcessor
{
	std::unique_ptr<FrameStrategy>		frame;
	struct UncerAnalyzer				uncerAnalyzer;

	UncerPlotProcessor() {uncerAnalyzer.Handle();}
	void SaveAs(TString fileName){	frame->SaveAs(fileName);	};

	//Take all the AnalysisData from UncerAnalyzer and plot them using ShiftPlot function using the Sigma_Log_Strategy
	void CompareSigma()
	{
		frame = std::make_unique<Sigma_Log_Strategy>();
		frame->Apply();
		delete frame->legends[0];
		frame->legends[0] = new TLegend(0.60, 0.15, 0.90, 0.45);
		frame->legends[0]->SetTextSize(0.020);

		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_AnaData,				"W",	"Sigma",	0,	47,	1,	frame->legends[0],	"Default"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly4_AnaData,				"W",	"Sigma",	0,	24,	1,	frame->legends[0],	"CB_Poly4"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_FixCBAN_Poly3_AnaData,		"W",	"Sigma",	0,	24,	2,	frame->legends[0],	"CB_FixCBAN_Poly3"	);
		ShiftPlot::Plot(uncerAnalyzer.CBG_Poly3_AnaData,				"W",	"Sigma",	0,	24,	4,	frame->legends[0],	"CBG_Poly3"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_SdB_AnaData,			"W",	"Sigma",	-5,	24,	6,	frame->legends[0],	"QEDPtShape"		);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_NarrMass_AnaData,		"W",	"Sigma",	-5,	24,	8,	frame->legends[0],	"NarrMass"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_WideMass_AnaData,		"W",	"Sigma",	-5,	24,	9,	frame->legends[0],	"WideMass"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_temp_AnaData,			"W",	"Sigma",	-5,	24,	28,	frame->legends[0],	"PtFitting"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_CohPtCut_AnaData,		"W",	"Sigma",	5,	24,	30,	frame->legends[0],	"CohPtCut"			);

		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_looseHF_AnaData,		"W",	"Sigma",	5,	26,	1,	frame->legends[0],	"HFveto"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_fluxPM_AnaData,		"W",	"Sigma",	5,	26,	2,	frame->legends[0],	"PhotonFluxPM"		);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_fluxM_AnaData,			"W",	"Sigma",	5,	26,	4,	frame->legends[0],	"PhotonFluxM" 		);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_PU_AnaData,			"W",	"Sigma",	5,	26,	6,	frame->legends[0],	"n-PileUp"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_TnP_Low_AnaData,		"W",	"Sigma",	5,	26,	8,	frame->legends[0],	"TnP"				);
		frame->DrawLegend();

		SaveAs("SigmaCompared.pdf");
	}

	//Take all the AnalysisData from UncerAnalyzer and plot them using ShiftPlot function using the R_Strategy
	void CompareR()
	{
		frame = std::make_unique<R_Strategy>();
		frame->Apply();
		delete frame->legends[0];
		frame->legends[0] = new TLegend(0.60, 0.15, 0.90, 0.45);
		frame->legends[0]->SetTextSize(0.020);

		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_AnaData,				"X",	"R",	0,		47,	1,	frame->legends[0],	"Default"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly4_AnaData,				"X",	"R",	0,		24,	1,	frame->legends[0],	"CB_Poly4"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_FixCBAN_Poly3_AnaData,		"X",	"R",	0,		24,	2,	frame->legends[0],	"CB_FixCBAN_Poly3"	);
		ShiftPlot::Plot(uncerAnalyzer.CBG_Poly3_AnaData,				"X",	"R",	0,		24,	4,	frame->legends[0],	"CBG_Poly3"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_SdB_AnaData,			"X",	"R",	-4e-4,	24,	6,	frame->legends[0],	"QEDPtShape"		);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_NarrMass_AnaData,		"X",	"R",	-4e-4,	24,	8,	frame->legends[0],	"NarrMass"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_WideMass_AnaData,		"X",	"R",	-4e-4,	24,	9,	frame->legends[0],	"WideMass"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_temp_AnaData,			"X",	"R",	-4e-4,	24,	28,	frame->legends[0],	"PtFitting"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_CohPtCut_AnaData,		"X",	"R",	4e-4,	24,	30,	frame->legends[0],	"CohPtCut"			);

		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_looseHF_AnaData,		"X",	"R",	4e-4,	26,	1,	frame->legends[0],	"HFveto"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_fluxPM_AnaData,		"X",	"R",	4e-4,	26,	2,	frame->legends[0],	"PhotonFluxPM"		);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_fluxM_AnaData,			"X",	"R",	4e-4,	26,	4,	frame->legends[0],	"PhotonFluxM" 		);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_PU_AnaData,			"X",	"R",	4e-4,	26,	6,	frame->legends[0],	"n-PileUp"			);
		ShiftPlot::Plot(uncerAnalyzer.CB_Poly3_TnP_Low_AnaData,		"X",	"R",	4e-4,	26,	8,	frame->legends[0],	"TnP"				);
		frame->DrawLegend();

		SaveAs("RCompared.pdf");
	}
};

void UncerPlotProcessor()
{
	struct UncerPlotProcessor processor;
	
	processor.CompareSigma();
	processor.CompareR();
}