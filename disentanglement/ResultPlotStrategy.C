
#include "../common/PlotStrategy.h"

enum class ResultPlotStrategyList
{
	Empty,

	DSigmaDy_CMS,
	DSigmaDy_00_0X_XX_CMS,

	DSigmaDy_ALICE_2019,
	DSigmaDy_ALICE_2021,
	DSigmaDy_LHCb_2022,

	DSigmaDy_LTA_SS,
	DSigmaDy_LTA_WS,

	Sigma_CMS,
	Sigma_ALICE_2019,
	Sigma_ALICE_2021,
	Sigma_LHCb_2022,

	R_CMS,
	R_ALICE_2019,
	R_ALICE_2021,
	R_LHCb_2022
};

struct Empty_PlotStrategy : PlotStrategy
{
	Empty_PlotStrategy(const AnalysisData& data_) :     PlotStrategy{data_} {}

	void Apply(TLegend*	leg)	override
	{
		ge 		= new TGraphErrors		();
		leg->AddEntry(ge,	"",	"");
	}
};

struct Sigma_PlotStrategy : PlotStrategy
{
	std::vector<double> X_AXIS_ERR		= std::vector<double>(data.GetSize("W"), 3.8);

	Sigma_PlotStrategy(const AnalysisData& data_) :     PlotStrategy{data_} {}

	void Apply(TLegend*	leg)	override
	{
		std::vector<double> SysErrLow, SysErrHigh;

		if(data.IsMapKeyExist("Sigma_SysErrLow") && data.IsMapKeyExist("Sigma_SysErrHigh"))	//If the data has AsymmErr
		{
			SysErrLow  = data.Get("Sigma_SysErrLow");
			SysErrHigh = data.Get("Sigma_SysErrHigh");
		}
		else
		{
			SysErrLow  = data.Get("Sigma_SysErr");
			SysErrHigh = data.Get("Sigma_SysErr");
		}

		ge 		= new TGraphErrors		(data.GetSize("Sigma"),	data.Get("W").data(),	data.Get("Sigma").data(),	nullptr,	data.Get("Sigma_Err").data()	);
		gae 	= new TGraphAsymmErrors (data.GetSize("Sigma"),
										data.Get("W").data(),        	data.Get("Sigma").data(),
										X_AXIS_ERR.data(),  			X_AXIS_ERR.data(),
										SysErrLow.data(), 				SysErrHigh.data());
	}
};

struct R_PlotStrategy : PlotStrategy
{
	std::vector<double> X_AXIS_ERR		= std::vector<double>(data.GetSize("X"), 4.5e-4);

	R_PlotStrategy(const AnalysisData& data_) :     PlotStrategy{data_} {}

	void Apply(TLegend*	leg)	override
	{
		std::vector<double> SysErrLow, SysErrHigh;

		if(data.IsMapKeyExist("R_SysErrLow") && data.IsMapKeyExist("R_SysErrHigh"))	//If the data has AsymmErr
		{
			SysErrLow  = data.Get("R_SysErrLow");
			SysErrHigh = data.Get("R_SysErrHigh");
		}
		else
		{
			SysErrLow  = data.Get("R_SysErr");
			SysErrHigh = data.Get("R_SysErr");
		}

		ge 		= new TGraphErrors		(data.GetSize("R"),	data.Get("X").data(),	data.Get("R").data(),	nullptr,	data.Get("R_Err").data()	);
		gae 	= new TGraphAsymmErrors	(data.GetSize("R"),
										data.Get("X").data(),        	data.Get("R").data(),
										X_AXIS_ERR.data(),  			X_AXIS_ERR.data(),
										SysErrLow.data(), 				SysErrHigh.data());
	}
};

struct DSigmaDy_PlotStrategy : PlotStrategy
{
	TString NeuConfig;

	DSigmaDy_PlotStrategy(const AnalysisData& data_,	TString NeuConfig_ = "AnAn") :     PlotStrategy{data_},	NeuConfig{NeuConfig_}	{}

	void Apply(TLegend*	leg)	override
	{
		auto temp = data.Get("Dy");

		if(temp[0]>0)
		{
			for (int i = 0; i < temp.size(); ++i)
			{
				temp[i] *= -1;
			}	
		}

		std::vector<double> SysErrLow, SysErrHigh;
		TString TempName = "DSigmaDy_"+NeuConfig;

		if(data.IsMapKeyExist(TempName+"_SysErrLow") && data.IsMapKeyExist(TempName+"_SysErrHigh"))	//If the data has AsymmErr
		{
			SysErrLow  = data.Get(TempName+"_SysErrLow");
			SysErrHigh = data.Get(TempName+"_SysErrHigh");
		}
		else
		{
			SysErrLow  = data.Get(TempName+"_SysErr");
			SysErrHigh = data.Get(TempName+"_SysErr");
		}

		ge  	= new TGraphErrors		(temp.size(),				temp.data(),					data.Get(TempName).data(),
																	data.Get("Dy_Err").data(),		data.Get(TempName+"_Err").data());
		gae 	= new TGraphAsymmErrors	(data.GetSize("Dy"),		temp.data(),					data.Get(TempName).data(),
																	data.Get("Dy_Err").data(),		data.Get("Dy_Err").data(),
																	SysErrLow.data(), 				SysErrHigh.data());
	}
};

struct Sigma_CMS_PlotStrategy : Sigma_PlotStrategy
{
	Sigma_CMS_PlotStrategy(const AnalysisData& data_) :     Sigma_PlotStrategy{data_} {}

	void Apply(TLegend*	leg)	override
	{
		Sigma_PlotStrategy::Apply(leg);
		leg->AddEntry(ge,	"CMS",	"p");

		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");
		ge  ->SetMarkerStyle(20);
		ge  ->SetMarkerSize(1.5);
		ge  ->SetMarkerColor(2);
		ge  ->SetLineColor(2);
		ge  ->SetLineWidth(2);
		ge  ->Draw("pezsame");
	}
};

struct Sigma_ALICE_2019_Strategy : Sigma_PlotStrategy
{
	Sigma_ALICE_2019_Strategy(const AnalysisData& data_)	: Sigma_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		Sigma_PlotStrategy::Apply(leg);
		leg->AddEntry(ge,	"ALICE* (-4 < y < -3.5)",	"p");

		gae	->SetMarkerStyle(24);
		gae	->SetFillColorAlpha(16, 0.7);
		gae	->SetFillStyle(1001);
		gae	->Draw("2same");
		ge	->SetMarkerStyle(25);
		ge	->SetMarkerColor(4);
		ge	->SetMarkerSize(1.5);
		ge	->SetLineColor(4);
		ge	->SetLineWidth(2);
		ge	->Draw("pezsame");
	}
};

struct Sigma_ALICE_2021_Strategy : Sigma_PlotStrategy
{
	Sigma_ALICE_2021_Strategy(const AnalysisData& data_)	: Sigma_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		Sigma_PlotStrategy::Apply(leg);
		leg->AddEntry(ge,	"ALICE* (|y| < 0.15)",	"p");

		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");
		ge	->SetMarkerStyle(24);
		ge	->SetMarkerColor(4);
		ge	->SetMarkerSize(1.5);
		ge	->SetLineColor(4);
		ge	->SetLineWidth(2);
		ge	->Draw("pezsame");
	}
};

struct Sigma_LHCb_2022_Strategy : Sigma_PlotStrategy
{
	Sigma_LHCb_2022_Strategy(const AnalysisData& data_)	: Sigma_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		Sigma_PlotStrategy::Apply(leg);
		leg->AddEntry(ge,	"LHCb* (-4.5 < y < -3.5)",	"p");

		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");
		ge->SetMarkerStyle(26);
		ge->SetMarkerColor(4);
		ge->SetMarkerSize(1.5);
		ge->SetLineColor(4);
		ge->SetLineWidth(2);
		ge	->Draw("pezsame");
	}
};

struct DSigmaDy_CMS_Strategy : DSigmaDy_PlotStrategy
{
	DSigmaDy_CMS_Strategy(const AnalysisData& data_,	TString NeuConfig_ = "AnAn",	int index_ = 0)	: DSigmaDy_PlotStrategy{data_, NeuConfig_} {	index = index_;	}
	
	void Apply(TLegend*	leg)	override
	{
		DSigmaDy_PlotStrategy::Apply(leg);
		leg->AddEntry(gae,	"CMS",	"lpf");

		gae ->SetMarkerStyle(20);
		gae ->SetMarkerColor(2);
		gae ->SetLineColor(2);
		gae ->SetFillColorAlpha(16, 0.5);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");

		ge  ->SetMarkerStyle(20);
		ge  ->SetMarkerSize(1.5);
		ge  ->SetMarkerColor(2);
		ge  ->SetLineColor(2);
		ge  ->SetLineWidth(2);
		ge  ->Draw("pezsame");
	}
};


struct DSigmaDy_00_0X_XX_CMS_Strategy : DSigmaDy_PlotStrategy
{
	TString legName; int color; int mstyle;
	DSigmaDy_00_0X_XX_CMS_Strategy(const AnalysisData& data_,	TString NeuConfig_ = "AnAn", int index_ = 0, TString legName="", int color=1, int mstyle=1 )	: DSigmaDy_PlotStrategy{data_, NeuConfig_},legName{legName},color{color},mstyle{mstyle} {	index = index_;	}
	
	void Apply(TLegend*	leg)	override
	{
		DSigmaDy_PlotStrategy::Apply(leg);
		leg->AddEntry(gae,	legName,	"pf");

		gae ->SetMarkerStyle(mstyle);
		gae ->SetMarkerColor(color);
		gae ->SetLineColor(color);
		gae ->SetFillColorAlpha(16, 0.5);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");

		ge  ->SetMarkerStyle(mstyle);
		ge  ->SetMarkerSize(1.5);
		ge  ->SetMarkerColor(color);
		ge  ->SetLineColor(color);
		ge  ->SetLineWidth(2);
		ge  ->Draw("pezsame");
	}
};



struct DSigmaDy_ALICE_2019_Strategy : DSigmaDy_PlotStrategy
{
	DSigmaDy_ALICE_2019_Strategy(const AnalysisData& data_,	TString NeuConfig_ = "AnAn",	int index_ = 0)	: DSigmaDy_PlotStrategy{data_, NeuConfig_} {	index = index_;	}
	
	void Apply(TLegend*	leg)	override
	{
		DSigmaDy_PlotStrategy::Apply(leg);
		leg->AddEntry(gae,	"ALICE 2019",	"lpf");

		gae ->SetMarkerStyle(25);
		gae ->SetMarkerColor(4);
		gae ->SetLineColor(4);
		gae ->SetFillColor(kGray);
		gae ->SetFillStyle(3001);
		gae ->Draw("2same");

		ge  ->SetMarkerStyle(25);
		ge  ->SetMarkerSize(1.5);
		ge  ->SetMarkerColor(4);
		ge  ->SetLineColor(4);
		ge  ->Draw("pezsame");
	}
};

struct	DSigmaDy_ALICE_2021_Strategy : DSigmaDy_PlotStrategy
{
	DSigmaDy_ALICE_2021_Strategy(const AnalysisData& data_,	TString NeuConfig_ = "AnAn",	int index_ = 0)	: DSigmaDy_PlotStrategy{data_, NeuConfig_} {	index = index_;	}
	
	void Apply(TLegend*	leg)	override
	{
		DSigmaDy_PlotStrategy::Apply(leg);
		leg->AddEntry(gae,	"ALICE 2021",	"lpf");

		gae ->SetMarkerStyle(24);
		gae ->SetMarkerColor(4);
		gae ->SetLineColor(4);
		gae ->SetFillColor(kGray);
		gae ->SetFillStyle(3001);
		gae ->Draw("2same");

		ge  ->SetMarkerStyle(24);
		ge  ->SetMarkerSize(1.5);
		ge  ->SetMarkerColor(4);
		ge  ->SetLineColor(4);
		ge  ->Draw("pezsame");
	}
};

struct DSigmaDy_LHCb_2022_Strategy : DSigmaDy_PlotStrategy
{
	DSigmaDy_LHCb_2022_Strategy(const AnalysisData& data_,	TString NeuConfig_ = "AnAn",	int index_ = 0)	: DSigmaDy_PlotStrategy{data_, NeuConfig_} {	index = index_;	}
	
	void Apply(TLegend*	leg)	override
	{
		DSigmaDy_PlotStrategy::Apply(leg);
		leg->AddEntry(gae,	"LHCb  2022",	"lpf");


		gae ->SetMarkerStyle(26);
		gae ->SetFillColor(kGray);
		gae ->SetFillStyle(3001);
		gae ->Draw("2same");

		ge  ->SetMarkerStyle(26);
		ge  ->SetMarkerSize(1.5);
		ge  ->Draw("pezsame");
	}
};

struct R_CMS_Strategy : R_PlotStrategy
{
	R_CMS_Strategy(const AnalysisData& data_)	: R_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		X_AXIS_ERR = {8e-6,	7e-6,	5e-6,	4.5e-4,	4e-4,	3e-4};
		R_PlotStrategy::Apply(leg);
		leg->AddEntry(ge,	"CMS",	"p");
		
		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");
		ge->SetMarkerStyle(20);
		ge->SetMarkerColor(2);
		ge->SetMarkerSize(1.5);
		ge->SetLineColor(2);
		ge->SetLineWidth(2);
		ge->Draw("pezsame");
	}
};

struct R_ALICE_2019_Strategy : R_PlotStrategy
{
	R_ALICE_2019_Strategy(const AnalysisData& data_)	: R_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		X_AXIS_ERR = {2.5e-3,2e-3};
		R_PlotStrategy::Apply(leg);
		leg->AddEntry(ge,	"ALICE* (-4 < y < -3.5)",	"p");

		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");
		ge->SetMarkerStyle(25);
		ge->SetMarkerColor(4);
		ge->SetMarkerSize(1.5);
		ge->SetLineColor(4);
		ge->SetLineWidth(2);
		ge->Draw("pezsame");
	}
};

struct R_ALICE_2021_Strategy : R_PlotStrategy
{
	R_ALICE_2021_Strategy(const AnalysisData& data_)	: R_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		X_AXIS_ERR = {5e-5};
		R_PlotStrategy::Apply(leg);
		leg->AddEntry(ge,	"ALICE* (|y| < 0.15)",	"p");

		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");
		ge->SetMarkerStyle(24);
		ge->SetMarkerColor(4);
		ge->SetMarkerSize(1.5);
		ge->SetLineColor(4);
		ge->SetLineWidth(2);
		ge->Draw("pezsame");
	}
};

struct R_LHCb_2022_Strategy : R_PlotStrategy
{
	R_LHCb_2022_Strategy(const AnalysisData& data_)	: R_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		X_AXIS_ERR = {2.5e-3,2e-3};
		R_PlotStrategy::Apply(leg);
		leg->AddEntry(ge,	"LHCb* (-4.5 < y < -3.5)",	"p");

		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");
		ge->SetMarkerStyle(26);
		ge->SetMarkerColor(4);
		ge->SetMarkerSize(1.5);
		ge->SetLineColor(4);
		ge->SetLineWidth(2);
		ge->Draw("pezsame");
	}
};
