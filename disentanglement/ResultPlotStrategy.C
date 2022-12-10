#include "AnalysisData.C"

enum class PlotStrategyList
{
	DSigmaDy_CMS,
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

struct ResultPlotStrategy
{
	const AnalysisData& 					data;
	TGraphErrors* 		ge;
	TGraphAsymmErrors*	gae;

	ResultPlotStrategy(const AnalysisData& data_) : data{data_} {};

	virtual void Apply(TLegend*	leg) = 0;
	virtual ~ResultPlotStrategy() = default;
};

struct Sigma_PlotStrategy : ResultPlotStrategy
{
	std::vector<double> X_AXIS_ERR		= std::vector<double>(data.GetSize("W"), 3.8);

	Sigma_PlotStrategy(const AnalysisData& data_) :     ResultPlotStrategy{data_} {}

	void Apply(TLegend*	leg)	override
	{
		ge 		= new TGraphErrors		(data.GetSize("Sigma"),	data.Get("W").data(),	data.Get("Sigma").data(),	nullptr,	data.Get("Sigma_Err").data()	);
		// gae 	= std::make_unique< TGraphAsymmErrors >		(data.GetSize("Sigma"),
															// data.Get("W").data(),        	data.Get("Sigma").data(),
															// X_AXIS_ERR.data(),  			X_AXIS_ERR.data(),
															// data.Get("Sigma_SysErr").data(), data.Get("Sigma_SysErr").data());
	}
};

struct R_PlotStrategy : ResultPlotStrategy
{
	std::vector<double> X_AXIS_ERR		= std::vector<double>(data.GetSize("X"), 4.5e-4);

	R_PlotStrategy(const AnalysisData& data_) :     ResultPlotStrategy{data_} {}

	void Apply(TLegend*	leg)	override
	{
		ge 		= new TGraphErrors		(data.GetSize("R"),	data.Get("X").data(),	data.Get("R").data(),	nullptr,	data.Get("R_Err").data()	);
		// gae 	= std::make_unique< TGraphAsymmErrors >		(data.GetSize("R"),
															// data.Get("X").data(),        	data.Get("R").data(),
															// X_AXIS_ERR.data(),  			X_AXIS_ERR.data(),
															// data.Get("R_SysErr").data(), data.Get("R_SysErr").data());
	}
};

struct DSigmaDy_PlotStrategy : ResultPlotStrategy
{
	DSigmaDy_PlotStrategy(const AnalysisData& data_) :     ResultPlotStrategy{data_} {}

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

		ge  	= new TGraphErrors		(temp.size(),		temp.data(),	data.Get("DSigmaDy_AnAn").data(),
														data.Get("Dy_Err").data(),	data.Get("DSigmaDy_AnAn_Err").data());
		// gae 	= std::make_unique< TGraphAsymmErrors >	(data.GetSize("Dy"),		data.Get("Dy").data(),	data.Get("DSigmaDy_AnAn").data(),
														// data.Get("Dy_Err").data(),	data.Get("Dy_Err").data(),
														// data.Get("DSigmaDy_AnAn_SysErr").data(),	data.Get("DSigmaDy_AnAn_SysErr").data());
	}
};

struct Sigma_CMS_PlotStrategy : Sigma_PlotStrategy
{
	Sigma_CMS_PlotStrategy(const AnalysisData& data_) :     Sigma_PlotStrategy{data_} {}

	void Apply(TLegend*	leg)	override
	{
		Sigma_PlotStrategy::Apply(leg);
		leg->AddEntry((TObject *)ge,	"CMS",	"p");

		// gae ->SetMarkerStyle(24);
		// gae ->SetFillColorAlpha(16, 0.7);
		// gae ->SetFillStyle(1001);
		// gae ->Draw("2same");
		ge  ->SetMarkerStyle(20);
		ge  ->SetMarkerSize(1.6);
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
		leg->AddEntry((TObject *)ge,	"ALICE* (-4 < y < -3.5)",	"p");

		// gae	->SetMarkerStyle(24);
		// gae	->SetFillColorAlpha(16, 0.7);
		// gae	->SetFillStyle(1001);
		// gae	->Draw("2same");
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
		leg->AddEntry((TObject *)ge,	"ALICE* (|y| < 0.15)",	"p");

		// gae ->SetMarkerStyle(24);
		// gae ->SetFillColorAlpha(16, 0.7);
		// gae ->SetFillStyle(1001);
		// gae ->Draw("2same");
		ge	->SetMarkerStyle(24);
		ge	->SetMarkerColor(4);
		ge	->SetMarkerSize(1.6);
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
		leg->AddEntry((TObject *)ge,	"LHCb* (-4.5 < y < -3.5)",	"p");

		// gae ->SetMarkerStyle(24);
		// gae ->SetFillColorAlpha(16, 0.7);
		// gae ->SetFillStyle(1001);
		// gae ->Draw("2same");
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
	DSigmaDy_CMS_Strategy(const AnalysisData& data_)	: DSigmaDy_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		DSigmaDy_PlotStrategy::Apply(leg);
		leg->AddEntry((TObject *)ge,	"CMS",	"lpf");

		// gae ->SetMarkerStyle(20);
		// gae ->SetMarkerColor(2);
		// gae ->SetLineColor(2);
		// gae ->SetFillColorAlpha(16, 0.5);
		// gae ->SetFillStyle(1001);
		// gae ->Draw("2same");
		ge  ->SetMarkerStyle(20);
		ge  ->SetMarkerSize(1.);
		ge  ->SetMarkerColor(2);
		ge  ->SetLineColor(2);
		ge  ->SetLineWidth(2);
		ge  ->Draw("pezsame");
	}
};

struct DSigmaDy_ALICE_2019_Strategy : DSigmaDy_PlotStrategy
{
	DSigmaDy_ALICE_2019_Strategy(const AnalysisData& data_)	: DSigmaDy_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		DSigmaDy_PlotStrategy::Apply(leg);
		leg->AddEntry((TObject *)ge,	"ALICE 2019",	"lpf");

		// gae  ->SetMarkerStyle(21);
		// gae ->SetFillColor(8);
		// gae ->SetFillStyle(3001);
		// gae ->Draw("2same");
		ge  ->SetMarkerStyle(21);
		ge  ->SetMarkerSize(1.);
		ge  ->Draw("pezsame");
	}
};

struct	DSigmaDy_ALICE_2021_Strategy : DSigmaDy_PlotStrategy
{
	DSigmaDy_ALICE_2021_Strategy(const AnalysisData& data_)	: DSigmaDy_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		DSigmaDy_PlotStrategy::Apply(leg);
		leg->AddEntry((TObject *)ge,	"ALICE 2021",	"lpf");

		// gae  ->SetMarkerStyle(21);
		// gae ->SetFillColor(kGreen);
		// gae ->SetFillStyle(3001);
		// gae ->Draw("2same");
		ge  ->SetMarkerStyle(21);
		ge  ->SetMarkerSize(1.);
		ge  ->Draw("pezsame");
	}
};

struct DSigmaDy_LHCb_2022_Strategy : DSigmaDy_PlotStrategy
{
	DSigmaDy_LHCb_2022_Strategy(const AnalysisData& data_)	: DSigmaDy_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		DSigmaDy_PlotStrategy::Apply(leg);
		leg->AddEntry((TObject *)ge,	"LHCb 2022",	"lpf");

		// gae  ->SetMarkerStyle(24);
		// gae ->SetFillColor(kGray);
		// gae ->SetFillStyle(3001);
		// gae ->Draw("2same");
		ge  ->SetMarkerStyle(26);
		ge  ->SetMarkerSize(1.);
		ge  ->Draw("pezsame");
	}
};

struct R_CMS_Strategy : R_PlotStrategy
{
	R_CMS_Strategy(const AnalysisData& data_)	: R_PlotStrategy{data_} {}
	
	void Apply(TLegend*	leg)	override
	{
		std::vector<double> X_AXIS_ERR = {8e-6,3e-4,7e-6,4e-4,5e-6,4.5e-4};
		R_PlotStrategy::Apply(leg);
		leg->AddEntry((TObject *)ge,	"CMS",	"p");
		
		// gae ->SetMarkerStyle(24);
		// gae ->SetFillColorAlpha(16, 0.7);
		// gae ->SetFillStyle(1001);
		// gae ->Draw("2same");
		ge->SetMarkerStyle(20);
		ge->SetMarkerColor(2);
		ge->SetMarkerSize(1.6);
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
		leg->AddEntry((TObject *)ge,	"ALICE* (-4 < y < -3.5)",	"p");

		// gae ->SetMarkerStyle(24);
		// gae ->SetFillColorAlpha(16, 0.7);
		// gae ->SetFillStyle(1001);
		// gae ->Draw("2same");
		ge->SetMarkerStyle(25);
		ge->SetMarkerColor(4);
		ge->SetMarkerSize(1.6);
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
		leg->AddEntry((TObject *)ge,	"ALICE* (|y| < 0.15)",	"p");

		// gae ->SetMarkerStyle(24);
		// gae ->SetFillColorAlpha(16, 0.7);
		// gae ->SetFillStyle(1001);
		// gae ->Draw("2same");
		ge->SetMarkerStyle(24);
		ge->SetMarkerColor(4);
		ge->SetMarkerSize(1.6);
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
		leg->AddEntry((TObject *)ge,	"LHCb* (-4.5 < y < -3.5)",	"p");

		// gae ->SetMarkerStyle(24);
		// gae ->SetFillColorAlpha(16, 0.7);
		// gae ->SetFillStyle(1001);
		// gae ->Draw("2same");
		ge->SetMarkerStyle(26);
		ge->SetMarkerColor(4);
		ge->SetMarkerSize(1.6);
		ge->SetLineColor(4);
		ge->SetLineWidth(2);
		ge->Draw("pezsame");
	}
};