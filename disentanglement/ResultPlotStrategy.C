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
	std::unique_ptr< TGraphErrors > 		ge;
	std::unique_ptr< TGraphAsymmErrors >	gae;

	ResultPlotStrategy(const AnalysisData& data_) : data{data_} {};

	virtual void Apply() = 0;
	virtual ~ResultPlotStrategy() = default;
};

struct Sigma_PlotStrategy : ResultPlotStrategy
{
	std::vector<double> X_AXIS_ERR		= std::vector<double>(data.GetSize("W"), 3.8);

	Sigma_PlotStrategy(const AnalysisData& data_) :     ResultPlotStrategy{data_} {}

	void Apply()	override
	{
		ge 		= std::make_unique< TGraphErrors >		(data.GetSize("W"),	data.Get("W").data(),	data.Get("Sigma").data(),	nullptr,	data.Get("Sigma_Err").data()	);
		// gae 	= std::make_unique< TGraphAsymmErrors >		(data.GetSize("W"),
		// 										data.Get("W").data(),        data.Get("Sigma").data(),
		// 										X_AXIS_ERR.data(),  X_AXIS_ERR.data(),
		// 										Sigmas_TotalSysErr.data(), Sigmas_TotalSysErr.data());
	}
};

struct DSigmaDy_PlotStrategy : ResultPlotStrategy
{
	std::vector<double> X_AXIS_ERR		= std::vector<double>(data.GetSize("W"), 3.8);

	DSigmaDy_PlotStrategy(const AnalysisData& data_) :     ResultPlotStrategy{data_} {}

	void Apply()	override
	{
		auto temp = data.Get("Dy");

		if(temp[0]>0)
		{
			for (int i = 0; i < temp.size(); ++i)
			{
				temp[i] *= -1;
			}	
		}

		ge  	= std::make_unique< TGraphErrors >		(temp.size(),		temp.data(),	data.Get("DSigmaDy_AnAn").data(),
														data.Get("Dy_Err").data(),	data.Get("DSigmaDy_AnAn_Err").data());
		// gae 	= std::make_unique< TGraphAsymmErrors >	(data.GetSize("Dy"),		data.Get("Dy").data(),	data.Get("DSigmaDy_AnAn").data(),
		// 												data.Get("Dy_Err").data(),	data.Get("Dy_Err").data(),
		// 												data.Get("DSigmaDy_AnAn_SysErr").data(),	data.Get("DSigmaDy_AnAn_SysErr").data());
	}
};

struct Sigma_CMS_PlotStrategy : Sigma_PlotStrategy
{
	Sigma_CMS_PlotStrategy(const AnalysisData& data_) :     Sigma_PlotStrategy{data_} {}

	void Apply()	override
	{
		Sigma_PlotStrategy::Apply();

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
	
	void Apply()	override
	{
		Sigma_PlotStrategy::Apply();

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
	
	void Apply()	override
	{
		Sigma_PlotStrategy::Apply();

		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");
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
	
	void Apply()	override
	{
		Sigma_PlotStrategy::Apply();

		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		ge->SetMarkerStyle(26);
		ge->SetMarkerColor(4);
		ge->SetMarkerSize(1.5);
		ge->SetLineColor(4);
		ge->SetLineWidth(2);
	}
};

struct DSigmaDy_CMS_Strategy : DSigmaDy_PlotStrategy
{
	DSigmaDy_CMS_Strategy(const AnalysisData& data_)	: DSigmaDy_PlotStrategy{data_} {}
	
	void Apply()	override
	{
		DSigmaDy_PlotStrategy::Apply();
		
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
	
	void Apply()	override
	{
		DSigmaDy_PlotStrategy::Apply();

		// gae  ->SetMarkerStyle(21);
		// gae ->SetFillColor(8);
		// gae ->SetFillStyle(3001);
		// gae ->Draw("2same");
		ge  ->SetMarkerStyle(21);
		ge  ->SetMarkerSize(1.);
		ge  ->Draw("pezsame");
	}
};