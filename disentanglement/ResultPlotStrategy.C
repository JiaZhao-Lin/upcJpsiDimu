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

struct Sigma_CMS_PlotStrategy : ResultPlotStrategy
{
	std::vector<double> X_AXIS_ERR		= std::vector<double>(data.GetSize("W"), 3.8);

	Sigma_CMS_PlotStrategy(const AnalysisData& data_) :     ResultPlotStrategy{data_} {}

	void Apply()
	{
		ge 		= std::make_unique< TGraphErrors >		(data.GetSize("W"),	data.Get("W").data(),	data.Get("Sigma").data(),	nullptr,	data.Get("Sigma_Err").data()	);
		// gae 	= std::make_unique< TGraphAsymmErrors >		(data.GetSize("W"),
		// 										data.Get("W").data(),        data.Get("Sigma").data(),
		// 										X_AXIS_ERR.data(),  X_AXIS_ERR.data(),
		// 										Sigmas_TotalSysErr.data(), Sigmas_TotalSysErr.data());

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

struct Sigma_ALICE_2019_Strategy : ResultPlotStrategy
{
	std::vector<double> X_AXIS_ERR		= std::vector<double>(data.GetSize("W"), 3.8);

	Sigma_ALICE_2019_Strategy(const AnalysisData& data_)	: ResultPlotStrategy{data_} {}
	
	void Apply()
	{
		ge		= std::make_unique< TGraphErrors >		(data.GetSize("W"),	data.Get("W").data(),	data.Get("Sigma").data(),	nullptr,	data.Get("Sigma_Err").data()	);
		// gae 	= std::make_unique< TGraphAsymmErrors >	(data.GetSize("W"),	data.Get("W").data(),	data.Get("Sigma").data(),
		// 										X_AXIS_ERR.data(),	X_AXIS_ERR.data(),
		// 										data_SysErrLow.data(),	data_SysErrHig.data());

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

struct Sigma_ALICE_2021_Strategy : ResultPlotStrategy
{
	std::vector<double> X_AXIS_ERR		= std::vector<double>(data.GetSize("W"), 3.8);

	Sigma_ALICE_2021_Strategy(const AnalysisData& data_)	: ResultPlotStrategy{data_} {}
	
	void Apply()
	{
		ge		= std::make_unique< TGraphErrors >		(data.GetSize("W"),	data.Get("W").data(),	data.Get("Sigma").data(),	nullptr,	data.Get("Sigma_Err").data()	);
		// gae 	= std::make_unique< TGraphAsymmErrors >	(data.GetSize("W"),	data.Get("W").data(),	data.Get("Sigma").data(),
		// 										X_AXIS_ERR.data(),	X_AXIS_ERR.data(),
		// 										data_SysErrLow.data(),	data_SysErrHig.data());

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

struct Sigma_LHCb_2022_Strategy : ResultPlotStrategy
{
	std::vector<double> X_AXIS_ERR		= std::vector<double>(data.GetSize("W"), 3.8);

	Sigma_LHCb_2022_Strategy(const AnalysisData& data_)	: ResultPlotStrategy{data_} {}
	
	void Apply()
	{
		ge		= std::make_unique< TGraphErrors >		(data.GetSize("W"),	data.Get("W").data(),	data.Get("Sigma").data(),	nullptr,	data.Get("Sigma_Err").data()	);
		// gae 	= std::make_unique< TGraphAsymmErrors >	(data.GetSize("W"),	data.Get("W").data(),	data.Get("Sigma").data(),
		// 										X_AXIS_ERR.data(),	X_AXIS_ERR.data(),
		// 										data_SysErrLow.data(),	data_SysErrHig.data());

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