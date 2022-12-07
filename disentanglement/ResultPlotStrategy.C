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
	const AnalysisData& data;
	ResultPlotStrategy(const AnalysisData& data_) : data{data_} {};

	virtual void Apply() = 0;
	virtual ~ResultPlotStrategy() = default;
};

struct Sigma_CMS_PlotStrategy : ResultPlotStrategy
{
	TGraphErrors* ge
	TGraphAsymmErrors*gae 
	Sigma_CMS_PlotStrategy(const AnalysisData& data_) :     ResultPlotStrategy{data_} {}

	void Apply()
	{
		ge 		= new TGraphErrors(ShadowRatio_ParamsMap.at("Ws").size(),	&ShadowRatio_ParamsMap.at("Ws")[0],	&ShadowRatio_ParamsMap.at("Sigmas")[0],	0,	&ShadowRatio_ParamsMap.at("Sigmas_Err")[0]	);
		gae 	= new TGraphAsymmErrors(ShadowRatio_ParamsMap.at("Ws").size(),
												&ShadowRatio_ParamsMap.at("Ws")[0],        &ShadowRatio_ParamsMap.at("Sigmas")[0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&Sigmas_TotalSysErr[0], &Sigmas_TotalSysErr[0]);

		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");
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
	Sigma_ALICE_2019_Strategy(const AnalysisData& data_) : ResultPlotStrategy{data_} {}

	void Apply()
	{
		gae ->SetMarkerStyle(24);
		gae ->SetFillColorAlpha(16, 0.7);
		gae ->SetFillStyle(1001);
		gae ->Draw("2same");
		ge  ->SetMarkerStyle(25);
		ge  ->SetMarkerColor(4);
		ge  ->SetMarkerSize(1.5);
		ge  ->SetLineColor(4);
		ge  ->SetLineWidth(2);
		ge  ->Draw("pezsame");
	}
};

