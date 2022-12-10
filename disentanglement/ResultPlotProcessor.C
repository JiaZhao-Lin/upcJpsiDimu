#include "FigureFrameStrategy.C"
#include "ResultPlotStrategy.C"

struct ResultPlotProcessor
{
	std::unique_ptr<FigureFrameStrategy> frame;
	std::vector< std::unique_ptr<ResultPlotStrategy> > plots;

	void SaveAs(TString fileName){	frame->SaveAs(fileName);	};

	void Process()
	{
		frame->Apply();
		for(auto &plot : plots) {	plot->Apply(frame->leg);	}
		frame->DrawLegend();
	}

	void SetFigureFrame(FigureFrameStrategyList frame_)
	{
		switch (frame_)
		{
		case FigureFrameStrategyList::Sigma_Log:
			frame = std::make_unique<	Sigma_Log_Strategy		>();
			break;
		case FigureFrameStrategyList::Sigma_LogLog:
			frame = std::make_unique<	Sigma_LogLog_Strategy	>();
			break;
		case FigureFrameStrategyList::R:
			frame = std::make_unique<	R_Strategy				>();
			break;

		case FigureFrameStrategyList::DSigmaDy:
			frame = std::make_unique<	DSigmaDy_Strategy		>();
			break;
		case FigureFrameStrategyList::DSigmaDy_NeuConfig:
			frame = std::make_unique<	DSigmaDy_NeuConfig_Strategy	>();
			break;
		
		default:
			std::runtime_error("FigureFrameStrategy Not Found!");
			break;
		}
	}

	void AddPlot(const AnalysisData& data_,	PlotStrategyList plot_)
	{
		switch (plot_)
		{
		case PlotStrategyList::Sigma_CMS:
			plots.push_back( std::make_unique<	Sigma_CMS_PlotStrategy		>(data_) );
			break;
		case PlotStrategyList::Sigma_ALICE_2019:
			plots.push_back( std::make_unique<	Sigma_ALICE_2019_Strategy	>(data_) );
			break;
		case PlotStrategyList::Sigma_ALICE_2021:
			plots.push_back( std::make_unique<	Sigma_ALICE_2021_Strategy	>(data_) );
			break;
		case PlotStrategyList::Sigma_LHCb_2022:
			plots.push_back( std::make_unique<	Sigma_LHCb_2022_Strategy	>(data_) );
			break;

		case PlotStrategyList::R_CMS:
			plots.push_back( std::make_unique<	R_CMS_Strategy				>(data_) );
			break;
		case PlotStrategyList::R_ALICE_2019:
			plots.push_back( std::make_unique<	R_ALICE_2019_Strategy		>(data_) );
			break;
		case PlotStrategyList::R_ALICE_2021:
			plots.push_back( std::make_unique<	R_ALICE_2021_Strategy		>(data_) );
			break;
		case PlotStrategyList::R_LHCb_2022:
			plots.push_back( std::make_unique<	R_LHCb_2022_Strategy		>(data_) );
			break;

		case PlotStrategyList::DSigmaDy_CMS:
			plots.push_back( std::make_unique<	DSigmaDy_CMS_Strategy		>(data_) );
			break;
		case PlotStrategyList::DSigmaDy_ALICE_2019:
			plots.push_back( std::make_unique<	DSigmaDy_ALICE_2019_Strategy>(data_) );
			break;
		case PlotStrategyList::DSigmaDy_ALICE_2021:
			plots.push_back( std::make_unique<	DSigmaDy_ALICE_2021_Strategy>(data_) );
			break;
		case PlotStrategyList::DSigmaDy_LHCb_2022:
			plots.push_back( std::make_unique<	DSigmaDy_LHCb_2022_Strategy	>(data_) );
			break;
		default:
			throw std::runtime_error("PlotStrategyList Not Found!");
			break;
		}
	}
};

// void ResultPlotProcessor()
// {
// 	AnalysisData Data_CMS("CMS");
// 	AnalysisDataObserver obs;
// 	Data_CMS.Subscribe(&obs);

// 	Data_CMS.LoadHist("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins_NewCohJpsi.appliedTnP.root");

// 	struct ResultPlotProcessor p(Data_CMS);
// 	p.SetFigureFrame(FigureFrameStrategyList::Sigma_Log);
// 	p.AddPlot(PlotStrategyList::Sigma_CMS);
// 	p.Process();
// }