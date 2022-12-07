#include "FigureFrameStrategy.C"
#include "ResultPlotStrategy.C"

struct ResultPlotProcessor
{
	std::unique_ptr<FigureFrameStrategy> frame;
	std::vector< std::unique_ptr<ResultPlotStrategy> > plots;
	const AnalysisData& data;

	ResultPlotProcessor(const AnalysisData& data_) : data{data_} {};

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
		default:
			std::runtime_error("FigureFrameStrategy Not Found!");
			break;
		}
	}

	void AddPlot(PlotStrategyList plot_)
	{
		switch (plot_)
		{
		case PlotStrategyList::Sigma_CMS:
			plots.push_back( std::make_unique<	Sigma_CMS_PlotStrategy		>(data) );
			break;
		case PlotStrategyList::Sigma_ALICE_2019:
			plots.push_back( std::make_unique<	Sigma_ALICE_2019_Strategy	>(data) );
			break;
		default:
			std::runtime_error("PlotStrategyList Not Found!");
			break;
		}
	}
};

void ResultPlotProcessor()
{
	AnalysisData Data_CMS("CMS");
	AnalysisDataObserver obs;
	Data_CMS.Subscribe(&obs);

	Data_CMS.LoadHist("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins_NewCohJpsi.appliedTnP.root");

	struct ResultPlotProcessor p(Data_CMS);
	p.SetFigureFrame(FigureFrameStrategyList::Sigma_Log);
	p.AddPlot(PlotStrategyList::Sigma_CMS);
	
}