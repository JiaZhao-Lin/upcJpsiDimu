/*
Main analysis chain for disentangling cross sections
Use Analyzer for connecting each analysis step 
Use AnalysisData for convenient manipulation of the data 
Use AnalysisDataObserver to notify each data change

Example Usage:
-----------------------------------------
AnalysisData Data_CMS("CMS");
AnalysisDataObserver obs;
Data_CMS.Subscribe(&obs);

Data_CMS.LoadHist(".root");

Analyzer                    AnaRoot                 (Data_CMS);
ConversionAnalyzer          cConversionAnalyzer     (Data_CMS);
-----------------------------------------

Design pattern used: Chain of Responsibility, Observer

Dec. 2022
JiaZhao Lin
*/

#include "DisentanglementAnalyzer.C"
#include "FwdRapAnalyzer.C"
#include "ResultPlotProcessor.C"

void Disentanglement()
{

	AnalysisData Data_CMS("CMS");
	AnalysisData Data_ALICE_2019("ALICE_2019");
	AnalysisDataObserver obs;
	Data_CMS.Subscribe(&obs);
	Data_ALICE_2019.Subscribe(&obs);

	Data_CMS.LoadHist("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins_NewCohJpsi.appliedTnP.root");
	Data_ALICE_2019.LoadFile("inFiles/DSigmaDy_ALICE_2019.txt");

	DisentanglementAnalyzer AnaRoot_CMS(Data_CMS,	"../simulation/flux/",      "_SigNN68p3R6p67a0p56");
	FwdRapAnalyzer AnaRoot_ALICE_2019(Data_ALICE_2019,	"../simulation/flux/",      "_SigNN68p3R6p67a0p56");

	AnaRoot_CMS.Handle();
	AnaRoot_ALICE_2019.Handle();

	// Data_CMS.Print();
	// Data_ALICE_2019.Print();

	

	struct ResultPlotProcessor p;
	p.SetFigureFrame(FigureFrameStrategyList::DSigmaDy);
	p.AddPlot		(Data_CMS,	PlotStrategyList::DSigmaDy_CMS);
	p.AddPlot		(Data_ALICE_2019,	PlotStrategyList::DSigmaDy_ALICE_2019);
	p.Process		();
	p.SaveAs		("./outFigures/test.pdf");
}