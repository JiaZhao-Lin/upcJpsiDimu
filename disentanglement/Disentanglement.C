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
#include "MidRapAnalyzer.C"
#include "ResultPlotProcessor.C"

void Disentanglement()
{

	AnalysisData Data_CMS("CMS");
	AnalysisData Data_ALICE_2019("ALICE_2019");
	AnalysisData Data_ALICE_2021("ALICE_2021");
	AnalysisData Data_LHCb_2022("LHCb_2022");

	AnalysisDataObserver obs;
	Data_CMS.Subscribe(&obs);
	// Data_ALICE_2019.Subscribe(&obs);
	// Data_ALICE_2021.Subscribe(&obs);
	// Data_LHCb_2022.Subscribe(&obs);

	Data_CMS.LoadHist("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins_NewCohJpsi.appliedTnP.root");
	Data_ALICE_2019.LoadFile("inFiles/DSigmaDy_ALICE_2019_Selected.txt");
	Data_ALICE_2021.LoadFile("inFiles/DSigmaDy_ALICE_2021_Selected.txt");
	Data_LHCb_2022.LoadFile("inFiles/DSigmaDy_LHCb_2022_Selected.txt");

	DisentanglementAnalyzer AnaRoot_CMS(Data_CMS,	"../simulation/flux/",      "_SigNN68p3R6p67a0p56");
	FwdRapAnalyzer AnaRoot_ALICE_2019(Data_ALICE_2019,	"../simulation/flux/",      "_SigNN68p3R6p67a0p56");
	MidRapAnalyzer AnaRoot_ALICE_2021(Data_ALICE_2021,	"../simulation/flux/",      "_SigNN68p3R6p67a0p56");
	FwdRapAnalyzer AnaRoot_LHCb_2022(Data_LHCb_2022,	"../simulation/flux/",      "_SigNN68p3R6p67a0p56");

	AnaRoot_CMS.Handle();
	AnaRoot_ALICE_2019.Handle();
	AnaRoot_ALICE_2021.Handle();
	AnaRoot_LHCb_2022.Handle();

	// Data_CMS.Print();
	// Data_ALICE_2019.Print();
	// Data_ALICE_2021.Print();
	// Data_LHCb_2022.Print();
	

	struct ResultPlotProcessor p;
	p.SetFigureFrame(FigureFrameStrategyList::R);
	p.AddPlot		(Data_CMS,			PlotStrategyList::R_CMS);
	p.AddPlot		(Data_ALICE_2019,	PlotStrategyList::R_ALICE_2019);
	p.AddPlot		(Data_ALICE_2021,	PlotStrategyList::R_ALICE_2021);
	p.AddPlot		(Data_LHCb_2022,	PlotStrategyList::R_LHCb_2022);
	p.Process		();
	p.SaveAs		("./outFigures/test.pdf");
}