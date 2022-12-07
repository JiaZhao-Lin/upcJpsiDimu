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

#include "ConversionAnalyzer.C"
#include "Fit2DAnalyzer.C"
#include "ShadowingRatioAnalyzer.C"

void Disentanglement()
{

	AnalysisData Data_CMS("CMS");
	AnalysisDataObserver obs;
	Data_CMS.Subscribe(&obs);

	Data_CMS.LoadHist("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins_NewCohJpsi.appliedTnP.root");

	Analyzer                    AnaRoot                 (Data_CMS);
	ConversionAnalyzer          cConversionAnalyzer     (Data_CMS);
	Fit2DAnalyzer               cFit2DAnalyzer          (Data_CMS,      "../simulation/flux/",      "_SigNN68p3R6p67a0p56");
	ShadowingRatioAnalyzer      cShadowingRatioAnalyzer (Data_CMS);

	AnaRoot.Add(&cConversionAnalyzer);
	AnaRoot.Add(&cFit2DAnalyzer);
	AnaRoot.Add(&cShadowingRatioAnalyzer);

	AnaRoot.Handle();
}