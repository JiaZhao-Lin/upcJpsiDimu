#ifndef DisentanglementAnalyzer_H
#define DisentanglementAnalyzer_H

#include "ConversionAnalyzer.C"
#include "Fit2DAnalyzer.C"
#include "ShadowingRatioAnalyzer.C"

struct DisentanglementAnalyzer : Analyzer
{
	ConversionAnalyzer		cConversionAnalyzer;
	Fit2DAnalyzer			cFit2DAnalyzer;
	ShadowingRatioAnalyzer	cShadowingRatioAnalyzer;

    DisentanglementAnalyzer(AnalysisData& data_,	TString inFileDir_,	TString subCase_) : Analyzer{data_},	cConversionAnalyzer{data_},	
																							cFit2DAnalyzer{data_, inFileDir_,	subCase_},	cShadowingRatioAnalyzer{data_}	{};

	void PrintHandling() const
	{
		cout << endl << std::setfill('+') << std::setw(100) << Form("%s::Handle Is Handling The Business... ", typeid(this).name()) << endl;
		cout << setfill(' ');
	}

    void Handle() override
    {
        PrintHandling();
        
		this->Add(&cConversionAnalyzer);
		this->Add(&cFit2DAnalyzer);
		this->Add(&cShadowingRatioAnalyzer);

        Analyzer::Handle();
    }
};


#endif