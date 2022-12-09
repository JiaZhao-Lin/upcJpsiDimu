#ifndef DisentangleFwdRapSigmaAnalyzer_H
#define DisentangleFwdRapSigmaAnalyzer_H

struct DisentangleFwdRapSigmaAnalyzer: PhotonFluxAnalyzer
{
	bool    print               = false;
	double  fluxScaleFactor    	= 1.0;

	DisentangleFwdRapSigmaAnalyzer(AnalysisData& data_, TString inFileDir_, TString subCase_) : PhotonFluxAnalyzer(data_, inFileDir_, subCase_) {};

	void PrintHandling() const
	{
		cout << endl << std::setfill('+') << std::setw(100) << Form("%s::Handle Is Handling The Business... ", typeid(this).name()) << endl;
		cout << setfill(' ');
	}

	void Calculate()
	{
		InterpolateFlux();

		std::vector<double> Sigma, Sigma_Err, Flux_dNdy_AnAn;
		Flux_dNdy_AnAn  = data.Get("dNdy_AnAn");

		for (int i = 0; i < data.GetSize("Rap"); ++i)
		{
			Sigma.push_back(	data.Get("DSigmaDy_AnAn", i) / data.Get("dNdy_AnAn", i)	);
		}

		for (int i = 0; i < data.GetSize("Rap"); ++i)
		{
			Sigma_Err.push_back(	Sigma[i] * data.Get("DSigmaDy_AnAn_Err", i) / data.Get("DSigmaDy_AnAn", i)	);
		}		

		data.Add("Sigma"        ,Sigma);
		data.Add("Sigma_Err"    ,Sigma_Err);
	}

	void Handle()	override
	{
		PrintHandling();

		Calculate();

		Analyzer::Handle();
	}
};

#endif