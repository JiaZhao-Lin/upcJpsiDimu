#ifndef DisentangleMidRapSigmaAnalyzer_H
#define DisentangleMidRapSigmaAnalyzer_H

struct DisentangleMidRapSigmaAnalyzer: PhotonFluxAnalyzer
{
	bool    print               = false;
	double  fluxScaleFactor    	= 1.0;

	DisentangleMidRapSigmaAnalyzer(AnalysisData& data_, TString inFileDir_, TString subCase_) : PhotonFluxAnalyzer(data_, inFileDir_, subCase_) {};

	void PrintHandling() const
	{
		cout << endl << std::setfill('+') << std::setw(100) << Form("%s::Handle Is Handling The Business... ", typeid(this).name()) << endl;
		cout << setfill(' ');
	}

	void Calculate()
	{
		InterpolateFlux();

		std::vector<double> Sigma, Sigma_Err, Sigma_SysErr, Flux_dNdy_AnAn;
		Flux_dNdy_AnAn  = data.Get("dNdy_AnAn");

		for (int i = 0; i < data.GetSize("Dy"); ++i)
		{
			Sigma.push_back(	data.Get("DSigmaDy_AnAn", i) / data.Get("dNdy_AnAn", i) / 2.0	);
		}

		for (int i = 0; i < data.GetSize("Dy"); ++i)
		{
			Sigma_Err	.push_back(	Sigma[i] * data.Get("DSigmaDy_AnAn_Err", i) / data.Get("DSigmaDy_AnAn", i)	);
			Sigma_SysErr.push_back(	Sigma[i] * TMath::Hypot(	data.Get("DSigmaDy_AnAn_Err", 	i)	/ data.Get("DSigmaDy_AnAn", i),	
																0.05	)	);
		}		

		data.Add("Sigma"        ,Sigma);
		data.Add("Sigma_Err"    ,Sigma_Err);
		data.Add("Sigma_SysErr" ,Sigma_SysErr);
	}

	void Handle()	override
	{
		PrintHandling();

		Calculate();

		Analyzer::Handle();
	}
};

#endif