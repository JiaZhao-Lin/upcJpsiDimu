#ifndef DisentangleFwdRapSigmaAnalyzer_H
#define DisentangleFwdRapSigmaAnalyzer_H

//For the disentanglement of the forward rapidity cross section
struct DisentangleFwdRapSigmaAnalyzer: Analyzer
{
	bool    print               = false;
	double  fluxScaleFactor    	= 1.0;
	class PhotonFluxReader cPhotonFluxReader;

	DisentangleFwdRapSigmaAnalyzer(AnalysisData& data_, TString inFileDir_, TString subCase_) : Analyzer{data_}, cPhotonFluxReader{data_, inFileDir_, subCase_} {};

	void PrintHandling() const
	{
		cout << endl << std::setfill('+') << std::setw(100) << Form("%s::Handle Is Handling The Business... ", typeid(this).name()) << endl;
		cout << setfill(' ');
	}

	void Calculate()
	{
		std::vector<double> Sigma, Sigma_Err, Sigma_SysErrLow, Sigma_SysErrHigh, Flux_dNdy_AnAn;
		Flux_dNdy_AnAn  = data.Get("dNdy_AnAn");

		for (int i = 0; i < data.GetSize("Dy"); ++i)
		{
			Sigma.push_back(	data.Get("DSigmaDy_AnAn", i) / data.Get("dNdy_AnAn", i)	);
		}

		for (int i = 0; i < data.GetSize("Dy"); ++i)
		{
			Sigma_Err		.push_back(	Sigma[i] * data.Get("DSigmaDy_AnAn_Err", i) / data.Get("DSigmaDy_AnAn", i)	);
			Sigma_SysErrLow	.push_back(	Sigma[i] * TMath::Hypot(	data.Get("DSigmaDy_AnAn_SysErrLow", i) / data.Get("DSigmaDy_AnAn", i),
																	data.Get("dNdy_AnAn_Uncer", i) / 100.0	)	);	// flux err from dNdy Uncer, Converting from % to fraction
			Sigma_SysErrHigh.push_back(	Sigma[i] * TMath::Hypot(	data.Get("DSigmaDy_AnAn_SysErrHigh", i) / data.Get("DSigmaDy_AnAn", i),
																	data.Get("dNdy_AnAn_Uncer", i) / 100.0	)	);	// flux err from dNdy Uncer, Converting from % to fraction
		}		

		data.Add("Sigma"        	,Sigma);
		data.Add("Sigma_Err"    	,Sigma_Err);
		data.Add("Sigma_SysErrLow" 	,Sigma_SysErrLow);
		data.Add("Sigma_SysErrHigh" ,Sigma_SysErrHigh);
	}

	void Handle()	override
	{
		PrintHandling();

		cPhotonFluxReader.SetFluxErr(true);
		cPhotonFluxReader.Handle();
		Calculate();

		Analyzer::Handle();
	}
};

#endif