#include "LoadSignal.h"
#include "constants.h"

struct LoadMvsPtvsRap_NeuDir : public LoadSignal
{
	//Priviate but not so priviate members-----------------------------------------------
	TH3D * hMvsPtvsRap;	//For AnAn
	TH3D * hMvsPtvsRap_NeuDir[nNeus][nNeus];
	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	LoadMvsPtvsRap_NeuDir(TString infileDir) : LoadSignal{infileDir} {	Read();	}
	virtual ~LoadMvsPtvsRap_NeuDir() = default;
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual Bool_t Read() override
	{	

		if( nNeus !=2 ) 
		{
			throw std::runtime_error("LoadMvsPtvsRap_NeuDir ----> Not running the correct 0nXn!!!");
			return kFALSE;
		}

		cout<<"------>START Reading LoadMvsPtvsRap_NeuDir: "<<infile->GetName()<<endl;

		//Getting AnAn
		hMvsPtvsRap = (TH3D*) infile->Get("hMvsPtvsRap");
		//Getting 0n0n, 0nXn, Xn0n, and XnXn
		for (int ip = 0; ip < nNeus; ip++) 
		{
			for (int im = 0; im < nNeus; im++) 
			{
				hMvsPtvsRap_NeuDir[ip][im]     = (TH3D*)infile->Get( Form("hMvsPtvsRap_NeuDir%dp%dm",     ip, im) );
				
				cout<<"Reading: "<<hMvsPtvsRap_NeuDir[ip][im]->GetName()     <<endl;
				cout<<"Entries: "<<hMvsPtvsRap_NeuDir[ip][im]->GetEntries() <<endl;
			}
		}

		cout<<"------>DONE  Reading LoadMvsPtvsRap_NeuDir: "<<infile->GetName()<<endl;
		return kTRUE;
	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	TH3D * GetHist(const int idx_1, const int idx_2) const {	return hMvsPtvsRap_NeuDir[idx_1][idx_2];	}
	TH3D * GetHist() 								 const {	return hMvsPtvsRap;	}
};


struct LoadEfficiency : public LoadSignal
{
	//Priviate but not so priviate members-----------------------------------------------
	TH1D* H_EffVsY_CohJpsi;
	TH1D* H_EffVsY_CohPsi;
	TH1D* H_EffVsY_CohPsi2Jpsi;
	const Bool_t SymmetricRapBin;
	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	LoadEfficiency(TString infileDir, const Bool_t SymmetricRapBin) : LoadSignal{infileDir}, 
																	SymmetricRapBin{SymmetricRapBin} {	Read();	}
	virtual ~LoadEfficiency() = default;
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual Bool_t Read() override
	{	
		cout<<"------>START Reading LoadEfficiency: "<<infile->GetName()<<endl;

		if(SymmetricRapBin){
			H_EffVsY_CohJpsi     = (TH1D*) infile->Get("hEffvsRap_Symm_CohJpsi");        //hEffvsRap_CohJpsi_0n0n, 0nXn, XnXn
			H_EffVsY_CohPsi      = (TH1D*) infile->Get("hEffvsRap_Symm_CohPsi2S");
			H_EffVsY_CohPsi2Jpsi = (TH1D*) infile->Get("hEffvsRap_Symm_CohPsi2SFeeddown");
		}
		else
		{
			H_EffVsY_CohJpsi     = (TH1D*) infile->Get("hEffvsRap_CohJpsi");        //hEffvsRap_CohJpsi_0n0n, 0nXn, XnXn
			H_EffVsY_CohPsi      = (TH1D*) infile->Get("hEffvsRap_CohPsi2S");
			H_EffVsY_CohPsi2Jpsi = (TH1D*) infile->Get("hEffvsRap_CohPsi2SFeeddown");
		}

		cout<<"------>DONE  Reading LoadEfficiency: "<<infile->GetName()<<endl;
		return kTRUE;
	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	TH1D* GetCohJpsi()		{if(!H_EffVsY_CohJpsi) 		std::runtime_error("LoadEfficiency: Empty!!!"); return H_EffVsY_CohJpsi;}
	TH1D* GetCohPsi()		{if(!H_EffVsY_CohPsi) 		std::runtime_error("LoadEfficiency: Empty!!!"); return H_EffVsY_CohPsi; }
	TH1D* GetCohPsi2Jpsi()	{if(!H_EffVsY_CohPsi2Jpsi) 	std::runtime_error("LoadEfficiency: Empty!!!"); return H_EffVsY_CohPsi2Jpsi;}
};


struct LoadAcceptance : public LoadSignal
{
	//Priviate but not so priviate members-----------------------------------------------
	TH1D* H_AccVsY_CohJpsi;

	const Bool_t SymmetricRapBin;
	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	LoadAcceptance(TString infileDir, const Bool_t SymmetricRapBin) : LoadSignal{infileDir}, 
																	SymmetricRapBin{SymmetricRapBin} {	Read();	}
	virtual ~LoadAcceptance() = default;
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual Bool_t Read() override
	{	
		cout<<"------>START Reading LoadAcceptance: "<<infile->GetName()<<endl;

		if(SymmetricRapBin)		H_AccVsY_CohJpsi     = (TH1D*) infile->Get("hAccvsRap_Symm_CohJpsi");
		else 					H_AccVsY_CohJpsi     = (TH1D*) infile->Get("hAccvsRap_CohJpsi");

		cout<<"------>DONE  Reading LoadAcceptance: "<<infile->GetName()<<endl;
		return kTRUE;
	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	TH1D* GetAcceptance()	{if(!H_AccVsY_CohJpsi)	std::runtime_error("LoadAcceptance: Empty!!!"); return H_AccVsY_CohJpsi;}
};


struct LoadJpsiXsec : public LoadSignal
{
	//Priviate but not so priviate members-----------------------------------------------
	TH1D* V_JpsiXsec[6];
	std::map<TString, std::vector<double>> Map_Xsec;
	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	LoadJpsiXsec(TString infileDir) : LoadSignal{infileDir}	{	Read();	}
	virtual ~LoadJpsiXsec() = default;
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual Bool_t Read() override
	{	
		cout<<"------>START Reading JpsiXsec: "<<infile->GetName()<<endl;

		V_JpsiXsec[0] = (TH1D*) infile->Get("hAnAn") ;
		V_JpsiXsec[1] = (TH1D*) infile->Get("hOnOn") ;
		V_JpsiXsec[2] = (TH1D*) infile->Get("hOnXn") ;
		V_JpsiXsec[3] = (TH1D*) infile->Get("hXnOn") ;
		V_JpsiXsec[4] = (TH1D*) infile->Get("hOnXnSum") ;
		V_JpsiXsec[5] = (TH1D*) infile->Get("hXnXn") ;

		std::vector<double> Xsec_0n0n, Xsec_0nXnSum, Xsec_XnXn, Rap;
		std::vector<double> XsecErr_0n0n, XsecErr_0nXnSum, XsecErr_XnXn, RapErr;
		for (int iy = nDiffRapBins/2 + 2 ; iy < nDiffRapBins + 2 ; ++iy)
		{
			Rap.push_back(	V_JpsiXsec[5]->GetBinCenter(iy)	); 				RapErr.push_back(		V_JpsiXsec[5]->GetBinWidth(iy)/2	);
			Xsec_0n0n.push_back(	V_JpsiXsec[1]->GetBinContent(iy)	);	XsecErr_0n0n.push_back(		V_JpsiXsec[1]->GetBinError(iy)	);
			Xsec_0nXnSum.push_back(	V_JpsiXsec[4]->GetBinContent(iy)	);	XsecErr_0nXnSum.push_back(	V_JpsiXsec[4]->GetBinError(iy)	);
			Xsec_XnXn.push_back(	V_JpsiXsec[5]->GetBinContent(iy)	);	XsecErr_XnXn.push_back(		V_JpsiXsec[5]->GetBinError(iy)	);
		}
		Map_Xsec["Rap"] = Rap; 						Map_Xsec["RapErr"] = RapErr;
		Map_Xsec["Xsec_0n0n"] = Xsec_0n0n; 			Map_Xsec["XsecErr_0n0n"] = XsecErr_0n0n;
		Map_Xsec["Xsec_0nXnSum"] = Xsec_0nXnSum; 	Map_Xsec["XsecErr_0nXnSum"] = XsecErr_0nXnSum;
		Map_Xsec["Xsec_XnXn"] = Xsec_XnXn; 			Map_Xsec["XsecErr_XnXn"] = XsecErr_XnXn;

		cout<<"------>DONE  Reading JpsiXsec: "<<infile->GetName()<<endl;
		return kTRUE;
	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	TH1D* GetJpsiXsec(const int i)	{if(!V_JpsiXsec[i])	std::runtime_error("LoadJpsiXsec: Empty!!!"); return V_JpsiXsec[i];}
	std::map<TString, std::vector<double>> GetMap()	const {if(!Map_Xsec.size())	std::runtime_error("LoadJpsiXsec: Empty!!!"); 	return Map_Xsec;}
	std::vector<double> GetMapElement(const TString n)	{if(!Map_Xsec.size())	std::runtime_error("LoadJpsiXsec: Empty!!!");	return Map_Xsec[n];}
	double GetMapElementVal(const TString n, const int i)	{if(!Map_Xsec.size())	std::runtime_error("LoadJpsiXsec: Empty!!!");	return Map_Xsec[n][i];}
};

struct LoadUncertainty : public LoadSignal
{
	//Priviate but not so priviate members-----------------------------------------------
	TH1D* V_TotalSysUncer[4];
	std::map<TString, std::vector<double>> Map_Uncer;
	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	LoadUncertainty(TString infileDir) : LoadSignal{infileDir}	{	Read();	}
	virtual ~LoadUncertainty() = default;
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual Bool_t Read() override
	{	
		cout<<"------>START Reading Uncertainty: "<<infile->GetName()<<endl;

		V_TotalSysUncer[0] = (TH1D*) infile->Get("hTotalSysSigmas") ;
		V_TotalSysUncer[1] = (TH1D*) infile->Get("hTotalSysSigmasUncerPerc") ;
		V_TotalSysUncer[2] = (TH1D*) infile->Get("hTotalSysR") ;
		V_TotalSysUncer[3] = (TH1D*) infile->Get("hTotalSysRUncerPerc") ;

		std::vector<double> Sigmas, R;
		std::vector<double> TotalSysSigmas_Err, TotalSysR_Err;
		for (int i = 0 ; i < V_TotalSysUncer[0]->GetNbinsX() ; ++i)
		{
			Sigmas 				.push_back( V_TotalSysUncer[0]->GetBinContent(i+1) );
			TotalSysSigmas_Err	.push_back( V_TotalSysUncer[0]->GetBinError(i+1) );
			R 					.push_back( V_TotalSysUncer[2]->GetBinContent(i+1) );
			TotalSysR_Err 		.push_back( V_TotalSysUncer[2]->GetBinError(i+1) );
		}
		Map_Uncer["Sigmas"] = Sigmas; 					Map_Uncer["Sigmas_Err"] = TotalSysSigmas_Err;
		Map_Uncer["R"] 		= R; 						Map_Uncer["R_Err"] 	  = TotalSysR_Err;

		cout<<"------>DONE  Reading Uncertainty: "<<infile->GetName()<<endl;
		return kTRUE;
	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	TH1D* GetUncertainty(const int i)	{if(!V_TotalSysUncer[i])	std::runtime_error("LoadUncertainty: Empty!!!"); return V_TotalSysUncer[i];}
	std::map<TString, std::vector<double>> GetMap()	const {if(!Map_Uncer.size())	std::runtime_error("LoadUncertainty: Empty!!!"); 	return Map_Uncer;}
	std::vector<double> GetMapElement(const TString n)	{if(!Map_Uncer.size())	std::runtime_error("LoadUncertainty: Empty!!!");	return Map_Uncer[n];}
	double GetMapElementVal(const TString n, const int i)	{if(!Map_Uncer.size())	std::runtime_error("LoadUncertainty: Empty!!!");	return Map_Uncer[n][i];}
};