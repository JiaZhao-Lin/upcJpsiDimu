#include "PdfFactory.h"
#include "HistWorker.h"

Double_t DoubleCrystalBall(Double_t x, Double_t N, Double_t mu, Double_t sigma, 
						Double_t cbNL, Double_t cbAlphaL, Double_t cbNR, Double_t cbAlphaR)
{

	Double_t A = TMath::Power(cbNL/fabs(cbAlphaL), cbNL) * TMath::Exp(-cbAlphaL*cbAlphaL/2.);
	Double_t B = cbNL/fabs(cbAlphaL) - fabs(cbAlphaL);

	Double_t C = TMath::Power(cbNR/fabs(cbAlphaR), cbNR) * TMath::Exp(-cbAlphaR*cbAlphaR/2.);
	Double_t D = cbNR/fabs(cbAlphaR) - fabs(cbAlphaR);

	Double_t norm = (x-mu)/sigma;

	if(norm < -cbAlphaL) 
	{
		return N * A * TMath::Power(B-norm, -cbNL);
	}
	else if(norm < cbAlphaR) 
	{
		return N * TMath::Exp(-0.5*norm*norm);
	}
	else 
	{
		return N * C * TMath::Power(D+norm, -cbNR);
	}
}

Double_t DoubleCrystalBall_Asym(Double_t x, Double_t N, Double_t mu, Double_t sigmaL, Double_t sigmaR,
						Double_t cbNL, Double_t cbAlphaL, Double_t cbNR, Double_t cbAlphaR)
{

	Double_t A = TMath::Power(cbNL/fabs(cbAlphaL), cbNL) * TMath::Exp(-cbAlphaL*cbAlphaL/2.);
	Double_t B = cbNL/fabs(cbAlphaL) - fabs(cbAlphaL);

	Double_t C = TMath::Power(cbNR/fabs(cbAlphaR), cbNR) * TMath::Exp(-cbAlphaR*cbAlphaR/2.);
	Double_t D = cbNR/fabs(cbAlphaR) - fabs(cbAlphaR);

	Double_t normL = (x-mu)/sigmaL;
	Double_t normR = (x-mu)/sigmaR;

	if(normL < -cbAlphaL) 
	{
		return N * A * TMath::Power(B-normL, -cbNL);
	}
	else if (normL < 0)
	{
		return N * TMath::Exp(-0.5*normL*normL);
	}
	else if(normR < cbAlphaR) 
	{
		return N * TMath::Exp(-0.5*normR*normR);
	}
	else 
	{
		return N * C * TMath::Power(D+normR, -cbNR);
	}
}


struct QEDPdf : public PdfFactory
{
	//Priviate but not so priviate members-----------------------------------------------
	enum QEDfuction { Cubic };
	TF1 * fQED;
	const double massLow4Fit, massHig4Fit;
	const int FitN;
	RooRealVar  & mMass;
	RooConstVar mP0, mP1, mP2, mP3;
	// RooRealVar  mP0, mP1, mP2, mP3;
	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	QEDPdf(TH1D* Hist, RooRealVar & mMass, const double massLow4Fit, const double massHig4Fit, const int FitN) : 
			PdfFactory{Hist}, mMass{mMass}, massLow4Fit{massLow4Fit}, massHig4Fit{massHig4Fit}, FitN{FitN}	{}
	virtual ~QEDPdf() = default;
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual RooGenericPdf * GetPdf()		{if(!Pdf)	throw std::runtime_error("QEDPdf ----> No PDF!!!");	return Pdf;	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	void Init( QEDfuction f = QEDfuction::Cubic )
	{
		switch(f)
		{
			case Cubic 		:	fQED = new TF1("fQED", fReject4QED, massLow4Fit, massHig4Fit, 4);	break;
			default			:	fQED = new TF1("fQED", fReject4QED, massLow4Fit, massHig4Fit, 4);	break;
		}

		for (int i = 0; i < FitN; ++i)
		{
			Hist ->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters
		}

		mP0	=	RooConstVar(  "mP0", "mP0",  fQED->GetParameter(0));
		mP1	=	RooConstVar(  "mP1", "mP1",  fQED->GetParameter(1));
		mP2	=	RooConstVar(  "mP2", "mP2",  fQED->GetParameter(2));
		mP3	=	RooConstVar(  "mP3", "mP3",  fQED->GetParameter(3));
		// mP0	=	RooRealVar(  "mP0", "mP0",  fQED->GetParameter(0), -1.e8,  1.e8);
		// mP1	=	RooRealVar(  "mP1", "mP1",  fQED->GetParameter(1), 0.,     1.e5);
		// mP2	=	RooRealVar(  "mP2", "mP2",  fQED->GetParameter(2), -1.e5,    0.);
		// mP3	=	RooRealVar(  "mP3", "mP3",  fQED->GetParameter(3), 0.,    1.e4);

		Pdf = new RooGenericPdf("qedPdf", "qedPdf", "mP0 + mP1*mMass + mP2*mMass*mMass + mP3*mMass*mMass*mMass", 
									RooArgSet(mP0, mP1, mP2, mP3, mMass));
	}
	double GetInitN(const double BinLow, const double BinHigh)
	{
		const int tem_QEDBinLow = HistWorker::FindBin(Hist, BinLow, 0);
		const int tem_QEDBinHig = HistWorker::FindBin(Hist, BinHigh, 1);
		double N = Hist->Integral(tem_QEDBinLow, tem_QEDBinHig)*(massHig4Fit-massLow4Fit)/(BinHigh-BinLow);

		return N;
	}
	//-----------------------------------------------------------------------------------
	protected:
	static double fReject4QED(double *x, double *par)
	{
		if((x[0]>2.80&& x[0]<3.35) || (x[0]>3.50 && x[0]<3.85))
		{
			TF1::RejectPoint();
			return 0;
		}

		return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3);
	}
};

struct JpsiPdf : public PdfFactory
{	
	//Priviate but not so priviate members-----------------------------------------------
	const double massLow4Fit, massHig4Fit;
	RooRealVar  & mMass;
	RooAddPdf 		*jpsiRooCrystalBallPdf;
	RooGenericPdf 	*jpsiDoubleCrystalBallPdf;
	RooGenericPdf	*jpsiDoubleCrystalBall_AsymPdf;

	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	JpsiPdf(TH1D* Hist, RooRealVar & mMass, const double massLow4Fit, const double massHig4Fit) : 
			PdfFactory{Hist}, mMass{mMass}, massLow4Fit{massLow4Fit}, massHig4Fit{massHig4Fit}	{}
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual RooGenericPdf * GetPdf()		{if(!Pdf)	throw std::runtime_error("JpsiPdf ----> No PDF!!!");	return Pdf;	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	void Init(RooRealVar& cbAlpha, RooRealVar& cbN, RooRealVar& jpsiSigma, RooRealVar& jpsiMu)
	{
		Pdf = new RooGenericPdf("jpsiPdf", "jpsiPdf",
					"ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma,jpsiMu)", 
					RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, jpsiMu));
	}

	void InitDoubleCrystalBall(RooRealVar& jpsiN, RooRealVar& jpsiMu, RooRealVar& jpsiSigma, RooRealVar& sigmaRatio, RooRealVar& cbNL, RooRealVar& cbAlphaL, RooRealVar& cbNR, RooRealVar& cbAlphaR)
	{
		jpsiDoubleCrystalBallPdf = new RooGenericPdf("jpsiDoubleCrystalBallPdf", "jpsiDoubleCrystalBallPdf",
					"DoubleCrystalBall(mMass, jpsiN, jpsiMu, jpsiSigma, cbNL, cbAlphaL, cbNR, cbAlphaR)", 
					RooArgSet(mMass, jpsiN, jpsiMu, jpsiSigma, cbNL, cbAlphaL, cbNR, cbAlphaR));
	}

	void InitDoubleCrystalBall_Asym(RooRealVar& jpsiN, RooRealVar& jpsiMu, RooRealVar& jpsiSigmaL, RooRealVar& jpsiSigmaR, RooRealVar& sigmaRatio, RooRealVar& cbNL, RooRealVar& cbAlphaL, RooRealVar& cbNR, RooRealVar& cbAlphaR)
	{
		jpsiDoubleCrystalBall_AsymPdf = new RooGenericPdf("jpsiDoubleCrystalBall_AsymPdf", "jpsiDoubleCrystalBall_AsymPdf",
					"DoubleCrystalBall_Asym(mMass, jpsiN, jpsiMu, jpsiSigmaL, jpsiSigmaR, cbNL, cbAlphaL, cbNR, cbAlphaR)", 
					RooArgSet(mMass, jpsiN, jpsiMu, jpsiSigmaL, jpsiSigmaR, cbNL, cbAlphaL, cbNR, cbAlphaR));
	}

	// RooAddPdf* 		GetRooCrystalBallPdf() 			{if(!jpsiRooCrystalBallPdf)			throw std::runtime_error("JpsiPdf ----> No jpsiRooCrystalBallPdf!!!");	return jpsiRooCrystalBallPdf;}
	RooGenericPdf* 	GetDoubleCrystalBallPdf()		{if(!jpsiDoubleCrystalBallPdf)		throw std::runtime_error("JpsiPdf ----> No DoubleCrystalBallPdf!!!");	return jpsiDoubleCrystalBallPdf;}
	RooGenericPdf* 	GetDoubleCrystalBall_AsymPdf()	{if(!jpsiDoubleCrystalBall_AsymPdf)	throw std::runtime_error("JpsiPdf ----> No DoubleCrystalBall_Asym!!!");	return jpsiDoubleCrystalBall_AsymPdf;}


	double GetInitN(const double BinLow, const double BinHigh, const double nQED4Init)
	{
		const int    tem_JpsiBinLow = HistWorker::FindBin(Hist, BinLow,  0);
		const int    tem_JpsiBinHig = HistWorker::FindBin(Hist, BinHigh, 1);

		//const double nJpsi4Init     = Hist->Integral(tem_JpsiBinLow, tem_JpsiBinHig)*0.80; //- nQED4Init*(3.30-2.80)/(massHig4Fit-massLow4Fit);
		double N 	= Hist->Integral(tem_JpsiBinLow, tem_JpsiBinHig) - nQED4Init*(BinHigh-BinLow)/(massHig4Fit-massLow4Fit);

		return N;
	}
	//-----------------------------------------------------------------------------------
};

struct PsiPdf : public PdfFactory
{	
	//Priviate but not so priviate members-----------------------------------------------
	const double massLow4Fit, massHig4Fit;
	RooRealVar  & mMass;
	RooConstVar massRatio = RooConstVar(  "massRatio",   "massRatio",  mPsi_PDG/mJpsi_PDG);
	RooAddPdf 		*psiRooCrystalBallPdf;
	RooGenericPdf 	*psiDoubleCrystalBallPdf;
	RooGenericPdf 	*psiDoubleCrystalBall_AsymPdf;
	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	PsiPdf(TH1D* Hist, RooRealVar & mMass, const double massLow4Fit, const double massHig4Fit) : 
			PdfFactory{Hist}, mMass{mMass}, massLow4Fit{massLow4Fit}, massHig4Fit{massHig4Fit}	{}
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual RooGenericPdf * GetPdf()		{if(!Pdf)	throw std::runtime_error("PsiPdf ----> No PDF!!!");	return Pdf;	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	void Init(RooRealVar& cbAlpha, RooRealVar& cbN, RooRealVar& jpsiSigma, RooRealVar& jpsiMu)
	{
		Pdf = new RooGenericPdf("psiPdf",  "psiPdf",
					"ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*massRatio,jpsiMu*massRatio)", 
					RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, jpsiMu, massRatio)); // psiMu = jpsiMu * massRatio; psiSigma = jpsiSigma * massRatio
	}

	void InitDoubleCrystalBall(RooRealVar& psiN, RooRealVar& jpsiMu, RooRealVar& jpsiSigma, RooRealVar& sigmaRatio, RooRealVar& cbNL, RooRealVar& cbAlphaL, RooRealVar& cbNR, RooRealVar& cbAlphaR)
	{
		psiDoubleCrystalBallPdf = new RooGenericPdf("psiDoubleCrystalBallPdf", "psiDoubleCrystalBallPdf",
					"DoubleCrystalBall(mMass, psiN, jpsiMu*massRatio, jpsiSigma*massRatio, cbNL, cbAlphaL, cbNR, cbAlphaR)", 
					RooArgSet(mMass, psiN, jpsiMu, jpsiSigma, cbNL, cbAlphaL, cbNR, cbAlphaR, massRatio));
	}

	void InitDoubleCrystalBall_Asym(RooRealVar& psiN, RooRealVar& jpsiMu, RooRealVar& jpsiSigmaL, RooRealVar& jpsiSigmaR, RooRealVar& sigmaRatio, RooRealVar& cbNL, RooRealVar& cbAlphaL, RooRealVar& cbNR, RooRealVar& cbAlphaR)
	{
		//psiDoubleCrystalBall_AsymPdf = new RooGenericPdf("psiDoubleCrystalBall_AsymPdf", "psiDoubleCrystalBall_AsymPdf",
		//			"DoubleCrystalBall_Asym(mMass, psiN, jpsiMu*massRatio, jpsiSigmaL*massRatio, jpsiSigmaR*massRatio, cbNL, cbAlphaL, cbNR, cbAlphaR)", 
		//			RooArgSet(mMass, psiN, jpsiMu, jpsiSigmaL, jpsiSigmaR, cbNL, cbAlphaL, cbNR, cbAlphaR, massRatio));
	}

	// RooAddPdf* 		GetPdfRooCrystalBall()			{if(!psiRooCrystalBallPdf)			throw std::runtime_error("PsiPdf ----> No RooCrystalBallPdf!!!");			return psiRooCrystalBallPdf;}
	RooGenericPdf* 	GetDoubleCrystalBallPdf()		{if(!psiDoubleCrystalBallPdf)		throw std::runtime_error("PsiPdf ----> No DoubleCrystalBallPdf!!!");		return psiDoubleCrystalBallPdf;}
	RooGenericPdf* 	GetDoubleCrystalBall_AsymPdf()	{if(!psiDoubleCrystalBall_AsymPdf)	throw std::runtime_error("PsiPdf ----> No psiDoubleCrystalBall_AsymPdf!!!");return psiDoubleCrystalBall_AsymPdf;}
	//-----------------------------------------------------------------------------------
};

struct TotalMassPdf: public PdfFactory
{
	//Priviate but not so priviate members-----------------------------------------------
	//-----------------------------------------------------------------------------------
	//Constructor------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	//Virtual Functions------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	//Free Functions---------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
};

struct PtPdf : public PdfFactory
{
	//Priviate but not so priviate members-----------------------------------------------
	//-----------------------------------------------------------------------------------
	//Constructor------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	//Virtual Functions------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	//Free Functions---------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
};
