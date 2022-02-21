#include "PdfFactory.h"
#include "HistWorker.h"


struct QEDPdf : public PdfFactory
{
	//Priviate but not so priviate members-----------------------------------------------
	enum QEDfuction { Cubic };
	TF1 * fQED;
	const double massLow4Fit, massHig4Fit;
	const int FitN;
	RooRealVar  & mMass;
	RooConstVar mP0, mP1, mP2, mP3;
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
		//RooRealVar mP0(  "mP0", "mP0",  fQED->GetParameter(0), -1.e8,  1.e8);
		//RooRealVar mP1(  "mP1", "mP1",  fQED->GetParameter(1), 0.,     1.e5);
		//RooRealVar mP2(  "mP2", "mP2",  fQED->GetParameter(2), -1.e5,    0.);
		//RooRealVar mP3(  "mP3", "mP3",  fQED->GetParameter(3), 0.,    1.e4);

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
	RooRealVar  cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, gausN;

	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	JpsiPdf(TH1D* Hist, RooRealVar & mMass, const double massLow4Fit, const double massHig4Fit) : 
			PdfFactory{Hist}, mMass{mMass}, massLow4Fit{massLow4Fit}, massHig4Fit{massHig4Fit}	{}
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual RooGenericPdf * GetPdf()		{if(!Pdf)	throw std::runtime_error("JpsiPdf ----> No PDF!!!");	return Pdf;	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	void Init(RooRealVar cbAlpha_, RooRealVar cbN_, RooRealVar jpsiSigma_, RooRealVar sigmaRatio_, RooRealVar jpsiMu_, RooRealVar gausN_)
	{
		cbAlpha = cbAlpha_; cbN = cbN_; jpsiSigma = jpsiSigma_; sigmaRatio = sigmaRatio_; jpsiMu = jpsiMu_; gausN = gausN_;
		Pdf = new RooGenericPdf("jpsiPdf", "jpsiPdf",
					"ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*sigmaRatio,jpsiMu) + gausN*TMath::Gaus(mMass, jpsiMu, jpsiSigma)", 
					RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, gausN));
	}
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
	RooRealVar  cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, gausN;
	RooConstVar massRatio = RooConstVar(  "massRatio",   "massRatio",  mPsi_PDG/mJpsi_PDG);
	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	PsiPdf(TH1D* Hist, RooRealVar & mMass, const double massLow4Fit, const double massHig4Fit) : 
			PdfFactory{Hist}, mMass{mMass}, massLow4Fit{massLow4Fit}, massHig4Fit{massHig4Fit}	{}
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual RooGenericPdf * GetPdf()		{if(!Pdf)	throw std::runtime_error("PsiPdf ----> No PDF!!!");	return Pdf;	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	void Init(RooRealVar cbAlpha_, RooRealVar cbN_, RooRealVar jpsiSigma_, RooRealVar sigmaRatio_, RooRealVar jpsiMu_, RooRealVar gausN_)
	{
		cbAlpha = cbAlpha_; cbN = cbN_; jpsiSigma = jpsiSigma_; sigmaRatio = sigmaRatio_; jpsiMu = jpsiMu_; gausN = gausN_;
		Pdf = new RooGenericPdf("psiPdf",  "psiPdf",
					"ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*sigmaRatio*massRatio,jpsiMu*massRatio) + gausN*TMath::Gaus(mMass, jpsiMu*massRatio, jpsiSigma*massRatio)", 
					RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, massRatio, gausN)); // psiMu = jpsiMu * massRatio; psiSigma = jpsiSigma * massRatio
	}
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