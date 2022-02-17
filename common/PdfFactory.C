#include "PdfFactory.h"
#include "HistWorker.h"


struct QEDPdf : public PdfFactory
{
	//Priviate but not so priviate members-----------------------------------------------
	enum QEDfuction { Cubic };
	TF1 * fQED;
	const double massLow4Fit, massHig4Fit;
	const int FitN;
	RooRealVar  mMass;
	RooConstVar mP0, mP1, mP2, mP3;
	//-----------------------------------------------------------------------------------

	//Constructor------------------------------------------------------------------------
	QEDPdf(TH1D* Hist, const double massLow4Fit, const double massHig4Fit, const int FitN) : 
			PdfFactory{Hist}, massLow4Fit{massLow4Fit}, massHig4Fit{massHig4Fit}, FitN{FitN}	{	Init();	}
	virtual ~QEDPdf() = default;
	//-----------------------------------------------------------------------------------

	//Virtual Functions------------------------------------------------------------------
	virtual RooGenericPdf * GetPdf()		{if(!Pdf)	throw std::runtime_error("QEDPdf ----> No PDF!!!");	return Pdf;	}
	virtual double			GetN() const	{if(!N)		throw std::runtime_error("QEDPdf ----> No N  !!!");	return N;	}
	//-----------------------------------------------------------------------------------

	//Free Functions---------------------------------------------------------------------
	void Init( QEDfuction f = QEDfuction::Cubic )
	{
		switch(f)
		{
			case Cubic 		:	fQED = new TF1("fQED", fReject4QED, massLow4Fit, massHig4Fit, 4);	break;
			default			:	fQED = new TF1("fQED", fReject4QED, massLow4Fit, massHig4Fit, 4);	break;
		}

		mMass = RooRealVar{"mMass", "m_{#mu#mu} (GeV)", massLow4Fit, massHig4Fit};

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

		const int tem_QEDBinLow = HistWorker::FindBin(Hist, 3.30, 0);
		const int tem_QEDBinHig = HistWorker::FindBin(Hist, 3.50, 1);
		N = Hist->Integral(tem_QEDBinLow, tem_QEDBinHig)*(massHig4Fit-massLow4Fit)/(3.50-3.30);

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
	//-----------------------------------------------------------------------------------
	//Constructor------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	//Virtual Functions------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	//Free Functions---------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
};

struct PsiPdf : public PdfFactory
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