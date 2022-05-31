#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"
#include "../common/LoadSignal.C"
#include "../common/HistWorker.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "RooChi2Var.h"
#include "RooMinimizer.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooFitResult.h"
using namespace RooFit;

int    mTextFont    = 42;
double mTextSize    = 0.045;
int    mTextColor   = 1;

double mMarkerStyle = 20;
double mMarkerSize  = 0.8;

double mTitleSize    = 0.06;
double mXTitleOffset = 0.95;
double mYTitleOffset = 0.95;
double mLabelSize    = 0.05;
double mTickLength   = 0.02;
int    mXNdivisions  = 210;
int    mYNdivisions  = 208;

int    mLineWidth = 2;
int    cohJpsiColor       = kBlue,      cohJpsiStyle      = 1;
int    incohJpsiColor     = kViolet-1,  incohJpsiStyle    = 2;
int    dissoJpsiColor     = kRed,       dissoJpsiStyle    = 2;
int    feeddownJpsiColor  = kAzure+10,  feeddownJpsiStyle = 1;
int    qedColor = kGreen-3, qedStyle = 1;
//------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------
TH1D* hMinY_0nXn[nDiffRapBins];
TH1D* hJpsiPtinY_0nXn[nDiffRapBins];
TH1D* hMinY_Xn0n[nDiffRapBins];
TH1D* hJpsiPtinY_Xn0n[nDiffRapBins];

TH1D* hMinY_0nXnXn0n[nDiffRapBins];
TH1D* hJpsiPtinY_0nXnXn0n[nDiffRapBins];

TH1D* hMinY_SameDir[nDiffRapBins/2];
TH1D* hMinY_OpstDir[nDiffRapBins/2];
TH1D* hJpsiPtinY_SameDir[nDiffRapBins/2];
TH1D* hJpsiPtinY_OpstDir[nDiffRapBins/2];
TH1D* hPureJpsiPtinY_SameDir[nDiffRapBins/2];
TH1D* hPureJpsiPtinY_OpstDir[nDiffRapBins/2];
TH1D* hMAllY_SameDir;
TH1D* hMAllY_OpstDir;
TH1D* hJpsiPtAllY_SameDir;
TH1D* hJpsiPtAllY_OpstDir;
TH1D* hPureJpsiPtAllY_SameDir;
TH1D* hPureJpsiPtAllY_OpstDir;

TH1D* hRatioinY_OpstDir2SameDir[nDiffRapBins/2];
TH1D* hRatioAllY_OpstDir2SameDir;
TH1D* hRatioinY_OpstDir2SameDir_Reb[nDiffRapBins/2];
TH1D* hRatioAllY_OpstDir2SameDir_Reb;

TH1D* hSideBandPtinY_0nXn[nDiffRapBins];
TH1D* hSideBandPtinY_Xn0n[nDiffRapBins];
TH1D* hSideBandPtinY_0nXnXn0n[nDiffRapBins];

TH1D* hSideBandPtinY_SameDir[nDiffRapBins/2];
TH1D* hSideBandPtinY_OpstDir[nDiffRapBins/2];
TH1D* hSideBandPtAllY_SameDir;
TH1D* hSideBandPtAllY_OpstDir;

const int nDirCases = 2;
const TString DirNames[nDirCases] = {"Same direction", "Opposite direction"};
double nQEDinJpsi_iDir_iY[nDirCases][nDiffRapBins/2];
double nQEDinJpsi_SameDir_AllY;
double nQEDinJpsi_OpstDir_AllY;

void prepareData( );
void fitMass4QEDYield( );
void compPtShape( );
void getPtRatios( );

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
void getPtRatio_JpsiIn0nXn( )
{
	prepareData( );

	// fitMass4QEDYield( );
	
	// compPtShape( );
	
	// getPtRatios( );
}
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
void prepareData()
{
	TString inFileDir 		  = "../anaData/jpsiHistos/rawSig.root";

	LoadMvsPtvsRap_NeuDir * hMvsPtvsRap_NeuDir = new LoadMvsPtvsRap_NeuDir(inFileDir);

	TH3D *hMvsPtvsRap_0nXn    = (TH3D*) hMvsPtvsRap_NeuDir->GetHist(1,0)->Clone( "hMvsPtvsRap_0nXn" );
	TH3D *hMvsPtvsRap_Xn0n    = (TH3D*) hMvsPtvsRap_NeuDir->GetHist(0,1)->Clone( "hMvsPtvsRap_Xn0n" );
	
	TH1D* hPtinY_LowMBand_0nXn[nDiffRapBins];
	TH1D* hPtinY_HigMBand_0nXn[nDiffRapBins];
	TH1D* hPtinY_LowMBand_Xn0n[nDiffRapBins];
	TH1D* hPtinY_HigMBand_Xn0n[nDiffRapBins];

	//TH3D *hMvsPtvsRap_0nXnSum = (TH3D*) hMvsPtvsRap_NeuDir[1][0]->Clone( "hMvsPtvsRap_0nXnSum" ); //there is no meaning to add in this way for our interested signals
	//hMvsPtvsRap_0nXnSum       -> Add( hMvsPtvsRap_NeuDir[0][1] );

	for(int iy=0; iy<nDiffRapBins; iy++)
	{
		cout<<"iy: "<<iy<<" "<<mDiffRapLow[iy]<<" <y< "<<mDiffRapHi[iy]<<endl;

		//int ptBinLow    = hMvsPtvsRap_0nXn->GetYaxis()->FindBin( 1.5 + mTinyNum);
		//int nptBinsMax  = hMvsPtvsRap_0nXn->GetYaxis()->FindBin( 3.0 - mTinyNum);
		int ptBinLow    = 1;
		int nptBinsMax  = hMvsPtvsRap_0nXn->GetNbinsY();
		
		int rapBinLow   = HistWorker::FindXBin(hMvsPtvsRap_0nXn, mDiffRapLow[iy], 0);
		int rapBinHi    = HistWorker::FindXBin(hMvsPtvsRap_0nXn, mDiffRapHi[iy] , 1);

		int mJpsiBinLow = HistWorker::FindZBin(hMvsPtvsRap_0nXn, mJpsiMassLow   , 0);
		int mJpsiBinHi  = HistWorker::FindZBin(hMvsPtvsRap_0nXn, mJpsiMassHi    , 1);
		
		hMinY_0nXn[iy]         = (TH1D *)hMvsPtvsRap_0nXn   ->ProjectionZ( Form("hMass_0nXn_iy%d",      iy), rapBinLow, rapBinHi, ptBinLow,    nptBinsMax );
		hJpsiPtinY_0nXn[iy]    = (TH1D *)hMvsPtvsRap_0nXn   ->ProjectionY( Form("hJpsiPt_0nXn_iy%d",    iy), rapBinLow, rapBinHi, mJpsiBinLow, mJpsiBinHi );

		hMinY_Xn0n[iy]         = (TH1D *)hMvsPtvsRap_Xn0n   ->ProjectionZ( Form("hMass_Xn0n_iy%d",      iy), rapBinLow, rapBinHi, ptBinLow,    nptBinsMax );
		hJpsiPtinY_Xn0n[iy]    = (TH1D *)hMvsPtvsRap_Xn0n   ->ProjectionY( Form("hJpsiPt_Xn0n_iy%d",    iy), rapBinLow, rapBinHi, mJpsiBinLow, mJpsiBinHi );

		//------------------------------------------------------------------------------------------------------------------------------------------
		//get sideband
		//------------------------------------------------------------------------------------------------------------------------------------------
		int leftSideMBinLow    = HistWorker::FindZBin(hMvsPtvsRap_0nXn, mLowMassBandLow , 0 );
		int leftSideMBinHi     = HistWorker::FindZBin(hMvsPtvsRap_0nXn, mLowMassBandHi  , 1 );
		int rightSideMBinLow   = HistWorker::FindZBin(hMvsPtvsRap_0nXn, mHiMassBandLow  , 0 );
		int rightSideMBinHi    = HistWorker::FindZBin(hMvsPtvsRap_0nXn, mHiMassBandHi   , 1 );
		
		hPtinY_LowMBand_0nXn[iy] = (TH1D *) hMvsPtvsRap_0nXn->ProjectionY( Form("hPtinY_LowMBand_0nXn_iy%d", iy), rapBinLow, rapBinHi, leftSideMBinLow,  leftSideMBinHi  );
		hPtinY_HigMBand_0nXn[iy] = (TH1D *) hMvsPtvsRap_0nXn->ProjectionY( Form("hPtinY_HigMBand_0nXn_iy%d", iy), rapBinLow, rapBinHi, rightSideMBinLow, rightSideMBinHi );
		hPtinY_LowMBand_Xn0n[iy] = (TH1D *) hMvsPtvsRap_Xn0n->ProjectionY( Form("hPtinY_LowMBand_Xn0n_iy%d", iy), rapBinLow, rapBinHi, leftSideMBinLow,  leftSideMBinHi  );
		hPtinY_HigMBand_Xn0n[iy] = (TH1D *) hMvsPtvsRap_Xn0n->ProjectionY( Form("hPtinY_HigMBand_Xn0n_iy%d", iy), rapBinLow, rapBinHi, rightSideMBinLow, rightSideMBinHi );

		hSideBandPtinY_0nXn[iy] = (TH1D *) hPtinY_LowMBand_0nXn[iy] ->Clone( Form("hSideBandPtinY_0nXn_iy%d", iy) );
		hSideBandPtinY_0nXn[iy] ->Add( hPtinY_HigMBand_0nXn[iy] );
		hSideBandPtinY_Xn0n[iy] = (TH1D *) hPtinY_LowMBand_Xn0n[iy] ->Clone( Form("hSideBandPtinY_Xn0n_iy%d", iy) );
		hSideBandPtinY_Xn0n[iy] ->Add( hPtinY_HigMBand_Xn0n[iy] );
	}//iy

	//combine pt distribution in symmetric y-bins of 0nXn and Xn0n
	for(int iy=0; iy<nDiffRapBins; iy++)
	{
		hMinY_0nXnXn0n[iy]      = (TH1D*) hMinY_0nXn[iy]->Clone( Form("hMinY_0nXnXn0n_iy%d",iy));
		hMinY_0nXnXn0n[iy]      ->Add( hMinY_Xn0n[nDiffRapBins-1-iy] );

		hJpsiPtinY_0nXnXn0n[iy] = (TH1D*) hJpsiPtinY_0nXn[iy]->Clone( Form("hJpsiPtinY_0nXnXn0n_iy%d",iy));
		hJpsiPtinY_0nXnXn0n[iy] ->Add( hJpsiPtinY_Xn0n[nDiffRapBins-1-iy] );
		
		hSideBandPtinY_0nXnXn0n[iy] = (TH1D*) hSideBandPtinY_0nXn[iy]->Clone( Form("hSideBandPtinY_0nXnXn0n_iy%d",iy));
		hSideBandPtinY_0nXnXn0n[iy] ->Add( hSideBandPtinY_Xn0n[nDiffRapBins-1-iy] );
	}
	
	//classify into same direction and opposite direction cases
	for(int iy=0; iy<nDiffRapBins; iy++)
	{
		if( iy<nDiffRapBins/2 ) //first half, opposite direction
		{
			hMinY_OpstDir[iy]          = (TH1D*) hMinY_0nXnXn0n[iy]          -> Clone( Form("hMinY_OpstDir_iy%d",         iy) );
			hJpsiPtinY_OpstDir[iy]     = (TH1D*) hJpsiPtinY_0nXnXn0n[iy]     -> Clone( Form("hJpsiPtinY_OpstDir_iy%d",    iy) );
			hSideBandPtinY_OpstDir[iy] = (TH1D*) hSideBandPtinY_0nXnXn0n[iy] -> Clone( Form("hSideBandPtinY_OpstDir_iy%d",iy) );
		}
		else                   //second half, same direction 
		{
			int iy_idx1 = (int)(iy - nDiffRapBins/2.);
			int iy_idx2 = (int)(-iy + (3./2)*nDiffRapBins -1);

			//cout<<"iy: "<<iy<<" iy_idx1: "<<iy_idx1<<" iy_idx2: "<<iy_idx2<<endl;

			hMinY_SameDir[iy_idx1]          = (TH1D*) hMinY_0nXnXn0n[iy_idx2]          -> Clone( Form("hMinY_SameDir_iy%d",          iy_idx1) );
			hJpsiPtinY_SameDir[iy_idx1]     = (TH1D*) hJpsiPtinY_0nXnXn0n[iy_idx2]     -> Clone( Form("hJpsiPtinY_SameDir_iy%d",     iy_idx1) );
			hSideBandPtinY_SameDir[iy_idx1] = (TH1D*) hSideBandPtinY_0nXnXn0n[iy_idx2] -> Clone( Form("hSideBandPtinY_SameDir_iy%d", iy_idx1) );
		}
	}

	//combine sub y-bins for same directions, or for opposite directions
	for(int iy=0; iy<nDiffRapBins/2; iy++)
	{
		if(iy==0) 
		{
			hMAllY_SameDir          = (TH1D*) hMinY_SameDir[iy]         ->Clone("hMAllY_SameDir");
			hMAllY_OpstDir          = (TH1D*) hMinY_OpstDir[iy]         ->Clone("hMAllY_OpstDir");
			hJpsiPtAllY_SameDir     = (TH1D*) hJpsiPtinY_SameDir[iy]    ->Clone("hJpsiPtAllY_SameDir");
			hJpsiPtAllY_OpstDir     = (TH1D*) hJpsiPtinY_OpstDir[iy]    ->Clone("hJpsiPtAllY_OpstDir");
			hSideBandPtAllY_SameDir = (TH1D*) hSideBandPtinY_SameDir[iy]->Clone("hSideBandPtAllY_SameDir");
			hSideBandPtAllY_OpstDir = (TH1D*) hSideBandPtinY_OpstDir[iy]->Clone("hSideBandPtAllY_OpstDir");
		}
		else
		{
			hMAllY_SameDir          -> Add( hMinY_SameDir[iy]      );
			hMAllY_OpstDir          -> Add( hMinY_OpstDir[iy]      );
			hJpsiPtAllY_SameDir     -> Add( hJpsiPtinY_SameDir[iy] );
			hJpsiPtAllY_OpstDir     -> Add( hJpsiPtinY_OpstDir[iy] );
			hSideBandPtAllY_SameDir -> Add( hSideBandPtinY_SameDir[iy] );
			hSideBandPtAllY_OpstDir -> Add( hSideBandPtinY_OpstDir[iy] );
		}
	}

	//Rebin
	for(int iy=0; iy<nDiffRapBins; iy++)
	{
		hJpsiPtinY_0nXn[iy]                ->Rebin(2);
		hJpsiPtinY_Xn0n[nDiffRapBins-iy-1] ->Rebin(2);
	}

	for(int iy=0; iy<nDiffRapBins/2; iy++)
	{
		hMinY_SameDir[iy] ->Rebin(2);
		hMinY_OpstDir[iy] ->Rebin(2);
	}

	//Rebin for Physics
}

//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
void fitMass4QEDYield( )
{
	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
	setPad(0.12, 0.08, 0.07, 0.13);
	c1->SetLogy(0);
	
	const double massLow4Fit=2.6; const double massHig4Fit=3.40;
	RooRealVar  mMass("mMass", "m_{#mu#mu} (GeV)", massLow4Fit, massHig4Fit);

	TFile *inf_Temps = TFile::Open("../simulation/out4effAndTemp/MassPtTemp_AllSpecs_massWindow_2.95_3.25_.root");
	TF1 *fQED = new TF1("fQED", fReject4QED, massLow4Fit, massHig4Fit, 4);
	TF1 *fCohJpsiTemp;
	
	for(int idir=0; idir<2; idir++)
	{
		for(int iy=0; iy<nDiffRapBins/2; iy++)
		{
			
			TH1D* hMass_inWork;

			if(idir==0) //same direction
			{
				hMass_inWork = (TH1D*) hMinY_SameDir[iy] -> Clone("hMass_inWork");
			}
			else
			{

				hMass_inWork = (TH1D*) hMinY_OpstDir[iy] -> Clone("hMass_inWork");
			}

			//------------------------------------------------------------------------------------------------------------
			//1. Fit mass in full pt range to get the QED yield value
			//------------------------------------------------------------------------------------------------------------

			//------------------------------------------------------------------------------------------------------------
			fCohJpsiTemp = (TF1 *) inf_Temps->Get( "fCohJpsiTemp" );
			//------------------------------------------------------------------------------------------------------------

			//------------------------------------------------------------------------------------------------------------
			RooRealVar  cbAlpha(    "cbAlpha",     "cbAlpha",    fCohJpsiTemp->GetParameter(1), 0.0, 20);
			RooRealVar  cbN(        "cbN",         "cbN",        fCohJpsiTemp->GetParameter(2), 0.0, 20);
			RooRealVar  sigmaRatio( "sigmaRatio",  "sigmaRatio", fCohJpsiTemp->GetParameter(3), 0.1, 10);
			RooRealVar  jpsiMu(     "jpsiMu",      "jpsiMu",     3.096, 3.0, 3.2  );
			RooRealVar  gausN(      "gausN",       "gausN",      3.6,   0,   10   );
			RooRealVar  jpsiSigma(  "jpsiSigma",   "jpsiSigma",  0.045, 0,  0.15  );

			RooGenericPdf *jpsiPdf = new RooGenericPdf("jpsiPdf", "jpsiPdf", 
					"ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*sigmaRatio,jpsiMu) + gausN*TMath::Gaus(mMass, jpsiMu, jpsiSigma)", 
					RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, gausN));

			//------------------------------------------------------------------------------------------------------------
			hMass_inWork ->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters
			hMass_inWork ->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters
			hMass_inWork ->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters

			RooRealVar mP0(  "mP0", "mP0",  fQED->GetParameter(0), -1.e5, 1.e5);
			RooRealVar mP1(  "mP1", "mP1",  fQED->GetParameter(1), -1.e5, 1.e5);
			RooRealVar mP2(  "mP2", "mP2",  fQED->GetParameter(2), -1.e5, 1.e5);
			RooRealVar mP3(  "mP3", "mP3",  fQED->GetParameter(3), -1.e5, 1.e5);
			RooGenericPdf *qedPdf = new RooGenericPdf("qedPdf", "qedPdf", "mP0 + mP1*mMass + mP2*mMass*mMass + mP3*mMass*mMass*mMass", RooArgSet(mP0, mP1, mP2, mP3, mMass));
			//------------------------------------------------------------------------------------------------------------

			//------------------------------------------------------------------------------------------------------------
			const int    tem_JpsiBinLow = HistWorker::FindBin(hMass_inWork, 2.95, 0);
			const int    tem_JpsiBinHig = HistWorker::FindBin(hMass_inWork, 3.25, 1);
			const double nJpsi4Init     = hMass_inWork->Integral(tem_JpsiBinLow, tem_JpsiBinHig);
			const int    tem_QEDBinLow  = HistWorker::FindBin(hMass_inWork, 2.60, 0);
			const int    tem_QEDBinHig  = HistWorker::FindBin(hMass_inWork, 2.90, 1);
			const double nQED4Init      = hMass_inWork->Integral(tem_QEDBinLow, tem_QEDBinHig)*(massHig4Fit-massLow4Fit)/(2.90-2.60);

			RooRealVar nJpsi("nJpsi", "nJpsi", nJpsi4Init*0.90,  0, 1.e6);
			RooRealVar nQED( "nQED",  "nQED",  nQED4Init*0.95,   0, 1.e6);
			RooAddPdf  totMassPdf("totMassPdf", "totMassPdf", RooArgList(*jpsiPdf, *qedPdf), RooArgList(nJpsi, nQED)); 

			RooDataHist dataMass("dataMass", "dataMass", mMass, hMass_inWork); 

			//------------------------------------------------------------------------------------------------------------
			//------------------------------------------------------------------------------------------------------------
			totMassPdf.fitTo( dataMass, Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE), Save());
			RooFitResult *ResFit = totMassPdf.fitTo( dataMass, Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE), Save());
			//------------------------------------------------------------------------------------------------------------
			//------------------------------------------------------------------------------------------------------------
			
			const double nJpsiValue  = nJpsi.getVal();
			const double nJpsiError  = nJpsi.getError();
					
			TF1* fQED4frac = new TF1("fQED4frac", "[0] + [1]*x + [2]*x*x +[3]*x*x*x", 0, 5);
			fQED4frac ->SetParameters( mP0.getVal(), mP1.getVal(), mP2.getVal(), mP3.getVal() );
			
			const double  fracQED  = fQED4frac->Integral(mJpsiMassLow, mJpsiMassHi) / fQED4frac->Integral(massLow4Fit, massHig4Fit);
						
			const double nQEDinJpsi     = nQED.getVal()  *fracQED;
			const double nQEDinJpsiErr  = nQED.getError()*fracQED;
			
			nQEDinJpsi_iDir_iY[idir][iy] = nQEDinJpsi;

			c1->cd();
			c1->SetLogy(0);

			int nFrameMBins  = (massHig4Fit - massLow4Fit)/hMass_inWork->GetBinWidth(1);
			RooPlot *frameMass = mMass.frame(Range(massLow4Fit, massHig4Fit), Title(""), Bins(nFrameMBins));
			frameMass ->SetMaximum( hMass_inWork->GetMaximum()*1.2 );
			frameMass ->GetYaxis()->SetTitleOffset(0.90);
			dataMass  .plotOn(frameMass, MarkerStyle(20), MarkerSize(1), MarkerColor(1), LineColor(1), LineWidth(2), DrawOption("pz"));
			totMassPdf.plotOn(frameMass, LineColor(2), LineStyle(1), LineWidth(2));
			totMassPdf.plotOn(frameMass, Components(RooArgSet(*jpsiPdf)), LineColor(kBlue),    LineStyle(5), LineWidth(2));
			totMassPdf.plotOn(frameMass, Components(RooArgSet(*qedPdf)),  LineColor(qedColor), LineStyle(2), LineWidth(3));

			//			cout<<endl;
			//			cout<<"******** Print frame ********"<<endl;
			//			frameMass->Print();
			//			cout<<"******** End ********"<<endl;
			//			cout<<endl;

			double chi2ndf = frameMass->chiSquare("totMassPdf_Norm[mMass]", "h_dataMass", 12);   //Jpsi+fit

			frameMass->Draw() ;
						
			TString yName = Form("%.1f<y<%.1f", fabs(mDiffRapHi[iy]), fabs(mDiffRapLow[iy]) );
			
			drawLatex(0.15, 0.86, "0nXn+Xn0n: J/#psi-Xn in "+DirNames[idir], mTextFont, 0.06, mTextColor);
			drawLatex(0.18, 0.80, yName,                mTextFont, 0.05, mTextColor);
			
			drawLatex(0.70, 0.75, Form("#chi^{2}/ndf = %1.1f", chi2ndf),                               mTextFont, mTextSize, 2);
			
			const double textDy = 0.05;
			drawLatex(0.25, 0.46+textDy, Form("N_{J/#psi} = %d #pm %d",   TMath::Nint(nJpsi.getVal()), TMath::Nint(nJpsi.getError())), mTextFont, mTextSize, mTextColor);
			drawLatex(0.25, 0.38+textDy, Form("N_{QED} = %d #pm %d",      TMath::Nint(nQED.getVal()),  TMath::Nint(nQED.getError()) ), mTextFont, mTextSize, mTextColor);
			drawLatex(0.25, 0.30+textDy, Form("N^{in J/#psi}_{QED} = %d #pm %d", (int)nQEDinJpsi, (int)nQEDinJpsiErr ), mTextFont, mTextSize, kGreen+2);

			//----------------------------------------------------------------------------------------------------------------------------------------------------
			c1->SaveAs( Form("out0nXn_4PtRatio/massSpec_0nXnXn0n_4RemoveQED_4PtRatio_idir%d_iy%d.png",  idir, iy) );
			c1->SaveAs( Form("out0nXn_4PtRatio/massSpec_0nXnXn0n_4RemoveQED_4PtRatio_idir%d_iy%d.pdf",  idir, iy) );
			//----------------------------------------------------------------------------------------------------------------------------------------------------

			delete ResFit;
			delete hMass_inWork;
			delete frameMass;
		}//iy
	}//idir

	for(int idir=0; idir<2; idir++)
	{
		//fit the Same Direction, all sub-ybin summed mass 
		TH1D* hMass_inWork;

		if(idir==0) //same direction
		{
			hMass_inWork = (TH1D*) hMAllY_SameDir -> Clone("hMass_inWork");
		}
		else
		{

			hMass_inWork = (TH1D*) hMAllY_OpstDir -> Clone("hMass_inWork");
		}

		//------------------------------------------------------------------------------------------------------------
		//1. Fit mass in full pt range to get the QED yield value
		//------------------------------------------------------------------------------------------------------------

		//------------------------------------------------------------------------------------------------------------
		fCohJpsiTemp = (TF1 *) inf_Temps->Get( "fCohJpsiTemp" );
		//------------------------------------------------------------------------------------------------------------

		//------------------------------------------------------------------------------------------------------------
		RooRealVar  cbAlpha(    "cbAlpha",     "cbAlpha",    fCohJpsiTemp->GetParameter(1), 0.0, 20);
		RooRealVar  cbN(        "cbN",         "cbN",        fCohJpsiTemp->GetParameter(2), 0.0, 20);
		RooRealVar  sigmaRatio( "sigmaRatio",  "sigmaRatio", fCohJpsiTemp->GetParameter(3), 0.1, 10);
		RooRealVar  jpsiMu(     "jpsiMu",      "jpsiMu",     3.096, 3.0, 3.2  );
		RooRealVar  gausN(      "gausN",       "gausN",      3.6,   0,   10   );
		RooRealVar  jpsiSigma(  "jpsiSigma",   "jpsiSigma",  0.045, 0,  0.15  );

		RooGenericPdf *jpsiPdf = new RooGenericPdf("jpsiPdf", "jpsiPdf", 
				"ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*sigmaRatio,jpsiMu) + gausN*TMath::Gaus(mMass, jpsiMu, jpsiSigma)", 
				RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, gausN));

		//------------------------------------------------------------------------------------------------------------
		hMass_inWork ->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters
		hMass_inWork ->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters
		hMass_inWork ->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters

		RooRealVar mP0(  "mP0", "mP0",  fQED->GetParameter(0), -1.e5, 1.e5);
		RooRealVar mP1(  "mP1", "mP1",  fQED->GetParameter(1), -1.e5, 1.e5);
		RooRealVar mP2(  "mP2", "mP2",  fQED->GetParameter(2), -1.e5, 1.e5);
		RooRealVar mP3(  "mP3", "mP3",  fQED->GetParameter(3), -1.e5, 1.e5);
		RooGenericPdf *qedPdf = new RooGenericPdf("qedPdf", "qedPdf", "mP0 + mP1*mMass + mP2*mMass*mMass + mP3*mMass*mMass*mMass", RooArgSet(mP0, mP1, mP2, mP3, mMass));
		//------------------------------------------------------------------------------------------------------------

		//------------------------------------------------------------------------------------------------------------
		const int    tem_JpsiBinLow = HistWorker::FindBin(hMass_inWork, 2.95, 0);
		const int    tem_JpsiBinHig = HistWorker::FindBin(hMass_inWork, 3.25, 1);
		const double nJpsi4Init     = hMass_inWork->Integral(tem_JpsiBinLow, tem_JpsiBinHig);
		const int    tem_QEDBinLow  = HistWorker::FindBin(hMass_inWork, 2.60, 0);
		const int    tem_QEDBinHig  = HistWorker::FindBin(hMass_inWork, 2.90, 1);
		const double nQED4Init      = hMass_inWork->Integral(tem_QEDBinLow, tem_QEDBinHig)*(massHig4Fit-massLow4Fit)/(2.90-2.60);

		RooRealVar nJpsi("nJpsi", "nJpsi", nJpsi4Init*0.90,  0, 1.e6);
		RooRealVar nQED( "nQED",  "nQED",  nQED4Init*0.95,   0, 1.e6);
		RooAddPdf  totMassPdf("totMassPdf", "totMassPdf", RooArgList(*jpsiPdf, *qedPdf), RooArgList(nJpsi, nQED)); 

		RooDataHist dataMass("dataMass", "dataMass", mMass, hMass_inWork); 

		//------------------------------------------------------------------------------------------------------------
		//------------------------------------------------------------------------------------------------------------
		totMassPdf.fitTo( dataMass, Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE), Save());
		RooFitResult *ResFit = totMassPdf.fitTo( dataMass, Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE), Save());
		//------------------------------------------------------------------------------------------------------------
		//------------------------------------------------------------------------------------------------------------

		const double nJpsiValue  = nJpsi.getVal();
		const double nJpsiError  = nJpsi.getError();

		TF1* fQED4frac = new TF1("fQED4frac", "[0] + [1]*x + [2]*x*x +[3]*x*x*x", 0, 5);
		fQED4frac ->SetParameters( mP0.getVal(), mP1.getVal(), mP2.getVal(), mP3.getVal() );

		const double  fracQED  = fQED4frac->Integral(mJpsiMassLow, mJpsiMassHi) / fQED4frac->Integral(massLow4Fit, massHig4Fit);

		const double nQEDinJpsi     = nQED.getVal()  *fracQED;
		const double nQEDinJpsiErr  = nQED.getError()*fracQED;

		if(idir==0) nQEDinJpsi_SameDir_AllY = nQEDinJpsi;
		else        nQEDinJpsi_OpstDir_AllY = nQEDinJpsi;
		
		c1->cd();
		c1->SetLogy(0);

		int nFrameMBins  = (massHig4Fit - massLow4Fit)/hMass_inWork->GetBinWidth(1);
		RooPlot *frameMass = mMass.frame(Range(massLow4Fit, massHig4Fit), Title(""), Bins(nFrameMBins));
		frameMass ->SetMaximum( hMass_inWork->GetMaximum()*1.2 );
		frameMass ->GetYaxis()->SetTitleOffset(0.90);
		dataMass  .plotOn(frameMass, MarkerStyle(20), MarkerSize(1), MarkerColor(1), LineColor(1), LineWidth(2), DrawOption("pz"));
		totMassPdf.plotOn(frameMass, LineColor(2), LineStyle(1), LineWidth(2));
		totMassPdf.plotOn(frameMass, Components(RooArgSet(*jpsiPdf)), LineColor(kBlue),    LineStyle(5), LineWidth(2));
		totMassPdf.plotOn(frameMass, Components(RooArgSet(*qedPdf)),  LineColor(qedColor), LineStyle(2), LineWidth(3));

		//			cout<<endl;
		//			cout<<"******** Print frame ********"<<endl;
		//			frameMass->Print();
		//			cout<<"******** End ********"<<endl;
		//			cout<<endl;

		double chi2ndf = frameMass->chiSquare("totMassPdf_Norm[mMass]", "h_dataMass", 12);   //Jpsi+fit

		frameMass->Draw() ;

		TString yName = Form("%.1f<y<%.1f", mDiffRapLow[nDiffRapBins/2],mDiffRapHi[nDiffRapBins-1] );

		drawLatex(0.15, 0.86, "0nXn+Xn0n: J/#psi-Xn in "+DirNames[idir], mTextFont, 0.06, mTextColor);
		drawLatex(0.18, 0.80, yName,                                     mTextFont, 0.05, mTextColor);

		drawLatex(0.70, 0.75, Form("#chi^{2}/ndf = %1.1f", chi2ndf),                               mTextFont, mTextSize, 2);

		const double textDy = 0.05;
		drawLatex(0.25, 0.46+textDy, Form("N_{J/#psi} = %d #pm %d",   TMath::Nint(nJpsi.getVal()), TMath::Nint(nJpsi.getError())), mTextFont, mTextSize, mTextColor);
		drawLatex(0.25, 0.38+textDy, Form("N_{QED} = %d #pm %d",      TMath::Nint(nQED.getVal()),  TMath::Nint(nQED.getError()) ), mTextFont, mTextSize, mTextColor);
		drawLatex(0.25, 0.30+textDy, Form("N^{in J/#psi}_{QED} = %d #pm %d", (int)nQEDinJpsi, (int)nQEDinJpsiErr ), mTextFont, mTextSize, kGreen+2);

		//----------------------------------------------------------------------------------------------------------------------------------------------------
		if(idir==0)
		{
			c1->SaveAs( "out0nXn_4PtRatio/massSpec_0nXnXn0n_4RemoveQED_4PtRatio_SameDir_AllY.png" );
			c1->SaveAs( "out0nXn_4PtRatio/massSpec_0nXnXn0n_4RemoveQED_4PtRatio_SameDir_AllY.pdf" );
		}
		else
		{
			c1->SaveAs( "out0nXn_4PtRatio/massSpec_0nXnXn0n_4RemoveQED_4PtRatio_OpstDir_AllY.png" );
			c1->SaveAs( "out0nXn_4PtRatio/massSpec_0nXnXn0n_4RemoveQED_4PtRatio_OpstDir_AllY.pdf" );
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------

	}//idir
}
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
void compPtShape( )
{

	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
	setPad(0.12, 0.08, 0.07, 0.13);
	c1->SetLogy(1);
	
	//compare 0nXn-ybin vs Xn0n-ybin-mirror
	for(int iy=0; iy<nDiffRapBins; iy++)
	{
		c1->cd();

		hJpsiPtinY_0nXn[iy]                 -> SetTitle("");
		hJpsiPtinY_0nXn[iy]                 -> SetYTitle("Entries");
		hJpsiPtinY_0nXn[iy]                 -> SetXTitle("p_{T} GeV/c");
		hJpsiPtinY_0nXn[iy]                 -> SetMarkerStyle(20);
		hJpsiPtinY_0nXn[iy]                 -> SetMarkerSize(0.6);
		hJpsiPtinY_0nXn[iy]                 -> SetMarkerColor(1);
		hJpsiPtinY_0nXn[iy]                 -> SetLineColor(1);
		hJpsiPtinY_0nXn[iy]                 -> SetAxisRange(-0.05, 3.0, "x");
		hJpsiPtinY_0nXn[iy]                 -> SetAxisRange(0.5, 2.e3, "y");
		hJpsiPtinY_0nXn[iy]                 -> Draw("pe");
		hJpsiPtinY_Xn0n[nDiffRapBins-iy-1]  -> SetMarkerStyle(24);
		hJpsiPtinY_Xn0n[nDiffRapBins-iy-1]  -> SetMarkerSize(0.6);
		hJpsiPtinY_Xn0n[nDiffRapBins-iy-1]  -> SetMarkerColor(4);
		hJpsiPtinY_Xn0n[nDiffRapBins-iy-1]  -> SetLineColor(4);
		hJpsiPtinY_Xn0n[nDiffRapBins-iy-1]  -> Draw("pesame");

		TLegend* leg0 = new TLegend(0.30, 0.70, 0.88, 0.89);
		leg0->SetBorderSize(0);
		leg0->SetFillColor(0);
		leg0->SetTextSize(0.045);
		leg0->AddEntry( hJpsiPtinY_0nXn[iy], Form("0nXn: %.2f<M<%.2f, %.1f<y< %.1f", mJpsiMassLow,mJpsiMassHi,mDiffRapLow[iy],               mDiffRapHi[iy]               ), "lp");
		leg0->AddEntry( hJpsiPtinY_Xn0n[iy], Form("Xn0n: %.2f<M<%.2f, %.1f<y< %.1f", mJpsiMassLow,mJpsiMassHi,mDiffRapLow[nDiffRapBins-iy-1],mDiffRapHi[nDiffRapBins-iy-1]), "lp");
		leg0->Draw();

		c1->SaveAs( Form("out0nXn_4PtRatio/hPt_inYBin_0nXn_vs_Xn0n_iy%d.png", iy) );
		c1->SaveAs( Form("out0nXn_4PtRatio/hPt_inYBin_0nXn_vs_Xn0n_iy%d.pdf", iy) );
	}//iy


	//compare the same-direction and opposite direction pt distribution for the combined 0nXnXn0n data
	for(int iy=0; iy<nDiffRapBins/2; iy++)
	{
		c1->cd();
		hJpsiPtinY_SameDir[iy]  -> SetTitle("");
		hJpsiPtinY_SameDir[iy]  -> SetYTitle("Entries");
		hJpsiPtinY_SameDir[iy]  -> SetXTitle("p_{T}^{#mu#mu} GeV/c");
		hJpsiPtinY_SameDir[iy]  -> SetMarkerStyle(20);
		hJpsiPtinY_SameDir[iy]  -> SetMarkerSize(0.6);
		hJpsiPtinY_SameDir[iy]  -> SetMarkerColor(1);
		hJpsiPtinY_SameDir[iy]  -> SetLineColor(1);
		hJpsiPtinY_SameDir[iy]  -> SetAxisRange(-0.05, 3.0, "x");
		hJpsiPtinY_SameDir[iy]  -> SetAxisRange(0.5, 2.e3, "y");
		hJpsiPtinY_SameDir[iy]  ->GetYaxis()->SetTitleSize(0.07);
		hJpsiPtinY_SameDir[iy]  ->GetYaxis()->SetTitleOffset(0.70);
		hJpsiPtinY_SameDir[iy]  ->GetYaxis()->SetLabelSize(0.05);
		hJpsiPtinY_SameDir[iy]  ->GetXaxis()->SetTitleSize(0.05);
		hJpsiPtinY_SameDir[iy]  ->GetXaxis()->SetTitleOffset(0.95);
		hJpsiPtinY_SameDir[iy]  ->GetXaxis()->SetLabelSize(0.03);
		hJpsiPtinY_SameDir[iy]  ->GetXaxis()->SetLabelFont(40);
		hJpsiPtinY_SameDir[iy]  -> Draw("pe");

		hJpsiPtinY_OpstDir[iy]  -> SetMarkerStyle(24);
		hJpsiPtinY_OpstDir[iy]  -> SetMarkerSize(0.6);
		hJpsiPtinY_OpstDir[iy]  -> SetMarkerColor(4);
		hJpsiPtinY_OpstDir[iy]  -> SetLineColor(4);
		hJpsiPtinY_OpstDir[iy]  -> Draw("pesame");
		
		TLegend* leg0 = new TLegend(0.20, 0.70, 0.88, 0.89);
		leg0->SetHeader( Form("0nXn+Xn0n: %.2f<M<%.2f GeV/c^{2}, %.1f<y<%.1f", mJpsiMassLow,mJpsiMassHi, fabs(mDiffRapHi[iy]), fabs(mDiffRapLow[iy])) );
		leg0->SetBorderSize(0);
		leg0->SetFillColor(0);
		leg0->SetTextSize(0.045);
		leg0->AddEntry( hJpsiPtinY_SameDir[iy], "J/#psi-Xn Same Direction",     "lp");
		leg0->AddEntry( hJpsiPtinY_OpstDir[iy], "J/#psi-Xn Opposite Direction", "lp");
		leg0->Draw();

		c1->SaveAs( Form("out0nXn_4PtRatio/hPt_inYBin_0nXnXn0n_SameOpstDirection_iy%d.png", iy) );
		c1->SaveAs( Form("out0nXn_4PtRatio/hPt_inYBin_0nXnXn0n_SameOpstDirection_iy%d.pdf", iy) );

		//---------------------------------------------------------------------------------------------------
		//compare inclusive Pt vs. QED and subtract QED
		//---------------------------------------------------------------------------------------------------
		hJpsiPtinY_SameDir[iy]      -> Draw("pe");
		hSideBandPtinY_SameDir[iy]  -> Scale(1./hSideBandPtinY_SameDir[iy]->Integral());
		hSideBandPtinY_SameDir[iy]  -> Scale( nQEDinJpsi_iDir_iY[0][iy] );
		hSideBandPtinY_SameDir[iy]  -> SetMarkerStyle(24);
		hSideBandPtinY_SameDir[iy]  -> SetMarkerSize(0.6);
		hSideBandPtinY_SameDir[iy]  -> SetMarkerColor(4);
		hSideBandPtinY_SameDir[iy]  -> SetLineColor(4);
		hSideBandPtinY_SameDir[iy]  -> Draw("pesame");
	
		hPureJpsiPtinY_SameDir[iy]  = (TH1D*) hJpsiPtinY_SameDir[iy]->Clone(Form("hPureJpsiPtinY_SameDir_iy%d",iy));
		hPureJpsiPtinY_SameDir[iy]  -> Add(hSideBandPtinY_SameDir[iy], -1);
		
		hPureJpsiPtinY_SameDir[iy]  -> SetMarkerStyle(20);
		hPureJpsiPtinY_SameDir[iy]  -> SetMarkerSize(0.6);
		hPureJpsiPtinY_SameDir[iy]  -> SetMarkerColor(2);
		hPureJpsiPtinY_SameDir[iy]  -> SetLineColor(2);
		hPureJpsiPtinY_SameDir[iy]  -> Draw("pesame");

		TLegend* leg1 = new TLegend(0.40, 0.60, 0.88, 0.80);
		leg1->SetHeader( Form("%.2f<M<%.2f GeV/c^{2}, %.1f<y<%.1f", mJpsiMassLow,mJpsiMassHi, fabs(mDiffRapHi[iy]), fabs(mDiffRapLow[iy])) );
		leg1->SetBorderSize(0);
		leg1->SetFillColor(0);
		leg1->SetTextSize(0.045);
		leg1->AddEntry( hJpsiPtinY_SameDir[iy],     "J/#psi + #gamma#gamma",                                             "lp");
		leg1->AddEntry( hSideBandPtinY_SameDir[iy], Form("#gamma#gamma fixed yield ( %.0f )", nQEDinJpsi_iDir_iY[0][iy]),   "lp");
		leg1->AddEntry( hPureJpsiPtinY_SameDir[iy], "J/#psi + #gamma#gamma - #gamma#gamma",                              "lp");
		leg1->Draw();

		drawLatex(0.15, 0.86, "0nXn+Xn0n: J/#psi-Xn in "+DirNames[0], mTextFont, 0.06, mTextColor);

		//---------------------------------------------------------------------------------------------------
		c1->SaveAs( Form("out0nXn_4PtRatio/hPt_inYBin_0nXnXn0n_SameDirection_comp2QEDandRm_iy%d.png", iy) );
		c1->SaveAs( Form("out0nXn_4PtRatio/hPt_inYBin_0nXnXn0n_SameDirection_comp2QEDandRm_iy%d.pdf", iy) );
		//---------------------------------------------------------------------------------------------------
		
		hJpsiPtinY_OpstDir[iy]  -> SetTitle("");
		hJpsiPtinY_OpstDir[iy]  -> SetYTitle("Entries");
		hJpsiPtinY_OpstDir[iy]  -> SetXTitle("p_{T}^{#mu#mu} GeV/c");
		hJpsiPtinY_OpstDir[iy]  -> SetMarkerStyle(20);
		hJpsiPtinY_OpstDir[iy]  -> SetMarkerSize(0.6);
		hJpsiPtinY_OpstDir[iy]  -> SetMarkerColor(1);
		hJpsiPtinY_OpstDir[iy]  -> SetLineColor(1);
		hJpsiPtinY_OpstDir[iy]  -> SetAxisRange(-0.05, 3.0, "x");
		hJpsiPtinY_OpstDir[iy]  -> SetAxisRange(0.5, 2.e3, "y");
		hJpsiPtinY_OpstDir[iy]  ->GetYaxis()->SetTitleSize(0.07);
		hJpsiPtinY_OpstDir[iy]  ->GetYaxis()->SetTitleOffset(0.70);
		hJpsiPtinY_OpstDir[iy]  ->GetYaxis()->SetLabelSize(0.05);
		hJpsiPtinY_OpstDir[iy]  ->GetXaxis()->SetTitleSize(0.05);
		hJpsiPtinY_OpstDir[iy]  ->GetXaxis()->SetTitleOffset(0.95);
		hJpsiPtinY_OpstDir[iy]  ->GetXaxis()->SetLabelSize(0.03);
		hJpsiPtinY_OpstDir[iy]  ->GetXaxis()->SetLabelFont(40);
		hJpsiPtinY_OpstDir[iy]  -> Draw("pe");

		hSideBandPtinY_OpstDir[iy]  -> Scale(1./hSideBandPtinY_OpstDir[iy]->Integral());
		hSideBandPtinY_OpstDir[iy]  -> Scale( nQEDinJpsi_iDir_iY[1][iy] );
		hSideBandPtinY_OpstDir[iy]  -> SetMarkerStyle(24);
		hSideBandPtinY_OpstDir[iy]  -> SetMarkerSize(0.6);
		hSideBandPtinY_OpstDir[iy]  -> SetMarkerColor(4);
		hSideBandPtinY_OpstDir[iy]  -> SetLineColor(4);
		hSideBandPtinY_OpstDir[iy]  -> Draw("pesame");
	
		hPureJpsiPtinY_OpstDir[iy]  = (TH1D*) hJpsiPtinY_OpstDir[iy]->Clone(Form("hPureJpsiPtinY_OpstDir_iy%d",iy));
		hPureJpsiPtinY_OpstDir[iy]  -> Add(hSideBandPtinY_OpstDir[iy], -1);

		hPureJpsiPtinY_OpstDir[iy]  -> SetMarkerStyle(20);
		hPureJpsiPtinY_OpstDir[iy]  -> SetMarkerSize(0.6);
		hPureJpsiPtinY_OpstDir[iy]  -> SetMarkerColor(2);
		hPureJpsiPtinY_OpstDir[iy]  -> SetLineColor(2);
		hPureJpsiPtinY_OpstDir[iy]  -> Draw("pesame");

		TLegend* leg2 = new TLegend(0.40, 0.60, 0.88, 0.80);
		leg2->SetHeader( Form("%.2f<M<%.2f GeV/c^{2}, %.1f<y<%.1f", mJpsiMassLow,mJpsiMassHi, fabs(mDiffRapHi[iy]), fabs(mDiffRapLow[iy])) );
		leg2->SetBorderSize(0);
		leg2->SetFillColor(0);
		leg2->SetTextSize(0.045);
		leg2->AddEntry( hJpsiPtinY_OpstDir[iy],     "J/#psi + #gamma#gamma",                                             "lp");
		leg2->AddEntry( hSideBandPtinY_OpstDir[iy], Form("#gamma#gamma fixed yield ( %.0f )", nQEDinJpsi_iDir_iY[1][iy]), "lp");
		leg2->AddEntry( hPureJpsiPtinY_OpstDir[iy], "J/#psi + #gamma#gamma - #gamma#gamma",                              "lp");
		leg2->Draw();

		drawLatex(0.15, 0.86, "0nXn+Xn0n: J/#psi-Xn in "+DirNames[1], mTextFont, 0.06, mTextColor);

		c1->SaveAs( Form("out0nXn_4PtRatio/hPt_inYBin_0nXnXn0n_OpstDirection_comp2QEDandRm_iy%d.png", iy) );
		c1->SaveAs( Form("out0nXn_4PtRatio/hPt_inYBin_0nXnXn0n_OpstDirection_comp2QEDandRm_iy%d.pdf", iy) );

	}//iy

	//plot all sub-y bins added, same direction vs opposite direction 
	c1->cd();
	hJpsiPtAllY_SameDir  -> SetTitle("");
	hJpsiPtAllY_SameDir  -> SetYTitle("Entries");
	hJpsiPtAllY_SameDir  -> SetXTitle("p_{T}^{#mu#mu} GeV/c");
	hJpsiPtAllY_SameDir  -> SetMarkerStyle(20);
	hJpsiPtAllY_SameDir  -> SetMarkerSize(0.6);
	hJpsiPtAllY_SameDir  -> SetMarkerColor(1);
	hJpsiPtAllY_SameDir  -> SetLineColor(1);
	hJpsiPtAllY_SameDir  -> SetAxisRange(-0.05, 3.0, "x");
	hJpsiPtAllY_SameDir  -> SetAxisRange(0.5, 2.e3, "y");
	hJpsiPtAllY_SameDir  ->GetYaxis()->SetTitleSize(0.07);
	hJpsiPtAllY_SameDir  ->GetYaxis()->SetTitleOffset(0.70);
	hJpsiPtAllY_SameDir  ->GetYaxis()->SetLabelSize(0.05);
	hJpsiPtAllY_SameDir  ->GetXaxis()->SetTitleSize(0.05);
	hJpsiPtAllY_SameDir  ->GetXaxis()->SetTitleOffset(0.95);
	hJpsiPtAllY_SameDir  ->GetXaxis()->SetLabelSize(0.03);
	hJpsiPtAllY_SameDir  ->GetXaxis()->SetLabelFont(40);

	hJpsiPtAllY_SameDir  -> Draw("pe");
	hJpsiPtAllY_OpstDir  -> SetMarkerStyle(24);
	hJpsiPtAllY_OpstDir  -> SetMarkerSize(0.6);
	hJpsiPtAllY_OpstDir  -> SetMarkerColor(4);
	hJpsiPtAllY_OpstDir  -> SetLineColor(4);
	hJpsiPtAllY_OpstDir  -> Draw("pesame");

	TLegend* leg0 = new TLegend(0.20, 0.70, 0.88, 0.89);
	leg0->SetHeader( Form("0nXn+Xn0n: %.2f<M<%.2f GeV/c^{2}, %.1f<y<%.1f", mJpsiMassLow,mJpsiMassHi, mDiffRapLow[nDiffRapBins/2],mDiffRapHi[nDiffRapBins-1]) );
	leg0->SetBorderSize(0);
	leg0->SetFillColor(0);
	leg0->SetTextSize(0.045);
	leg0->AddEntry( hJpsiPtAllY_SameDir, "J/#psi-Xn Same Direction",     "lp");
	leg0->AddEntry( hJpsiPtAllY_OpstDir, "J/#psi-Xn Opposite Direction", "lp");
	leg0->Draw();

	c1->SaveAs( "out0nXn_4PtRatio/hPt_AllY_0nXnXn0n_SameOpstDirection.png" );
	c1->SaveAs( "out0nXn_4PtRatio/hPt_AllY_0nXnXn0n_SameOpstDirection.pdf" );
	
	//compare inclusive Pt in Jpsi mass window to the QED contributions (SideBandShape*nQEDinJpsiM)
	//---------------------------------------------------------------------------------------------------
	//compare inclusive Pt vs. QED and subtract QED
	//---------------------------------------------------------------------------------------------------
	hJpsiPtAllY_SameDir      -> Draw("pe");
	hSideBandPtAllY_SameDir  -> Scale(1./hSideBandPtAllY_SameDir->Integral());
	hSideBandPtAllY_SameDir  -> Scale( nQEDinJpsi_SameDir_AllY );
	hSideBandPtAllY_SameDir  -> SetMarkerStyle(24);
	hSideBandPtAllY_SameDir  -> SetMarkerSize(0.6);
	hSideBandPtAllY_SameDir  -> SetMarkerColor(4);
	hSideBandPtAllY_SameDir  -> SetLineColor(4);
	hSideBandPtAllY_SameDir  -> Draw("pesame");

	hPureJpsiPtAllY_SameDir  = (TH1D*) hJpsiPtAllY_SameDir->Clone( "hPureJpsiPtinY_SameDir" );
	hPureJpsiPtAllY_SameDir  -> Add(hSideBandPtAllY_SameDir, -1);

	hPureJpsiPtAllY_SameDir  -> SetMarkerStyle(20);
	hPureJpsiPtAllY_SameDir  -> SetMarkerSize(0.6);
	hPureJpsiPtAllY_SameDir  -> SetMarkerColor(2);
	hPureJpsiPtAllY_SameDir  -> SetLineColor(2);
	hPureJpsiPtAllY_SameDir  -> Draw("pesame");

	TLegend* leg1 = new TLegend(0.40, 0.60, 0.88, 0.80);
	leg1->SetHeader( Form("%.2f<M<%.2f GeV/c^{2}, %.1f<y<%.1f", mJpsiMassLow,mJpsiMassHi, mDiffRapLow[nDiffRapBins/2],mDiffRapHi[nDiffRapBins-1]) );
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->SetTextSize(0.045);
	leg1->AddEntry( hJpsiPtAllY_SameDir,     "J/#psi + #gamma#gamma",                                             "lp");
	leg1->AddEntry( hSideBandPtAllY_SameDir, Form("#gamma#gamma fixed yield ( %.0f )", nQEDinJpsi_SameDir_AllY),   "lp");
	leg1->AddEntry( hPureJpsiPtAllY_SameDir, "J/#psi + #gamma#gamma - #gamma#gamma",                              "lp");
	leg1->Draw();

	drawLatex(0.15, 0.86, "0nXn+Xn0n: J/#psi-Xn in "+DirNames[0], mTextFont, 0.06, mTextColor);

	//---------------------------------------------------------------------------------------------------
	c1->SaveAs( Form("out0nXn_4PtRatio/hPt_AllY_0nXnXn0n_SameDirection_comp2QEDandRm.png") );
	c1->SaveAs( Form("out0nXn_4PtRatio/hPt_AllY_0nXnXn0n_SameDirection_comp2QEDandRm.pdf") );
	//---------------------------------------------------------------------------------------------------
	hJpsiPtAllY_OpstDir  -> SetTitle("");
	hJpsiPtAllY_OpstDir  -> SetYTitle("Entries");
	hJpsiPtAllY_OpstDir  -> SetXTitle("p_{T}^{#mu#mu} GeV/c");
	hJpsiPtAllY_OpstDir  -> SetMarkerStyle(20);
	hJpsiPtAllY_OpstDir  -> SetMarkerSize(0.6);
	hJpsiPtAllY_OpstDir  -> SetMarkerColor(1);
	hJpsiPtAllY_OpstDir  -> SetLineColor(1);
	hJpsiPtAllY_OpstDir  -> SetAxisRange(-0.05, 3.0, "x");
	hJpsiPtAllY_OpstDir  -> SetAxisRange(0.5, 2.e3, "y");
	hJpsiPtAllY_OpstDir  ->GetYaxis()->SetTitleSize(0.07);
	hJpsiPtAllY_OpstDir  ->GetYaxis()->SetTitleOffset(0.70);
	hJpsiPtAllY_OpstDir  ->GetYaxis()->SetLabelSize(0.05);
	hJpsiPtAllY_OpstDir  ->GetXaxis()->SetTitleSize(0.05);
	hJpsiPtAllY_OpstDir  ->GetXaxis()->SetTitleOffset(0.95);
	hJpsiPtAllY_OpstDir  ->GetXaxis()->SetLabelSize(0.03);
	hJpsiPtAllY_OpstDir  ->GetXaxis()->SetLabelFont(40);
	hJpsiPtAllY_OpstDir  -> Draw("pe");

	hSideBandPtAllY_OpstDir  -> Scale(1./hSideBandPtAllY_OpstDir->Integral());
	hSideBandPtAllY_OpstDir  -> Scale( nQEDinJpsi_OpstDir_AllY );
	hSideBandPtAllY_OpstDir  -> SetMarkerStyle(24);
	hSideBandPtAllY_OpstDir  -> SetMarkerSize(0.6);
	hSideBandPtAllY_OpstDir  -> SetMarkerColor(4);
	hSideBandPtAllY_OpstDir  -> SetLineColor(4);
	hSideBandPtAllY_OpstDir  -> Draw("pesame");

	hPureJpsiPtAllY_OpstDir  = (TH1D*) hJpsiPtAllY_OpstDir->Clone( "hPureJpsiPtinY_OpstDir" );
	hPureJpsiPtAllY_OpstDir  -> Add(hSideBandPtAllY_OpstDir, -1);

	hPureJpsiPtAllY_OpstDir  -> SetMarkerStyle(20);
	hPureJpsiPtAllY_OpstDir  -> SetMarkerSize(0.6);
	hPureJpsiPtAllY_OpstDir  -> SetMarkerColor(2);
	hPureJpsiPtAllY_OpstDir  -> SetLineColor(2);
	hPureJpsiPtAllY_OpstDir  -> Draw("pesame");

	TLegend* leg2 = new TLegend(0.40, 0.60, 0.88, 0.80);
	leg2->SetHeader( Form("%.2f<M<%.2f GeV/c^{2}, %.1f<y<%.1f", mJpsiMassLow,mJpsiMassHi, mDiffRapLow[nDiffRapBins/2],mDiffRapHi[nDiffRapBins-1]) );
	leg2->SetBorderSize(0);
	leg2->SetFillColor(0);
	leg2->SetTextSize(0.045);
	leg2->AddEntry( hJpsiPtAllY_OpstDir,     "J/#psi + #gamma#gamma",                                             "lp");
	leg2->AddEntry( hSideBandPtAllY_OpstDir, Form("#gamma#gamma fixed yield ( %.0f )", nQEDinJpsi_OpstDir_AllY),   "lp");
	leg2->AddEntry( hPureJpsiPtAllY_OpstDir, "J/#psi + #gamma#gamma - #gamma#gamma",                              "lp");
	leg2->Draw();

	drawLatex(0.15, 0.86, "0nXn+Xn0n: J/#psi-Xn in "+DirNames[1], mTextFont, 0.06, mTextColor);

	//---------------------------------------------------------------------------------------------------
	c1->SaveAs( Form("out0nXn_4PtRatio/hPt_AllY_0nXnXn0n_OpstDirection_comp2QEDandRm.png") );
	c1->SaveAs( Form("out0nXn_4PtRatio/hPt_AllY_0nXnXn0n_OpstDirection_comp2QEDandRm.pdf") );
	//---------------------------------------------------------------------------------------------------

	delete c1;
}
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void getPtRatios( ) //take ratio for the pure Jpsi pt shape
{
	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
	setPad(0.12, 0.08, 0.07, 0.13);
	c1->SetLogy(1);
	
	const double x1=-0.05; const double x2=2.9;
	const double y1=4.e-2; const double y2=2.2;
	//const double y1=3.e-2; const double y2=3.0;

	//TH1D* hRatioinY_OpstDir2SameDir[nDiffRapBins/2];
	//TH1D* hRatioAllY_OpstDir2SameDir;
	
	for(int iy=0; iy<nDiffRapBins/2; iy++)
	{
		hRatioinY_OpstDir2SameDir[iy] = (TH1D*) hPureJpsiPtinY_OpstDir[iy]  -> Clone( Form("hRatioinY_OpstDir2SameDir_iy%d",iy) );
		hRatioinY_OpstDir2SameDir[iy] -> Divide( hPureJpsiPtinY_SameDir[iy] );
		
		c1->cd();

		hRatioinY_OpstDir2SameDir[iy]  -> SetTitle("");
		hRatioinY_OpstDir2SameDir[iy]  -> SetYTitle("J/#psi-Xn direction: Opposite/Same");
		hRatioinY_OpstDir2SameDir[iy]  -> SetXTitle("p_{T}^{J/#psi} GeV/c");
		hRatioinY_OpstDir2SameDir[iy]  -> GetYaxis()->CenterTitle();
		hRatioinY_OpstDir2SameDir[iy]  -> GetYaxis()->SetTitleSize(0.06);
		hRatioinY_OpstDir2SameDir[iy]  -> GetYaxis()->SetTitleOffset(0.95);
		hRatioinY_OpstDir2SameDir[iy]  -> GetYaxis()->SetLabelSize(0.05);
		hRatioinY_OpstDir2SameDir[iy]  -> GetXaxis()->SetTitleSize(0.05);
		hRatioinY_OpstDir2SameDir[iy]  -> GetXaxis()->SetTitleOffset(0.95);
		hRatioinY_OpstDir2SameDir[iy]  -> GetXaxis()->SetLabelSize(0.04);
		hRatioinY_OpstDir2SameDir[iy]  -> GetXaxis()->SetLabelFont(40);
		hRatioinY_OpstDir2SameDir[iy]  -> SetMarkerStyle(20);
		hRatioinY_OpstDir2SameDir[iy]  -> SetMarkerSize(0.8);
		hRatioinY_OpstDir2SameDir[iy]  -> SetMarkerColor(1);
		hRatioinY_OpstDir2SameDir[iy]  -> SetLineColor(1);
		hRatioinY_OpstDir2SameDir[iy]  -> SetAxisRange( x1, x2, "x");
		hRatioinY_OpstDir2SameDir[iy]  -> SetAxisRange( y1, y2, "y");
		hRatioinY_OpstDir2SameDir[iy]  -> Draw("pe");

		drawLatex(0.15, 0.86, Form("0nXn+Xn0n, %.1f<y<%.1f", fabs(mDiffRapHi[iy]), fabs(mDiffRapLow[iy])), mTextFont, 0.06, mTextColor);

		c1->SaveAs( Form("out0nXn_4PtRatio/hPureJpsiPtRatio_0nXnXn0n_OpstDir2SameDir_iy%d.png",iy) );
		c1->SaveAs( Form("out0nXn_4PtRatio/hPureJpsiPtRatio_0nXnXn0n_OpstDir2SameDir_iy%d.pdf",iy) );
	}//iy

	c1->cd();
	hRatioAllY_OpstDir2SameDir  = (TH1D*) hPureJpsiPtAllY_OpstDir->Clone("hRatioAllY_OpstDir2SameDir");
	hRatioAllY_OpstDir2SameDir  -> Divide( hPureJpsiPtAllY_SameDir );

	hRatioAllY_OpstDir2SameDir  -> SetTitle("");
	hRatioAllY_OpstDir2SameDir  -> SetYTitle("J/#psi-Xn direction: Opposite/Same");
	hRatioAllY_OpstDir2SameDir  -> SetXTitle("p_{T}^{J/#psi} GeV/c");
	hRatioAllY_OpstDir2SameDir  -> GetYaxis()->CenterTitle();
	hRatioAllY_OpstDir2SameDir  -> GetYaxis()->SetTitleSize(0.06);
	hRatioAllY_OpstDir2SameDir  -> GetYaxis()->SetTitleOffset(0.95);
	hRatioAllY_OpstDir2SameDir  -> GetYaxis()->SetLabelSize(0.05);
	hRatioAllY_OpstDir2SameDir  -> GetXaxis()->SetTitleSize(0.05);
	hRatioAllY_OpstDir2SameDir  -> GetXaxis()->SetTitleOffset(0.95);
	hRatioAllY_OpstDir2SameDir  -> GetXaxis()->SetLabelSize(0.04);
	hRatioAllY_OpstDir2SameDir  -> GetXaxis()->SetLabelFont(40);
	hRatioAllY_OpstDir2SameDir  -> SetMarkerStyle(20);
	hRatioAllY_OpstDir2SameDir  -> SetMarkerSize(0.8);
	hRatioAllY_OpstDir2SameDir  -> SetMarkerColor(1);
	hRatioAllY_OpstDir2SameDir  -> SetLineColor(1);
	hRatioAllY_OpstDir2SameDir  -> SetAxisRange(x1, x2, "x");
	hRatioAllY_OpstDir2SameDir  -> SetAxisRange(y1, y2, "y");
	hRatioAllY_OpstDir2SameDir  -> Draw("pe");

	drawLatex(0.15, 0.86, Form("0nXn+Xn0n, %.1f<y<%.1f", mDiffRapLow[nDiffRapBins/2],mDiffRapHi[nDiffRapBins-1]), mTextFont, 0.06, mTextColor);

	c1->SaveAs( Form("out0nXn_4PtRatio/hPureJpsiPtRatio_AllY_0nXnXn0n_OpstDir2SameDir.png") );
	c1->SaveAs( Form("out0nXn_4PtRatio/hPureJpsiPtRatio_AllY_0nXnXn0n_OpstDir2SameDir.pdf") );

	//------------------------------------------------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------------------------------------------
	//Rebin and redo the ratio
	TH1D* hPureJpsiPtinY_OpstDir_Reb[nDiffRapBins/2];
	TH1D* hPureJpsiPtinY_SameDir_Reb[nDiffRapBins/2];

	//bins used for the first time CMS Forward meeting
	//	const int    nPtBins_4Ratio = 25;
	//	const double ptBds_4Ratio[nPtBins_4Ratio+1] =
	//	{
	//		0.0,  0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20,
	//		0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,  1.0,
	//		1.20, 1.40, 1.60, 1.80, 2.0,
	//		2.50, 3.0
	//	};
	
	const int    nPtBins_4Ratio = 25-5;
	const double ptBds_4Ratio[nPtBins_4Ratio+1] =
	{
		0.0, 0.04, 0.08, 0.12, 0.16, 0.20,
		0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,  1.0,
		1.20, 1.40, 1.60, 1.80, 2.0,
		2.50, 3.0
	};


	for(int iy=0; iy<nDiffRapBins/2; iy++)
	{
		hPureJpsiPtinY_OpstDir_Reb[iy] = (TH1D*) rebHisto(hPureJpsiPtinY_OpstDir[iy], Form("hPureJpsiPtinY_OpstDir_Reb_iy%d",iy), nPtBins_4Ratio, ptBds_4Ratio);
		hPureJpsiPtinY_SameDir_Reb[iy] = (TH1D*) rebHisto(hPureJpsiPtinY_SameDir[iy], Form("hPureJpsiPtinY_SameDir_Reb_iy%d",iy), nPtBins_4Ratio, ptBds_4Ratio);
	}

	TH1D* hPureJpsiPtAllY_OpstDir_Reb = (TH1D*) rebHisto( hPureJpsiPtAllY_OpstDir, "hPureJpsiPtAllY_OpstDir_Reb", nPtBins_4Ratio, ptBds_4Ratio);
	TH1D* hPureJpsiPtAllY_SameDir_Reb = (TH1D*) rebHisto( hPureJpsiPtAllY_SameDir, "hPureJpsiPtAllY_SameDir_Reb", nPtBins_4Ratio, ptBds_4Ratio);

	hRatioAllY_OpstDir2SameDir_Reb = (TH1D*) hPureJpsiPtAllY_OpstDir_Reb->Clone("hRatioAllY_OpstDir2SameDir_Reb");
	hRatioAllY_OpstDir2SameDir_Reb -> Divide( hPureJpsiPtAllY_SameDir_Reb );
	
	for(int iy=0; iy<nDiffRapBins/2; iy++)
	{
		hRatioinY_OpstDir2SameDir_Reb[iy] = (TH1D*) hPureJpsiPtinY_OpstDir_Reb[iy]  -> Clone( Form("hRatioinY_OpstDir2SameDir_Reb_iy%d",iy) );
		hRatioinY_OpstDir2SameDir_Reb[iy] -> Divide( hPureJpsiPtinY_SameDir_Reb[iy] );

		c1->cd();

		hRatioinY_OpstDir2SameDir_Reb[iy]  -> SetTitle("");
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> SetYTitle("J/#psi-Xn direction: Opposite/Same");
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> SetXTitle("p_{T}^{J/#psi} GeV/c");
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> GetYaxis()->CenterTitle();
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> GetYaxis()->SetTitleSize(0.06);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> GetYaxis()->SetTitleOffset(0.95);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> GetYaxis()->SetLabelSize(0.05);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> GetXaxis()->SetTitleSize(0.05);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> GetXaxis()->SetTitleOffset(0.99);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> GetXaxis()->SetLabelSize(0.04);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> GetXaxis()->SetLabelFont(40);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> SetMarkerStyle(20);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> SetMarkerSize(1.0);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> SetMarkerColor(1);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> SetLineColor(1);
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> SetAxisRange( x1, x2, "x");
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> SetAxisRange( y1, y2, "y");
		hRatioinY_OpstDir2SameDir_Reb[iy]  -> Draw("pe");

		drawLatex(0.15, 0.86, Form("0nXn+Xn0n, %.1f<y<%.1f", fabs(mDiffRapHi[iy]), fabs(mDiffRapLow[iy])) , mTextFont, 0.06, mTextColor);

		c1->SaveAs( Form("out0nXn_4PtRatio/hPureJpsiPtRatio_0nXnXn0n_OpstDir2SameDir_Reb_iy%d.png",iy) );
		c1->SaveAs( Form("out0nXn_4PtRatio/hPureJpsiPtRatio_0nXnXn0n_OpstDir2SameDir_Reb_iy%d.pdf",iy) );
	}

	//----------------------------------------------------------------------------------------------------
	//Overlay the two y-bins and compare the Ratio difference
	//----------------------------------------------------------------------------------------------------
	c1->cd();
	c1->SetLogy(0);
	hRatioinY_OpstDir2SameDir_Reb[0]  -> SetAxisRange( 0.00, 0.8, "x");
	hRatioinY_OpstDir2SameDir_Reb[0]  -> SetAxisRange( 0.05,  1.6, "y");
	//hRatioinY_OpstDir2SameDir_Reb[0]  -> SetAxisRange( 0.01, 1.60, "y");
	hRatioinY_OpstDir2SameDir_Reb[0]  -> Draw("pe");
	hRatioinY_OpstDir2SameDir_Reb[1]  -> SetMarkerColor(4);
	hRatioinY_OpstDir2SameDir_Reb[1]  -> SetLineColor(4);
	hRatioinY_OpstDir2SameDir_Reb[1]  -> Draw("pesame");
	
	TLegend* leg_comY = new TLegend(0.30, 0.70, 0.88, 0.89);
	leg_comY->SetBorderSize(0);
	leg_comY->SetFillColor(0);
	leg_comY->SetTextSize(0.045);
	leg_comY->AddEntry( hRatioinY_OpstDir2SameDir_Reb[0],  Form("0nXn+Xn0n, %.1f<y<%.1f", fabs(mDiffRapHi[0]), fabs(mDiffRapLow[0])));
	leg_comY->AddEntry( hRatioinY_OpstDir2SameDir_Reb[1],  Form("0nXn+Xn0n, %.1f<y<%.1f", fabs(mDiffRapHi[1]), fabs(mDiffRapLow[1])));
	leg_comY->Draw();

	//----------------------------------------------------------------------------------------------------
	c1->SaveAs( "out0nXn_4PtRatio/hPureJpsiPtRatio_0nXnXn0n_OpstDir2SameDir_Reb_comp_twoYBins.png" );
	c1->SaveAs( "out0nXn_4PtRatio/hPureJpsiPtRatio_0nXnXn0n_OpstDir2SameDir_Reb_comp_twoYBins.pdf" );
	//----------------------------------------------------------------------------------------------------

	//----------------------------------------------------------------------------------------------------
	c1->cd();
	c1->SetLogy(1);

	hRatioAllY_OpstDir2SameDir_Reb  -> SetTitle("");
	hRatioAllY_OpstDir2SameDir_Reb  -> SetYTitle("J/#psi-Xn direction: Opposite/Same");
	hRatioAllY_OpstDir2SameDir_Reb  -> SetXTitle("p_{T}^{J/#psi} GeV/c");
	hRatioAllY_OpstDir2SameDir_Reb  -> GetYaxis()->CenterTitle();
	hRatioAllY_OpstDir2SameDir_Reb  -> GetYaxis()->SetTitleSize(0.06);
	hRatioAllY_OpstDir2SameDir_Reb  -> GetYaxis()->SetTitleOffset(0.95);
	hRatioAllY_OpstDir2SameDir_Reb  -> GetYaxis()->SetLabelSize(0.05);
	hRatioAllY_OpstDir2SameDir_Reb  -> GetXaxis()->SetTitleSize(0.05);
	hRatioAllY_OpstDir2SameDir_Reb  -> GetXaxis()->SetTitleOffset(0.95);
	hRatioAllY_OpstDir2SameDir_Reb  -> GetXaxis()->SetLabelSize(0.04);
	hRatioAllY_OpstDir2SameDir_Reb  -> GetXaxis()->SetLabelFont(40);
	hRatioAllY_OpstDir2SameDir_Reb  -> SetMarkerStyle(20);
	hRatioAllY_OpstDir2SameDir_Reb  -> SetMarkerSize(1.0);
	hRatioAllY_OpstDir2SameDir_Reb  -> SetMarkerColor(1);
	hRatioAllY_OpstDir2SameDir_Reb  -> SetLineColor(1);
	hRatioAllY_OpstDir2SameDir_Reb  -> SetAxisRange(x1, x2, "x");
	hRatioAllY_OpstDir2SameDir_Reb  -> SetAxisRange(y1, y2, "y");
	hRatioAllY_OpstDir2SameDir_Reb  -> Draw("pe");

	drawLatex(0.15, 0.86, Form("0nXn+Xn0n, %.1f<y<%.1f", mDiffRapLow[nDiffRapBins/2],mDiffRapHi[nDiffRapBins-1]), mTextFont, 0.06, mTextColor);

	//----------------------------------------------------------------------------------------------------
	c1->SaveAs( Form("out0nXn_4PtRatio/hPureJpsiPtRatio_AllY_0nXnXn0n_OpstDir2SameDir_Reb.png") );
	c1->SaveAs( Form("out0nXn_4PtRatio/hPureJpsiPtRatio_AllY_0nXnXn0n_OpstDir2SameDir_Reb.pdf") );
	//----------------------------------------------------------------------------------------------------
}


