#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"

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

//------------------------------------------------------------------------------------------------------------
const Bool_t  mStorePDF = kFALSE;

const double mTinyNum = 1.e-6;
const double mOffSet  = 0.1;

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


const int NnCases = 6;
const TString nCasesName[NnCases] = {"AnAn", "OnOn", "0nXn", "Xn0n", "OnXnSum", "XnXn"};
TH1D* hCohMass_in_ny[NnCases][nDiffRapBins+1]; //last content is the sum of all y-bins, within coherent pt threshold
TH1D* hMass_in_ny[NnCases][nDiffRapBins+1];    //last content is the sum of all y-bins, for all pt 
TH1D* hPt_Jpsi_in_ny[NnCases][nDiffRapBins+1]; //last content is the sum of all y-bins
TH1D* H_EffVsY_CohJpsi;
TH1D* H_EffVsY_CohPsi;
TH1D* H_EffVsY_CohPsi2Jpsi;

double fD_inPtCut[NnCases][nDiffRapBins+1];    //fD within pt<0.20 GeV/c
double fDerr_inPtCut[NnCases][nDiffRapBins+1];    //fD within pt<0.20 GeV/c
double NJpsi_inMFit[NnCases][nDiffRapBins+1]; //# Jpsi within pt<0.20 GeV/c from mass fitting
double NerrJpsi_inMFit[NnCases][nDiffRapBins+1]; //# Jpsi within pt<0.20 GeV/c from mass fitting
double Eff_CohJpsi[NnCases][nDiffRapBins+1];     //efficiency of coherent jpsi
const double temAcc[nDiffRapBins+1] = {0.170, 0.348, 0.348, 0.170, 0.250};//need to be updated by corrected one later
//------------------------------------------------------------------------------------------------------------
void prepareData();
void loadEff();
void fitCohMass_4RNRfD( const double massLow4Fit=2.6, const double massHig4Fit=4.2);
void fitFullMassAndPt_4Decouple( const double massLow4Fit=2.6, const double massHig4Fit=4.2,  const double ptLow4Fit=0,     const double ptHig4Fit=3.5);
void saveFiles();
//------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------
void getJpsiPsi_nsns()
{
	prepareData();
	
	loadEff();
	
	fitCohMass_4RNRfD(2.6, 4.2);
	
	//fitFullMassAndPt_4Decouple(2.6,4.2, -0.01,3.0);
}
//------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------
void prepareData()
{
	TFile* infile = new TFile("../anaData/jpsiHistos/rawSig.root", "read");
	cout<<"readin: "<<infile->GetName()<<endl;

	TH3D *hMvsPtvsRap_NeuDir[nNeus][nNeus];
	TH3D *hMvsAsyPhivsRap_NeuDir[nNeus][nNeus];

	if( nNeus !=2 ) 
	{
		cout<<"we are not running the correct 0nXn!!!!!!!"<<endl;
		return;
	}

	for (int ip = 0; ip < nNeus; ip++) 
	{
		cout<<"iplus: "<<ip<<endl;

		for (int im = 0; im < nNeus; im++) 
		{
			cout<<"iminus: "<<im<<endl;
			
			hMvsPtvsRap_NeuDir[ip][im]     = (TH3D*)infile->Get( Form("hMvsPtvsRap_NeuDir%dp%dm",     ip, im) );
			
			cout<<"readin: "<<hMvsPtvsRap_NeuDir[ip][im]->GetName()     <<endl;
			cout<<"Entries: "<<hMvsPtvsRap_NeuDir[ip][im]->GetEntries() <<endl;
		}
	}
	
	TH3D *hMvsPtvsRap_AnAn    = (TH3D*) infile->Get("hMvsPtvsRap"); //All Sum, Inclusive
	TH3D *hMvsPtvsRap_0n0n    = (TH3D*) hMvsPtvsRap_NeuDir[0][0]->Clone( "hMvsPtvsRap_0n0n"    );
	TH3D *hMvsPtvsRap_XnXn    = (TH3D*) hMvsPtvsRap_NeuDir[1][1]->Clone( "hMvsPtvsRap_XnXn"    );

	TH3D *hMvsPtvsRap_0nXn    = (TH3D*) hMvsPtvsRap_NeuDir[1][0]->Clone( "hMvsPtvsRap_0nXn"    );
	TH3D *hMvsPtvsRap_Xn0n    = (TH3D*) hMvsPtvsRap_NeuDir[0][1]->Clone( "hMvsPtvsRap_Xn0n"    );
	
	TH3D *hMvsPtvsRap_0nXnSum = (TH3D*) hMvsPtvsRap_NeuDir[0][1]->Clone( "hMvsPtvsRap_0nXnSum" );
	hMvsPtvsRap_0nXnSum       -> Add(hMvsPtvsRap_NeuDir[1][0]);
	
	cout<<"hMvsPtvsRap_0nXn->GetEntries(): "<<hMvsPtvsRap_0nXn->GetEntries()<<endl;
	cout<<"hMvsPtvsRap_Xn0n->GetEntries(): "<<hMvsPtvsRap_Xn0n->GetEntries()<<endl;
	
	TH3D* hMvsPtvsRap_inWork;
	for(int i_ncase=0; i_ncase<NnCases; i_ncase++)
	{
		cout<<"i_ncase: "<<i_ncase<<endl;

		if(i_ncase==0) hMvsPtvsRap_inWork = (TH3D*) hMvsPtvsRap_AnAn    ->Clone( Form("hMvsPtvsRap_i_ncase%d",i_ncase) );
		if(i_ncase==1) hMvsPtvsRap_inWork = (TH3D*) hMvsPtvsRap_0n0n    ->Clone( Form("hMvsPtvsRap_i_ncase%d",i_ncase) );
		if(i_ncase==2) hMvsPtvsRap_inWork = (TH3D*) hMvsPtvsRap_0nXn    ->Clone( Form("hMvsPtvsRap_i_ncase%d",i_ncase) );
		if(i_ncase==3) hMvsPtvsRap_inWork = (TH3D*) hMvsPtvsRap_Xn0n    ->Clone( Form("hMvsPtvsRap_i_ncase%d",i_ncase) );
		if(i_ncase==4) hMvsPtvsRap_inWork = (TH3D*) hMvsPtvsRap_0nXnSum ->Clone( Form("hMvsPtvsRap_i_ncase%d",i_ncase) );
		if(i_ncase==5) hMvsPtvsRap_inWork = (TH3D*) hMvsPtvsRap_XnXn    ->Clone( Form("hMvsPtvsRap_i_ncase%d",i_ncase) );
		
		for(int iy=0; iy<nDiffRapBins; iy++)
		{
			cout<<"iy: "<<iy<<" "<<mDiffRapLow[iy]<<" <y< "<<mDiffRapHi[iy]<<endl;

			int rapBinLow   = hMvsPtvsRap_inWork->GetXaxis()->FindBin( mDiffRapLow[iy] + mTinyNum );
			int rapBinHi    = hMvsPtvsRap_inWork->GetXaxis()->FindBin( mDiffRapHi[iy]  - mTinyNum );
			
			int ptBinLow    = hMvsPtvsRap_inWork->GetYaxis()->FindBin( 0.00            + mTinyNum );
			int ptBinHi     = hMvsPtvsRap_inWork->GetYaxis()->FindBin( mPtCut4Coh      - mTinyNum ); //only look at pt<0.2 GeV/c for Coh signals, for fD
			int nptBinsMax  = hMvsPtvsRap_inWork->GetNbinsY();
			
			int mJpsiBinLow = hMvsPtvsRap_inWork->GetZaxis()->FindBin( mJpsiMassLow + mTinyNum );
			int mJpsiBinHi  = hMvsPtvsRap_inWork->GetZaxis()->FindBin( mJpsiMassHi  - mTinyNum );

			hCohMass_in_ny[i_ncase][iy] = (TH1D *)hMvsPtvsRap_inWork->ProjectionZ( Form("hCohMass_iNeuCase%d_iy%d", i_ncase, iy), rapBinLow, rapBinHi, ptBinLow,    ptBinHi    );
			
			hMass_in_ny[i_ncase][iy]    = (TH1D *)hMvsPtvsRap_inWork->ProjectionZ( Form("hMass_iNeuCase%d_iy%d",    i_ncase, iy), rapBinLow, rapBinHi, 1,           nptBinsMax );
			hPt_Jpsi_in_ny[i_ncase][iy] = (TH1D *)hMvsPtvsRap_inWork->ProjectionY( Form("hPt_Jpsi_iNeuCase%d_iy%d", i_ncase, iy), rapBinLow, rapBinHi, mJpsiBinLow, mJpsiBinHi );
			
			//			massBinLow = hMvsPtvsRap_inWork->GetZaxis()->FindBin(mLowMassBandLow + mTinyNum);
			//			massBinHi  = hMvsPtvsRap_inWork->GetZaxis()->FindBin(mLowMassBandHi  - mTinyNum);
			//			hLowMassBandPt_NeuDir[ip][im][iy] 
			//				= (TH1D *)hMvsPtvsRap_inWork->ProjectionY(Form("hLowMassBandPt_NeuDir%dp%dm_RapBin%d", ip, im, iy), rapBinLow, rapBinHi, massBinLow, massBinHi);
			//			hLowMassBandPt_NeuDir[ip][im][iy]->SetTitle(Form("%1.1f < y < %1.1f", mDiffRapLow[iy], mDiffRapHi[iy]));
			//
			//			massBinLow = hMvsPtvsRap_inWork->GetZaxis()->FindBin(mJpsiMassLow + mTinyNum);
			//			massBinHi  = hMvsPtvsRap_inWork->GetZaxis()->FindBin(mJpsiMassHi - mTinyNum);
			//			hJpsiPt_NeuDir[ip][im][iy] 
			//				= (TH1D *)hMvsPtvsRap_inWork->ProjectionY(Form("hJpsiPt_NeuDir%dp%dm_RapBin%d", ip, im, iy), rapBinLow, rapBinHi, massBinLow, massBinHi);
			//			hJpsiPt_NeuDir[ip][im][iy]->SetTitle(Form("%1.1f < y < %1.1f", mDiffRapLow[iy], mDiffRapHi[iy]));
			//
			//			massBinLow = hMvsPtvsRap_inWork->GetZaxis()->FindBin(mHiMassBandLow + mTinyNum);
			//			massBinHi  = hMvsPtvsRap_inWork->GetZaxis()->FindBin(mHiMassBandHi - mTinyNum);
			//			hHiMassBandPt_NeuDir[ip][im][iy] = (TH1D *)hMvsPtvsRap_inWork->ProjectionY(Form("hHiMassBandPt_NeuDir%dp%dm_RapBin%d", ip, im, iy), rapBinLow, rapBinHi, massBinLow, massBinHi);
			//			hHiMassBandPt_NeuDir[ip][im][iy]->SetTitle(Form("%1.1f < y < %1.1f", mDiffRapLow[iy], mDiffRapHi[iy]));
		
			if( iy==0 )
			{
				hCohMass_in_ny[i_ncase][nDiffRapBins]  = (TH1D *) hCohMass_in_ny[i_ncase][iy] ->Clone( Form("hCohMass_in_iNeuCase%d", i_ncase) );
				
				hMass_in_ny[i_ncase][nDiffRapBins]     = (TH1D *) hMass_in_ny[i_ncase][iy]    ->Clone( Form("hMass_in_iNeuCase%d",    i_ncase) );
				hPt_Jpsi_in_ny[i_ncase][nDiffRapBins]  = (TH1D *) hPt_Jpsi_in_ny[i_ncase][iy] ->Clone( Form("hPt_Jpsi_in_iNeuCase%d", i_ncase) );
				
				hCohMass_in_ny[i_ncase][nDiffRapBins]  -> SetTitle( Form("%1.1f < |y| < %1.1f", mDiffRapLow[nDiffRapBins/2], mDiffRapHi[nDiffRapBins-1]) );
				hMass_in_ny[i_ncase][nDiffRapBins]     -> SetTitle( Form("%1.1f < |y| < %1.1f", mDiffRapLow[nDiffRapBins/2], mDiffRapHi[nDiffRapBins-1]) );
				hPt_Jpsi_in_ny[i_ncase][nDiffRapBins]  -> SetTitle( Form("%1.1f < |y| < %1.1f", mDiffRapLow[nDiffRapBins/2], mDiffRapHi[nDiffRapBins-1]) );
			}
			else
			{
				hCohMass_in_ny[i_ncase][nDiffRapBins]  -> Add( hCohMass_in_ny[i_ncase][iy] );
				
				hMass_in_ny[i_ncase][nDiffRapBins]     -> Add( hMass_in_ny[i_ncase][iy] );
				hPt_Jpsi_in_ny[i_ncase][nDiffRapBins]  -> Add( hPt_Jpsi_in_ny[i_ncase][iy] );
			}
		}//iy
	}//i_ncase
}
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
void loadEff()
{
	TFile* infile_eff = new TFile("../simulation/out4effAndTemp/Efficiency_AllSpecs.root", "read");
	cout<<"readin: "<<infile_eff->GetName()<<endl;
	
	//use the same neutron configurations for all, due to the fact that the Jpsi reconstruction eff no much differences
	H_EffVsY_CohJpsi     = (TH1D*) infile_eff->Get("hEffvsRap_CohJpsi");        //hEffvsRap_CohJpsi_0n0n, 0nXn, XnXn
	H_EffVsY_CohPsi      = (TH1D*) infile_eff->Get("hEffvsRap_CohPsi2S");
	H_EffVsY_CohPsi2Jpsi = (TH1D*) infile_eff->Get("hEffvsRap_CohPsi2SFeeddown");
}
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
void fitCohMass_4RNRfD( const double massLow4Fit=2.6, const double massHig4Fit=4.2)
{
	cout<<"now, let's fit the pt spectras !!!"<<endl;
	
	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
	setPad(0.12, 0.08, 0.07, 0.13);
	c1->cd();
	c1->SetLogy(0);
	RooRealVar  mMass("mMass", "m_{#mu#mu} (GeV)", massLow4Fit, massHig4Fit);

	TFile *inf_Temps = TFile::Open("../simulation/out4effAndTemp/MassPtTemp_AllSpecs_massWindow_2.95_3.25_.root");
	TF1 *fQED;
	fQED = new TF1("fQED", fReject4QED, massLow4Fit, massHig4Fit, 4);
	TF1 *fCohJpsiTemp;
	//------------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------
	for(int i_ncase=0; i_ncase<NnCases; i_ncase++)
	{
		if(i_ncase != 5 ) continue; //tem skip
		
		cout<<"i_ncase: "<<i_ncase<<endl;

		for(int iy=0; iy<nDiffRapBins+1; iy++)
		{
			if(iy==nDiffRapBins) continue; //skip all y added fitting 
			//if(iy != 3)          continue;
			
			cout<<"iy: "<<iy<<" "<<mDiffRapLow[iy]<<" <y< "<<mDiffRapHi[iy]<<endl;
			
			//fCohJpsiTemp = (TF1  *) inf_Temps->Get( "fCohJpsiTemp" ); //if use total temp parms to initialize Jpsi sig parms
		
			if( i_ncase ==0 ) //AnAn
			{
				//if use temp parms from different y bins to initialize
				if(     iy==nDiffRapBins)fCohJpsiTemp = (TF1  *) inf_Temps->Get( "fCohJpsiTemp" ); 
				else if(iy<nDiffRapBins) fCohJpsiTemp = (TF1  *) inf_Temps->Get( Form("fCohJpsiTemp_RapBin%d", iy) );
			}
			else if(i_ncase>=2||i_ncase<=4)
			{
				if(     iy==nDiffRapBins)fCohJpsiTemp = (TF1  *) inf_Temps->Get( "fCohJpsi_0nXnTemp" ); 
				else if(iy<nDiffRapBins) fCohJpsiTemp = (TF1  *) inf_Temps->Get( Form("fCohJpsi_0nXnTemp_RapBin%d", iy) );
			}
			else
			{
				if(     iy==nDiffRapBins)fCohJpsiTemp = (TF1  *) inf_Temps->Get( "fCohJpsi_"+nCasesName[i_ncase]+"Temp" ); 
				else if(iy<nDiffRapBins) fCohJpsiTemp = (TF1  *) inf_Temps->Get( Form("fCohJpsi_"+nCasesName[i_ncase]+"Temp_RapBin%d", iy) );
			}
			
			double yMean = (mDiffRapLow[iy]+mDiffRapHi[iy])/2.;

			double eff_CohJpsi      = H_EffVsY_CohJpsi     ->GetBinContent(H_EffVsY_CohJpsi    ->FindBin(yMean));
			double eff_CohPsi       = H_EffVsY_CohPsi      ->GetBinContent(H_EffVsY_CohPsi     ->FindBin(yMean));
			double eff_CohPsi2Jpsi  = H_EffVsY_CohPsi2Jpsi ->GetBinContent(H_EffVsY_CohPsi2Jpsi->FindBin(yMean));
			
			Eff_CohJpsi[i_ncase][iy] = eff_CohJpsi;

			cout<<"eff_CohJpsi: "    <<eff_CohJpsi    <<endl;
			cout<<"eff_CohPsi: "     <<eff_CohPsi     <<endl;
			cout<<"eff_CohPsi2Jpsi: "<<eff_CohPsi2Jpsi<<endl;
			
			TH1D* hCohMass = (TH1D*) hCohMass_in_ny[i_ncase][iy]->Clone("hMass");
			
			if( i_ncase==0               ) hCohMass->Rebin(1);
			if( i_ncase==1 || i_ncase==2 ) hCohMass->Rebin(2);
			if( i_ncase==3 || i_ncase==4 ) hCohMass->Rebin(2);
			if( i_ncase==5               ) hCohMass->Rebin(3);
			//------------------------------------------------------------------------------------------------------------
			//------------------------------------------------------------------------------------------------------------
			
			const double Init_cbAlpha     = fCohJpsiTemp->GetParameter(1);
			const double Init_cbN         = fCohJpsiTemp->GetParameter(2);
			const double Init_sigmaRatio  = fCohJpsiTemp->GetParameter(3);

			RooRealVar  cbAlpha(    "cbAlpha",     "cbAlpha",    fCohJpsiTemp->GetParameter(1), Init_cbAlpha*0.0,    Init_cbAlpha*5.   );
			RooRealVar  cbN(        "cbN",         "cbN",        fCohJpsiTemp->GetParameter(2), Init_cbN*0.0,        Init_cbN*5.       );
			RooRealVar  sigmaRatio( "sigmaRatio",  "sigmaRatio", fCohJpsiTemp->GetParameter(3), Init_sigmaRatio*0.0, Init_sigmaRatio*5.);
			RooRealVar  jpsiMu(     "jpsiMu",      "jpsiMu",     3.096, 2.9,  3.3  );
			RooRealVar  gausN(      "gausN",       "gausN",      3.5,   0.00, 10   );
			RooRealVar  jpsiSigma(  "jpsiSigma",   "jpsiSigma",  0.045, 0,    0.15 );
			RooConstVar massRatio(  "massRatio",   "massRatio",  mPsi_PDG/mJpsi_PDG);

			RooGenericPdf *jpsiPdf = new RooGenericPdf("jpsiPdf", "jpsiPdf",
					"ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*sigmaRatio,jpsiMu) + gausN*TMath::Gaus(mMass, jpsiMu, jpsiSigma)", 
					RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, gausN));

			RooGenericPdf *psiPdf  = new RooGenericPdf("psiPdf",  "psiPdf",
					"ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*sigmaRatio*massRatio,jpsiMu*massRatio) + gausN*TMath::Gaus(mMass, jpsiMu*massRatio, jpsiSigma*massRatio)", 
					RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, massRatio, gausN)); // psiMu = jpsiMu * massRatio; psiSigma = jpsiSigma * massRatio

			
			hCohMass->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters
			hCohMass->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters

			RooConstVar mP0(  "mP0", "mP0",  fQED->GetParameter(0));
			RooConstVar mP1(  "mP1", "mP1",  fQED->GetParameter(1));
			RooConstVar mP2(  "mP2", "mP2",  fQED->GetParameter(2));
			RooConstVar mP3(  "mP3", "mP3",  fQED->GetParameter(3));
			//RooRealVar mP0(  "mP0", "mP0",  fQED->GetParameter(0), -1.e8,  1.e8);
			//RooRealVar mP1(  "mP1", "mP1",  fQED->GetParameter(1), 0.,     1.e5);
			//RooRealVar mP2(  "mP2", "mP2",  fQED->GetParameter(2), -1.e5,    0.);
			//RooRealVar mP3(  "mP3", "mP3",  fQED->GetParameter(3), 0.,    1.e4);
			RooGenericPdf *qedPdf = new RooGenericPdf("qedPdf", "qedPdf", "mP0 + mP1*mMass + mP2*mMass*mMass + mP3*mMass*mMass*mMass", RooArgSet(mP0, mP1, mP2, mP3, mMass));

			//// directly use QED template from simulation
			//int jpsiMBinLow = hQEDMassHistTemp->GetXaxis()->FindBin(massLow4Fit + mTinyNum);
			//int jpsiMBinHi  = hQEDMassHistTemp->GetXaxis()->FindBin(massHig4Fit - mTinyNum);
			//int jpsiMassBinLow = hQEDMassHistTemp->GetXaxis()->FindBin(mJpsiMassLow + mTinyNum);
			//int jpsiMassBinHi  = hQEDMassHistTemp->GetXaxis()->FindBin(mJpsiMassHi - mTinyNum);
			//double mQEDFrac = hQEDMassHistTemp->Integral(jpsiMassBinLow, jpsiMassBinHi)*1./hQEDMassHistTemp->Integral(jpsiMBinLow, jpsiMBinHi);

			//hQEDMassHistTemp->RebinX(5);
			//RooDataHist hQEDMassRooHist("hQEDMassRooHist", "hQEDMassRooHist", mMass, hQEDMassHistTemp);
			//RooHistPdf  qedPdf("qedPdf", "qedPdf", mMass, hQEDMassRooHist, 2); // RebinX and interpolation order to make the QED pdf smooth 
			//------------------------------------------------------------------------------------------------------------
			//------------------------------------------------------------------------------------------------------------
			
			const int    tem_QEDBinLow  = hCohMass->FindBin(3.30 + mTinyNum);
			const int    tem_QEDBinHig  = hCohMass->FindBin(3.50 - mTinyNum);
			const double nQED4Init      = hCohMass->Integral(tem_QEDBinLow, tem_QEDBinHig)*(massHig4Fit-massLow4Fit)/(3.50-3.30);

			const int    tem_JpsiBinLow = hCohMass->FindBin(2.80 + mTinyNum);
			const int    tem_JpsiBinHig = hCohMass->FindBin(3.30 - mTinyNum);
			//const double nJpsi4Init     = hCohMass->Integral(tem_JpsiBinLow, tem_JpsiBinHig)*0.80; //- nQED4Init*(3.30-2.80)/(massHig4Fit-massLow4Fit);
			const double nJpsi4Init     = hCohMass->Integral(tem_JpsiBinLow, tem_JpsiBinHig) - nQED4Init*(3.30-2.80)/(massHig4Fit-massLow4Fit);
			
			const double nPsi4Init      = nJpsi4Init*0.050;

			RooRealVar nJpsi("nJpsi", "nJpsi", nJpsi4Init*1.01,  nJpsi4Init*0.00, nJpsi4Init*10);
			RooRealVar nPsi( "nPsi",  "nPsi",  nPsi4Init*1.00,   nPsi4Init*0.00,  nPsi4Init*10 );
			RooRealVar nQED( "nQED",  "nQED",  nQED4Init,        nQED4Init*0.00,  nQED4Init*10 );

			RooAddPdf  totMassPdf("totMassPdf", "totMassPdf", RooArgList(*jpsiPdf, *psiPdf, *qedPdf), RooArgList(nJpsi, nPsi, nQED));

			RooDataHist dataMass("dataMass", "dataMass", mMass, hCohMass); 

			//------------------------------------------------------------------------------------------------------------
			//------------------------------------------------------------------------------------------------------------
			//totMassPdf.fitTo( dataMass, Range(2.70, 3.50), Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE));
			totMassPdf.fitTo( dataMass, Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE), Save());
			RooFitResult *ResFit = totMassPdf.fitTo( dataMass, Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE), Save());
			//totMassPdf.fitTo(dataMass,Range(massLow4Fit, massHig4Fit),Extended(kTRUE),SumW2Error(kTRUE),Hesse(kTRUE),Minos(kFALSE),Save());
			//------------------------------------------------------------------------------------------------------------

			//calculate ratio #Psi/#Jpsi and its uncertainty
			const TMatrixDSym &mtrx_cov = ResFit->covarianceMatrix();   //sigma_AB = rho_AB*sigmaA*sigmaB, rho_AB is the correlationCoefficient
			cout << "covariance matrix" << endl;
			mtrx_cov.Print();

			const double nJpsiValue  = nJpsi.getVal();
			const double nJpsiError  = nJpsi.getError();
			const double nPsiValue   = nPsi.getVal();
			const double nPsiError   = nPsi.getError();
			const double RN          = nPsiValue/nJpsiValue;
			const double RNErr       = RN*sqrt( pow(nJpsiError/nJpsiValue, 2) + pow(nPsiError/nPsiValue, 2) - (2.*mtrx_cov[7][8])/(nPsiValue*nJpsiValue) );

			cout<<"RN: "<<RN<<" +/- "<<RNErr<<endl;

			//------------------------------------------------------------------------------------------------------------
			//calculate R value and fD value
			//------------------------------------------------------------------------------------------------------------
			const double aa    = br_Jpsi2uu * eff_CohJpsi;
			const double bb    = br_Psi2uu  * eff_CohPsi;
			const double cc    = br_Psi2Jpsi*eff_CohPsi2Jpsi*br_Jpsi2uu;

			const double R     = ( RN * br_Jpsi2uu * eff_CohJpsi ) / ( br_Psi2uu*eff_CohPsi - RN*br_Psi2Jpsi*eff_CohPsi2Jpsi*br_Jpsi2uu  );
			const double fD    = R * (eff_CohPsi2Jpsi/eff_CohJpsi) * br_Psi2Jpsi;
			const double RErr  = ( aa/(bb-cc*RN) - (aa*cc*RN)/pow(bb-cc*RN,2) ) * RNErr;
			const double fDErr = (eff_CohPsi2Jpsi/eff_CohJpsi) * br_Psi2Jpsi * RErr;

			NJpsi_inMFit[i_ncase][iy]    = nJpsiValue;
			NerrJpsi_inMFit[i_ncase][iy] = nJpsiError;
			fD_inPtCut[i_ncase][iy]      = fD;
			fDerr_inPtCut[i_ncase][iy]   = fDErr;
			//------------------------------------------------------------------------------------------------------------

			c1->cd();
			c1->SetLogy(0);

			int nFrameMBins  = (massHig4Fit - massLow4Fit)/hCohMass->GetBinWidth(1);
			RooPlot *frameMass = mMass.frame(Range(massLow4Fit, massHig4Fit), Title(""), Bins(nFrameMBins));
			//frameMass ->GetYaxis()->SetTitleSize(0.10);
			frameMass ->GetYaxis()->SetTitleOffset(0.90);
			dataMass.plotOn(frameMass, MarkerStyle(20), MarkerSize(1), MarkerColor(1), LineColor(1), LineWidth(2), DrawOption("pz"));
			totMassPdf.plotOn(frameMass, LineColor(2), LineStyle(1), LineWidth(2));
			totMassPdf.plotOn(frameMass, Components(RooArgSet(*jpsiPdf)), LineColor(kBlue),    LineStyle(5), LineWidth(2));
			totMassPdf.plotOn(frameMass, Components(RooArgSet(*psiPdf)),  LineColor(kBlue+2),  LineStyle(6), LineWidth(2));
			totMassPdf.plotOn(frameMass, Components(RooArgSet(*qedPdf)),  LineColor(qedColor), LineStyle(2), LineWidth(3));

			//			cout<<endl;
			//			cout<<"******** Print frame ********"<<endl;
			//			frameMass->Print();
			//			cout<<"******** End ********"<<endl;
			//			cout<<endl;


			double chi2ndf = frameMass->chiSquare("totMassPdf_Norm[mMass]", "h_dataMass", 9);   //Jpsi+Psi fit

			frameMass->Draw() ;
			
			TString yName = "";
			if(iy<nDiffRapBins) yName = Form("%1.1f < y^{#mu#mu} < %1.1f",   mDiffRapLow[iy],             mDiffRapHi[iy]            );
			else                yName = Form("%1.1f < |y^{#mu#mu}| < %1.1f", mDiffRapLow[nDiffRapBins/2], mDiffRapHi[nDiffRapBins-1]);
			
			const TString ptName = Form("%1.0f < p_{T}^{#mu#mu} < %1.1f GeV/c",  0.0,   mPtCut4Coh);

			drawLatex(0.15, 0.86, nCasesName[i_ncase], mTextFont, 0.06, mTextColor);
			drawLatex(.55, 0.86, yName,                mTextFont, 0.05, mTextColor);
			drawLatex(.55, 0.80, ptName,               mTextFont, 0.05, mTextColor);
			drawLatex(.15, 0.75, Form("#chi^{2}/ndf = %1.1f", chi2ndf),                                 mTextFont, mTextSize, mTextColor);
			const double textDy = 0.05;
			drawLatex(.55, 0.66+textDy, Form("N_{J/#psi} = %d #pm %d",   TMath::Nint(nJpsi.getVal()), TMath::Nint(nJpsi.getError())), mTextFont, mTextSize, mTextColor);
			drawLatex(.55, 0.61+textDy, Form("N_{#psi(2S)} = %d #pm %d", TMath::Nint(nPsi.getVal()),  TMath::Nint(nPsi.getError()) ), mTextFont, mTextSize, mTextColor);
			drawLatex(.55, 0.55+textDy, Form("N_{QED} = %d #pm %d",      TMath::Nint(nQED.getVal()),  TMath::Nint(nQED.getError()) ), mTextFont, mTextSize, mTextColor);
			drawLatex(.45, 0.45+textDy, Form("R_{N} = #frac{N_{#psi(2S)}}{N_{J/#psi}} = %.3f #pm %.4f", nPsi.getVal()/nJpsi.getVal(), RNErr ), mTextFont, mTextSize, mTextColor);
			drawLatex(.45, 0.35+textDy, Form("R = #frac{#sigma_{#psi(2S)}}{#sigma_{J/#psi}} = %.3f #pm %.4f", R, RErr ), mTextFont, mTextSize, mTextColor);
			drawLatex(.45, 0.25+textDy, Form("f_{D} = #frac{FD J/#psi}{primary J/#psi} = %.3f #pm %.4f",      fD,fDErr),  mTextFont, mTextSize, mTextColor);

			//----------------------------------------------------------------------------------------------------------------------------------------------------
			c1->SaveAs( Form("outplots/massSpec_4JpsiPsi_"+nCasesName[i_ncase]+"_iy%d.png",  iy) );
			c1->SaveAs( Form("outplots/massSpec_4JpsiPsi_"+nCasesName[i_ncase]+"_iy%d.pdf",  iy) );
			//----------------------------------------------------------------------------------------------------------------------------------------------------
		
			delete ResFit;
			delete hCohMass;
			delete frameMass;
		}//iy
	}//incase

	delete c1;
	delete inf_Temps;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------
void fitFullMassAndPt_4Decouple( const double massLow4Fit=2.6, const double massHig4Fit=4.2, const double ptLow4Fit = 0.0, const double ptHig4Fit = 3.5 )
{
	cout<<"now, let's fit the pt spectras !!!"<<endl;
	TCanvas* c2 = new TCanvas("c2", "c2", 900, 900);
	
	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
	setPad(0.12, 0.08, 0.07, 0.13);
	c1->cd();
	c1->SetLogy(0);
	RooRealVar  mMass("mMass", "m_{#mu#mu} (GeV)", massLow4Fit, massHig4Fit);
	
	TFile *inf_Temps = TFile::Open("../simulation/out4effAndTemp/MassPtTemp_AllSpecs_massWindow_2.95_3.25_.root");
	TF1 *fQED = new TF1("fQED", fReject4QED, massLow4Fit, massHig4Fit, 4);
	TF1 *fCohJpsiTemp;
	
	for(int i_ncase=0; i_ncase<NnCases; i_ncase++)
	{
		cout<<"i_ncase: "<<i_ncase<<endl;
		if(i_ncase != 0 ) continue; //tem skip

		for(int iy=0; iy<nDiffRapBins+1; iy++)
		{
			if(iy==nDiffRapBins) continue; //temperory to skip all y added fitting 
			//if(iy != 0)          continue; //tem skip
			
			const double dY = mDiffRapHi[iy] - mDiffRapLow[iy]; //need to be updated if use for 1.6<|y|<2.4

			cout<<"iy: "<<iy<<" "<<mDiffRapLow[iy]<<" <y< "<<mDiffRapHi[iy]<<endl;
	
			TH1D* hMass = (TH1D*) hMass_in_ny[i_ncase][iy]    ->Clone("hMass");
			TH1D* hPt   = (TH1D*) hPt_Jpsi_in_ny[i_ncase][iy] ->Clone("hPt");

			if( i_ncase==0 )
			{
				hMass ->Rebin(2);
				hPt   ->Rebin(2);
			}
			else if ( i_ncase==2|| i_ncase==3 )
			{
				hMass ->Rebin(2);
				hPt   ->Rebin(2);
			}
			else if(  i_ncase==5 )
			{
				hMass ->Rebin(2);
				hPt   ->Rebin(2);
			}

			//------------------------------------------------------------------------------------------------------------
			//1. Fit mass in full pt range to get the QED yield value
			//------------------------------------------------------------------------------------------------------------
			
			//------------------------------------------------------------------------------------------------------------
			if(      iy==nDiffRapBins ) fCohJpsiTemp = (TF1 *) inf_Temps->Get( "fCohJpsiTemp" );
			else if( iy<nDiffRapBins  ) fCohJpsiTemp = (TF1 *) inf_Temps->Get( Form("fCohJpsiTemp_RapBin%d", iy) );
			//------------------------------------------------------------------------------------------------------------

			//------------------------------------------------------------------------------------------------------------
			RooRealVar  cbAlpha(    "cbAlpha",     "cbAlpha",    fCohJpsiTemp->GetParameter(1), 0,  20);
			RooRealVar  cbN(        "cbN",         "cbN",        fCohJpsiTemp->GetParameter(2), 0,  20);
			RooRealVar  sigmaRatio( "sigmaRatio",  "sigmaRatio", fCohJpsiTemp->GetParameter(3), 0.1,10);
			RooRealVar  jpsiMu(     "jpsiMu",      "jpsiMu",     3.096, 3.0, 3.2  );
			RooRealVar  gausN(      "gausN",       "gausN",      3.6,   0,   10   );
			RooRealVar  jpsiSigma(  "jpsiSigma",   "jpsiSigma",  0.045, 0,  0.15  );
			RooConstVar massRatio(  "massRatio",   "massRatio",  mPsi_PDG/mJpsi_PDG);

			RooGenericPdf *jpsiPdf = new RooGenericPdf("jpsiPdf", "jpsiPdf", 
					"ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*sigmaRatio,jpsiMu) + gausN*TMath::Gaus(mMass, jpsiMu, jpsiSigma)", 
					RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, gausN));

			RooGenericPdf *psiPdf  = new RooGenericPdf("psiPdf",  "psiPdf",  
					"ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*sigmaRatio*massRatio,jpsiMu*massRatio) + gausN*TMath::Gaus(mMass, jpsiMu*massRatio, jpsiSigma*massRatio)", 
					RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, massRatio, gausN)); // psiMu = jpsiMu * massRatio; psiSigma = jpsiSigma * massRatio

			//------------------------------------------------------------------------------------------------------------
			hMass ->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters
			hMass ->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters
			hMass ->Fit(fQED, "R", "",  massLow4Fit, massHig4Fit); //use side band to initialized QED parameters

			RooConstVar mP0(  "mP0", "mP0",  fQED->GetParameter(0));
			RooConstVar mP1(  "mP1", "mP1",  fQED->GetParameter(1));
			RooConstVar mP2(  "mP2", "mP2",  fQED->GetParameter(2));
			RooConstVar mP3(  "mP3", "mP3",  fQED->GetParameter(3));
			RooGenericPdf *qedPdf = new RooGenericPdf("qedPdf", "qedPdf", "mP0 + mP1*mMass + mP2*mMass*mMass + mP3*mMass*mMass*mMass", RooArgSet(mP0, mP1, mP2, mP3, mMass));
			//------------------------------------------------------------------------------------------------------------

			//------------------------------------------------------------------------------------------------------------
			
			const int    tem_JpsiBinLow = hMass->FindBin(2.95 + mTinyNum);
			const int    tem_JpsiBinHig = hMass->FindBin(3.25 - mTinyNum);
			const double nJpsi4Init     = hMass->Integral(tem_JpsiBinLow, tem_JpsiBinHig);
			const int    tem_QEDBinLow  = hMass->FindBin(3.25 + mTinyNum);
			const int    tem_QEDBinHig  = hMass->FindBin(3.50 - mTinyNum);
			const double nQED4Init      = hMass->Integral(tem_QEDBinLow, tem_QEDBinHig)*(massHig4Fit-massLow4Fit)/(3.50-3.25);

			RooRealVar nJpsi("nJpsi", "nJpsi", nJpsi4Init*0.80,  0, 1.e6);
			RooRealVar nPsi( "nPsi",  "nPsi",  nJpsi4Init*0.05,  0, 1.e4);
			RooRealVar nQED( "nQED",  "nQED",  nQED4Init*0.95,   0, 1.e6);
			RooAddPdf  totMassPdf("totMassPdf", "totMassPdf", RooArgList(*jpsiPdf, *psiPdf, *qedPdf), RooArgList(nJpsi, nPsi, nQED)); 

			RooDataHist dataMass("dataMass", "dataMass", mMass, hMass); 

			//------------------------------------------------------------------------------------------------------------
			//------------------------------------------------------------------------------------------------------------
			totMassPdf.fitTo( dataMass, Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE), Save());
			RooFitResult *ResFit = totMassPdf.fitTo( dataMass, Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE), Save());
			//------------------------------------------------------------------------------------------------------------
			//------------------------------------------------------------------------------------------------------------

			const double nJpsiValue  = nJpsi.getVal();
			const double nJpsiError  = nJpsi.getError();
			const double nPsiValue   = nPsi.getVal();
			const double nPsiError   = nPsi.getError();
		
			TF1* fQED4frac = new TF1("fQED4frac", "[0] + [1]*x + [2]*x*x +[3]*x*x*x", 0, 5);
			fQED4frac ->SetParameters( mP0.getVal(), mP1.getVal(), mP2.getVal(), mP3.getVal() );

			const double  fracQED  = fQED4frac->Integral(mJpsiMassLow, mJpsiMassHi) / fQED4frac->Integral(massLow4Fit, massHig4Fit);
			
			const double nQEDinJpsi     = nQED.getVal()  *fracQED;
			const double nQEDinJpsiErr  = nQED.getError()*fracQED;

			c1->cd();
			c1->SetLogy(0);

			int nFrameMBins  = (massHig4Fit - massLow4Fit)/hMass->GetBinWidth(1);
			RooPlot *frameMass = mMass.frame(Range(massLow4Fit, massHig4Fit), Title(""), Bins(nFrameMBins));
			//frameMass ->GetYaxis()->SetTitleSize(0.10);
			frameMass ->GetYaxis()->SetTitleOffset(0.90);
			dataMass  .plotOn(frameMass, MarkerStyle(20), MarkerSize(1), MarkerColor(1), LineColor(1), LineWidth(2), DrawOption("pz"));
			totMassPdf.plotOn(frameMass, LineColor(2), LineStyle(1), LineWidth(2));
			totMassPdf.plotOn(frameMass, Components(RooArgSet(*jpsiPdf)), LineColor(kBlue),    LineStyle(5), LineWidth(2));
			totMassPdf.plotOn(frameMass, Components(RooArgSet(*psiPdf)),  LineColor(kBlue+2),  LineStyle(6), LineWidth(2));
			totMassPdf.plotOn(frameMass, Components(RooArgSet(*qedPdf)),  LineColor(qedColor), LineStyle(2), LineWidth(3));

			//			cout<<endl;
			//			cout<<"******** Print frame ********"<<endl;
			//			frameMass->Print();
			//			cout<<"******** End ********"<<endl;
			//			cout<<endl;

			double chi2ndf = frameMass->chiSquare("totMassPdf_Norm[mMass]", "h_dataMass", 9);   //Jpsi+Psi fit

			frameMass->Draw() ;
			
			TString yName = "";
			if(iy<nDiffRapBins) yName = Form( "%1.1f < y^{#mu#mu} < %1.1f",   mDiffRapLow[iy],             mDiffRapHi[iy]             );
			else                yName = Form( "%1.1f < |y^{#mu#mu}| < %1.1f", mDiffRapLow[nDiffRapBins/2], mDiffRapHi[nDiffRapBins-1] );
			
			const TString ptName = Form("%1.0f < p_{T}^{#mu#mu} < %1.1f GeV/c",  fabs(ptLow4Fit), ptHig4Fit);

			drawLatex(0.15, 0.86, nCasesName[i_ncase], mTextFont, 0.06, mTextColor);
			drawLatex(.55, 0.86, yName,                mTextFont, 0.05, mTextColor);
			drawLatex(.55, 0.80, ptName,               mTextFont, 0.05, mTextColor);
			drawLatex(.15, 0.75, Form("#chi^{2}/ndf = %1.1f", chi2ndf),                               mTextFont, mTextSize, mTextColor);
			const double textDy = 0.05;
			drawLatex(.55, 0.66+textDy, Form("N_{J/#psi} = %d #pm %d",   TMath::Nint(nJpsi.getVal()), TMath::Nint(nJpsi.getError())), mTextFont, mTextSize, mTextColor);
			drawLatex(.55, 0.61+textDy, Form("N_{#psi(2S)} = %d #pm %d", TMath::Nint(nPsi.getVal()),  TMath::Nint(nPsi.getError()) ), mTextFont, mTextSize, mTextColor);
			drawLatex(.55, 0.55+textDy, Form("N_{QED} = %d #pm %d",      TMath::Nint(nQED.getVal()),  TMath::Nint(nQED.getError()) ), mTextFont, mTextSize, mTextColor);
			drawLatex(.55, 0.48+textDy, Form("N^{in J/#psi}_{QED} = %d #pm %d", (int)nQEDinJpsi, (int)nQEDinJpsiErr ), mTextFont, mTextSize, mTextColor);

			//----------------------------------------------------------------------------------------------------------------------------------------------------
			c1->SaveAs( Form("outplots/massSpec_4ptFitConstrain_"+nCasesName[i_ncase]+"_iy%d.png",  iy) );
			//----------------------------------------------------------------------------------------------------------------------------------------------------
		
			delete ResFit;
			delete hMass;
			delete frameMass;

			//------------------------------------------------------------------------------------------------------------
			//Fit PT Spectral to get fI factor, We only have CohJpsi for 0nXn, XnXn,
			//------------------------------------------------------------------------------------------------------------
			TH1D* hCohJpsiPtHist, *hInCohJpsiPtHist, *hFeeddownJpsiPtHist, *hQEDPtHist;
			if( i_ncase==0  )
			{
				if(iy==nDiffRapBins) hCohJpsiPtHist = (TH1D *)inf_Temps->Get( "hCohJpsiPt" );
				else                 hCohJpsiPtHist = (TH1D *)inf_Temps->Get( Form("hCohJpsiPt_RapBin%d", iy) );
			}
			else
			{
				if(iy==nDiffRapBins) hCohJpsiPtHist = (TH1D *)inf_Temps->Get(      "hCohJpsi_"+nCasesName[i_ncase]+"Pt");
				else                 hCohJpsiPtHist = (TH1D *)inf_Temps->Get( Form("hCohJpsi_"+nCasesName[i_ncase]+"Pt_RapBin%d", iy) );
			}
			
			if(iy==nDiffRapBins)
			{
				hInCohJpsiPtHist    = (TH1D *)inf_Temps->Get( "hInCohJpsiPt"         );
				hFeeddownJpsiPtHist = (TH1D *)inf_Temps->Get( "hCohPsi2SFeeddownPt"  );
				hQEDPtHist          = (TH1D *)inf_Temps->Get( "hLowMassGammaGammaPt" );
			}
			else
			{
				hInCohJpsiPtHist    = (TH1D *)inf_Temps->Get( Form("hInCohJpsiPt_RapBin%d",         iy) );
				hFeeddownJpsiPtHist = (TH1D *)inf_Temps->Get( Form("hCohPsi2SFeeddownPt_RapBin%d",  iy) );
				hQEDPtHist          = (TH1D *)inf_Temps->Get( Form("hLowMassGammaGammaPt_RapBin%d", iy) );
			}
			//--------------------------------------------------------------------------------------

			RooRealVar mPt("mPt", "p_{T} (GeV)", ptLow4Fit, ptHig4Fit);
			
			const int mRebPt = 1;
	
			//--------------------------------------------------------------------------------------
			hCohJpsiPtHist->RebinX(mRebPt);
			RooDataHist hCohJpsiPtRooHist(     "hCohJpsiPtRooHist",      "hCohJpsiPtRooHist",      mPt, hCohJpsiPtHist);
			RooHistPdf  cohJpsiPdf(            "cohJpsiPdf",             "cohJpsiPdf",             mPt, hCohJpsiPtRooHist,      0);
			//--------------------------------------------------------------------------------------

			//--------------------------------------------------------------------------------------
			hInCohJpsiPtHist->RebinX(mRebPt);
			RooDataHist hInCohJpsiPtRooHist(   "hInCohJpsiPtRooHist",    "hInCohJpsiPtRooHist",    mPt, hInCohJpsiPtHist);
			RooHistPdf  incohJpsiPdf(          "incohJpsiPdf",           "incohJpsiPdf",           mPt, hInCohJpsiPtRooHist,    0);
			//--------------------------------------------------------------------------------------

			//--------------------------------------------------------------------------------------
			hFeeddownJpsiPtHist->RebinX(mRebPt);
			RooDataHist hFeeddownJpsiPtRooHist("hFeeddownJpsiPtRooHist", "hFeeddownJpsiPtRooHist", mPt, hFeeddownJpsiPtHist);
			RooHistPdf  feeddownJpsiPdf(       "feeddownJpsiPdf",        "feeddownJpsiPdf",        mPt, hFeeddownJpsiPtRooHist, 0);
			//--------------------------------------------------------------------------------------

			////--------------------------------------------------------------------------------------
			////------------------------gammagamma-->mumu--------------------------------------------------------------
			////--------------------------------------------------------------------------------------
			//if use the simulated tempaltes for gammagamma-->mumu
			hQEDPtHist->RebinX(mRebPt);
			RooDataHist hQEDPtRooHist("hQEDPtRooHist", "hQEDPtRooHist", mPt, hQEDPtHist);
			RooHistPdf  qedPtPdf(     "qedPtPdf",      "qedPtPdf",      mPt, hQEDPtRooHist, 0);

			////if use the Side Bands as tempaltes for gammagamma-->mumu
			//RooDataHist hQEDPtRooHist("hQEDPtRooHist", "hQEDPtRooHist", mPt, hPt_SideBand);
			//RooHistPdf  qedPtPdf(     "qedPtPdf",      "qedPtPdf",      mPt, hQEDPtRooHist, 0);

			//RooConstVar bpd("bpd", "bpd", 1.79);
			//RooConstVar npd("npd", "npd", 3.58);
			RooRealVar bpd("bpd", "bpd", 1.79, 0, 5.0);
			RooRealVar npd("npd", "npd", 3.58, 0, 10);
			RooGenericPdf *dissoJpsiPdf = new RooGenericPdf("dissoJpsiPdf", "dissoJpsiPdf", "mPt*TMath::Power(1+(bpd/npd)*mPt*mPt, -npd)", RooArgSet(mPt, bpd, npd));
		
			const int    higBin_FDW = hFeeddownJpsiPtHist->FindBin( mPtCut4Coh );
			const int    higBin_Coh = hCohJpsiPtHist     ->FindBin( mPtCut4Coh );

			const double frac_FDW = hFeeddownJpsiPtHist  ->Integral(1, higBin_FDW) / hFeeddownJpsiPtHist->Integral(1, hFeeddownJpsiPtHist->GetNbinsX());
			const double frac_Coh = hCohJpsiPtHist       ->Integral(1, higBin_Coh) / hCohJpsiPtHist     ->Integral(1, hCohJpsiPtHist->GetNbinsX()     );

			const double fDValue = fD_inPtCut[i_ncase][iy]*(frac_Coh/frac_FDW); //need to rescale to full pt fD*(pdf_feeddown->Integral()/pdf_cohJpsi->Integral())

			RooConstVar fracPrim(     "fracPrim",       "fracPrim",      1./(1.+fDValue));
			RooRealVar  nCohJpsi_wFDW("nCohJpsi_wFDW",  "nCohJpsi_wFDW", nJpsiValue*0.80, 0, nJpsiValue);
			RooRealVar  nInCohJpsi(   "nInCohJpsi",     "nInCohJpsi",    nJpsiValue*0.04, 0, nJpsiValue*0.20);
			RooRealVar  nDissoJpsi(   "nDissoJpsi",     "nDissoJpsi",    nJpsiValue*0.06, 0, nJpsiValue*0.20);
			RooConstVar nQEDBg(       "nQEDBg",         "nQEDBg",        nQEDinJpsi    );
			
			RooDataHist dataPt("dataPt", "dataPt", mPt, hPt); 
			
			RooAddPdf Pdf_CohJpsi_wFDW("Pdf_CohJpsi_wFDW", "cohJpsiPdf+feeddownJpsiPdf", RooArgList(cohJpsiPdf,feeddownJpsiPdf), fracPrim);

			RooAddPdf totPtPdf("totPtPdf", "totPtPdf", 
					RooArgList( Pdf_CohJpsi_wFDW, incohJpsiPdf, *dissoJpsiPdf, qedPtPdf ),
					RooArgList( nCohJpsi_wFDW,    nInCohJpsi,    nDissoJpsi,   nQEDBg)  );
			
			totPtPdf.fitTo(dataPt,Extended(kTRUE),SumW2Error(kTRUE),Hesse(kTRUE),Minos(kFALSE),Save());
			
			//calcualte numbers
			const double N_CohJpsi       = nCohJpsi_wFDW.getVal()  *fracPrim.getVal();
			const double Nerr_CohJpsi    = nCohJpsi_wFDW.getError()*fracPrim.getVal();
			const double N_InCohJpsi     = nInCohJpsi.getVal();
			const double Nerr_InCohJpsi  = nInCohJpsi.getError();
			const double N_DissoJpsi     = nDissoJpsi.getVal();
			const double Nerr_DissoJpsi  = nDissoJpsi.getError();
			
			//calcualte number of coherent Jpsi within pt<0.20 GeV/c region
			const double N_CohJpsi_inPtCut    = N_CohJpsi   *frac_Coh;
			const double Nerr_CohJpsi_inPtCut = Nerr_CohJpsi*frac_Coh;

			//calculate number of incoherent Jpsi within pt<0.20 GeV/c
			const int    higBin_InCoh         = hInCohJpsiPtHist ->FindBin( mPtCut4Coh );
			const double frac_InCoh           = hInCohJpsiPtHist ->Integral(1, higBin_InCoh) / hInCohJpsiPtHist->Integral(1, hInCohJpsiPtHist->GetNbinsX());
			const double N_InCoh_inPtCut      = N_InCohJpsi    * frac_InCoh;
			const double Nerr_InCoh_inPtCut   = Nerr_InCohJpsi * frac_InCoh;
			
			//calculate number of incoherent Jpsi with n-disso within pt<0.20 GeV/c
			mPt.setRange("CohSignal", 0., mPtCut4Coh );
			RooAbsReal* frac_disso_inPtCut       = dissoJpsiPdf->createIntegral( mPt, NormSet(mPt), Range("CohSignal") );
			const double fracValue_disso_inPtCut = frac_disso_inPtCut->getVal();
			const double N_disso_inPtCut    = N_DissoJpsi   *fracValue_disso_inPtCut;
			const double Nerr_disso_inPtCut = Nerr_DissoJpsi*fracValue_disso_inPtCut;

			const double fI_Value = (N_InCoh_inPtCut + N_disso_inPtCut)/N_CohJpsi_inPtCut;
			const double fI_Error = sqrt( 
					pow(Nerr_InCoh_inPtCut + Nerr_disso_inPtCut,2) / pow(N_InCoh_inPtCut+N_disso_inPtCut,2) 
					+ pow(Nerr_CohJpsi_inPtCut/N_CohJpsi_inPtCut,2)
					)
				*fI_Value;
			

			//ALICE method, where 
			//Double_t f_I   =  (N_I     + N_diss     ) / (N_coh2);
			//Double_t ErrfI = sqrt((N_IError + N_dissError)*(N_IError + N_dissError) / ((N_I + N_diss)*(N_I + N_diss)) + (N_cohError*N_cohError)/(N_coh2*N_coh2) )*f_I;
			//const double 
			
			//calculate expected coherent Jpsi number = N Jpsi from mass fit / (1+fD+fI)
			const double NJpsi_Coh_cal    = NJpsi_inMFit[i_ncase][iy]    / (1. + fD_inPtCut[i_ncase][iy] + fI_Value);
			const double NerrJpsi_Coh_cal = NerrJpsi_inMFit[i_ncase][iy] / (1. + fD_inPtCut[i_ncase][iy] + fI_Value); //need to update with fD, fI uncertainties
			
			const double unit_ub2mb = 1000.;
			const double Lum =  mCMSLum*(unit_ub2mb); //in ub, need to use mb to compare to Alice
			
			//calculate the cross section: xsec = (NJpsiFromMfit/(1+fD+fI))*(1/eff)*(1/BR)*(1/Lum)*(1/dy)
			const double xsecValue = NJpsi_Coh_cal    * (1./Eff_CohJpsi[i_ncase][iy]) * (1./br_Jpsi2uu) * (1./Lum) * (1./dY) * (1./temAcc[iy]);
			const double xsecError = NerrJpsi_Coh_cal * (1./Eff_CohJpsi[i_ncase][iy]) * (1./br_Jpsi2uu) * (1./Lum) * (1./dY) * (1./temAcc[iy]);

			int nFramePtBins = (ptHig4Fit - ptLow4Fit)/hPt->GetBinWidth(1);
			cout<<nFramePtBins<<endl;

			RooPlot *framePt = mPt.frame(Range(ptLow4Fit, ptHig4Fit), Title(""), Bins(nFramePtBins));

			dataPt  .plotOn(framePt, MarkerStyle(20), MarkerSize(0.6), MarkerColor(1), LineColor(1), LineWidth(mLineWidth), DrawOption("pz"));
			totPtPdf.plotOn(framePt, LineColor(1), LineStyle(1), LineWidth(mLineWidth));
			totPtPdf.plotOn(framePt, Components(RooArgSet(cohJpsiPdf)),      LineColor(cohJpsiColor),      LineStyle(cohJpsiStyle),      LineWidth(mLineWidth));
			totPtPdf.plotOn(framePt, Components(RooArgSet(feeddownJpsiPdf)), LineColor(feeddownJpsiColor), LineStyle(feeddownJpsiStyle), LineWidth(mLineWidth));
			totPtPdf.plotOn(framePt, Components(RooArgSet(incohJpsiPdf)),    LineColor(incohJpsiColor),    LineStyle(incohJpsiStyle),    LineWidth(mLineWidth));
			totPtPdf.plotOn(framePt, Components(RooArgSet(*dissoJpsiPdf)),   LineColor(dissoJpsiColor),    LineStyle(dissoJpsiStyle),    LineWidth(mLineWidth));
			totPtPdf.plotOn(framePt, Components(RooArgSet(qedPtPdf)),        LineColor(qedColor),          LineStyle(qedStyle),          LineWidth(mLineWidth));

			cout<<endl;
			cout<<"******** Print frame ********"<<endl;
			framePt->Print();
			cout<<"******** End ********"<<endl;
			cout<<endl;
			
			chi2ndf = framePt->chiSquare("totPtPdf_Norm[mPt]", "h_dataPt", 5);
			////chi2ndf = framePt->chiSquare("totPtPdf_Norm[mPt]", "h_dataPt", 7);
	
			c2 -> Divide(1,2);
			c2 -> cd(1);
			gPad->SetPad(0.0,0.25,1.0,0.96);
			gPad->SetBottomMargin(0);
			gPad->SetRightMargin(0.05);
			gPad->SetTopMargin(0);
			//setPad(0.12, 0.08, 0.07, 0.13);

			//gPad->SetLogy(0);
			//framePt->GetYaxis()->SetRangeUser(0.5, hPt->GetMaximum()*1.3);
			gPad->SetLogy(1);
			gPad->SetLogx(0);
			framePt->GetYaxis()->SetRangeUser(0.5, hPt->GetMaximum()*5);
			framePt->Draw() ;
			drawLatex(0.15, 0.92, "CMS Pb-Pb #sqrt{s_{NN}} = 5.02 TeV UPC ( "+nCasesName[i_ncase]+" )",  42,        0.06,      mTextColor );
			drawLatex(0.18, 0.82, yName,                                                                 mTextFont, 0.05,      mTextColor );
			const double textDy2 = 0.05;
			drawLatex(0.25, 0.70+textDy2, Form("For p_{T}<%.2f GeV/c:",              mPtCut4Coh),                                   mTextFont, 0.035, mTextColor);
			drawLatex(0.25, 0.65+textDy2, Form("N^{Coh}_{J/#psi} = %d #pm %d",  (int)N_CohJpsi_inPtCut, (int)Nerr_CohJpsi_inPtCut), mTextFont, 0.035, mTextColor);
			drawLatex(0.25, 0.60+textDy2, Form("f_{I} = (N^{All}_{InCoh}/N^{Coh}_{J/#psi}) = %.3f #pm %.3f", fI_Value, fI_Error),     mTextFont, 0.035, mTextColor);
			drawLatex(0.25, 0.52+textDy2, Form("N^{in Mfit}_{J/#psi}/(1+f_{I}+f_{D}) = %d #pm %d",  (int)NJpsi_Coh_cal, (int)NerrJpsi_Coh_cal), mTextFont, 0.035, mTextColor);
			drawLatex(0.35, 0.45+textDy2, Form("#frac{#sigma^{Coh}_{J/#psi}}{dy} = %.3f #pm %.3f (mb)", xsecValue, xsecError ),                                     mTextFont, 0.035, mTextColor);

			TLegend  *leg =  new TLegend(0.62, 0.45, 0.88, 0.80);
			leg->SetFillStyle(0);
			leg->SetFillColor(0);
			leg->SetTextFont(mTextFont);
			leg->SetTextSize(0.035);
			leg->AddEntry(framePt->findObject("h_dataPt"),            "Data",        "p");
			leg->AddEntry(framePt->findObject("totPtPdf_Norm[mPt]"),   Form("Total fit: #chi^{2}/ndf = %1.1f", chi2ndf),   "l");
//(RooHist::h_dataPt,RooCurve::totPtPdf_Norm[mPt],RooCurve::totPtPdf_Norm[mPt]_Comp[cohJpsiPdf],RooCurve::totPtPdf_Norm[mPt]_Comp[feeddownJpsiPdf],RooCurve::totPtPdf_Norm[mPt]_Comp[incohJpsiPdf],RooCurve::totPtPdf_Norm[mPt]_Comp[dissoJpsiPdf],RooCurve::totPtPdf_Norm[mPt]_Comp[qedPtPdf])
			const TString curveName[5]  = {"totPtPdf_Norm[mPt]_Comp[cohJpsiPdf]","totPtPdf_Norm[mPt]_Comp[incohJpsiPdf]", "totPtPdf_Norm[mPt]_Comp[dissoJpsiPdf]", "totPtPdf_Norm[mPt]_Comp[feeddownJpsiPdf]","totPtPdf_Norm[mPt]_Comp[qedPtPdf]"};
			const TString curveTitle[5] = {"Coherent J/#psi", "Incoherent J/#psi", "Incoherent J/#psi with disso.", "Coherent #psi' #rightarrow J/#psi+X", "#gamma#gamma #rightarrow #mu#mu"};
			for(int icv=0; icv<5; icv++) 
			{
				//if(icv==0||icv==3) continue;
				leg->AddEntry(framePt->findObject( curveName[icv]), curveTitle[icv], "l");
			}

			leg->Draw("same");

			c2	->cd(2);
			gPad->SetPad(0.0,0.0,1.0,0.25);
			gPad->SetTopMargin(0);
			gPad->SetRightMargin(0.05);
			gPad->SetBottomMargin(0.45);

			//calculate the fiting pull: (Data-FitCurve)/#Sigma_{Data}
			RooHist *hpull_pt = framePt->pullHist("h_dataPt", "totPtPdf_Norm[mPt]");
			hpull_pt ->SetMarkerStyle(24);
			hpull_pt ->SetMarkerSize(0.6);
			hpull_pt ->SetMarkerColor(1);
			hpull_pt ->SetLineColor(1);
			hpull_pt ->SetLineWidth(1);

			RooPlot *framePtPull = mPt.frame( Range(ptLow4Fit, ptHig4Fit), Title(""), Bins(nFramePtBins) );
			framePtPull ->addPlotable(hpull_pt, "pz");
			framePtPull ->GetYaxis()->SetRangeUser(-14.0, 14.0);
			framePtPull ->SetYTitle("#frac{Data-Fit}{#sigma_{Data}}");
			framePtPull ->SetXTitle("p_{T} (GeV/c)");
			framePtPull ->GetYaxis()->CenterTitle();
			framePtPull ->GetYaxis()->SetNdivisions(6);
			framePtPull ->GetYaxis()->SetTitleSize(0.15);
			framePtPull ->GetYaxis()->SetTitleOffset(0.33);
			framePtPull ->GetYaxis()->SetLabelSize(0.10);
			//framePtPull ->GetYaxis()->SetLabelFont(20);
			framePtPull ->GetXaxis()->SetTitleSize(0.20);
			framePtPull ->GetXaxis()->SetTitleOffset(0.95);
			framePtPull ->GetXaxis()->SetLabelSize(0.16);
			framePtPull ->GetXaxis()->SetLabelFont(40);
			framePtPull ->SetTickLength(0.08);
			framePtPull ->Draw() ;
			
			//	TH2D* htem2d_4pt = new TH2D("htem2d_4pt", "", nFramePtBins, ptLow, ptHi, 10, -0.001, 0.001);
			//	htem2d_4pt->SetYTitle("Data/Fit");
			//	htem2d_4pt->SetXTitle("p_{T} (GeV/c)");
			//	htem2d_4pt->GetYaxis()->SetNdivisions(4);
			//	htem2d_4pt->GetYaxis()->SetTitleSize(0.22);
			//	htem2d_4pt->GetYaxis()->SetTitleOffset(0.25);
			//	htem2d_4pt->GetYaxis()->SetLabelSize(0.15);
			//	htem2d_4pt->GetYaxis()->SetLabelFont(40);
			//	htem2d_4pt->GetXaxis()->SetTitleSize(0.20);
			//	htem2d_4pt->GetXaxis()->SetTitleOffset(0.95);
			//	htem2d_4pt->GetXaxis()->SetLabelSize(0.16);
			//	htem2d_4pt->GetXaxis()->SetLabelFont(40);
			//	htem2d_4pt->SetTickLength(0.08);
			//	htem2d_4pt->Draw() ;
			//
			//	histPull_pt->Draw("lsame");
			//	histPull_pt->Draw("pesame");
			//
			c2->SaveAs( Form("outplots/ptSpec_4decouple_"+nCasesName[i_ncase]+"_iy%d.png",  iy) );
			c2->SaveAs( Form("outplots/ptSpec_4decouple_"+nCasesName[i_ncase]+"_iy%d.pdf",  iy) );

			c2 -> Clear();
		}//iy
	}//ixn
	delete c2;
	
	cout << "End of program !" << endl;
}
