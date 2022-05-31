#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"

Double_t CrystalBall2(Double_t *x, Double_t *par);

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
Int_t    mTextFont  = 42;
Double_t mTextSize  = 0.035;
Int_t    mTextColor = 1;

Double_t mMarkerStyle = 20;
Double_t mMarkerSize  = 0.8;

Double_t mTitleSize    = 0.06;
Double_t mXTitleOffset = 0.95;
Double_t mYTitleOffset = 0.95;
Double_t mLabelSize    = 0.05;
Double_t mTickLength   = 0.02;
Int_t    mXNdivisions  = 210;
Int_t    mYNdivisions  = 208;


void fitTemp()
{
	TFile *inf_Temps = new TFile("../simulation/effAndTemp/MassPtTemp_AllSpecs.root", "read");
	TH1D  *hPt_CohJpsi_Sim = (TH1D *)inf_Temps->Get("hCohJpsiPt"          );

	TFile* infile_data = new TFile("./outplots/outdata_pt_QED_CohJpsi.root", "read"); //Coherent Jpsi pt shape from tiny pt bin fitting, only use pt<0.10 GeV/c region to find stretch factor
	TH1D* hPt_CohJpsi_fromDataFit = (TH1D*) infile_data->Get("hPt_CohJpsi");
	
	gStyle->SetFrameLineWidth(3);
	TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);

	const double x1 = 0.;//0.005;
	const double x2 = 0.400;

	hPt_CohJpsi_Sim ->Scale(1./hPt_CohJpsi_Sim->Integral());
	hPt_CohJpsi_Sim ->GetYaxis()->SetTitleSize(0.05);
	hPt_CohJpsi_Sim ->GetYaxis()->SetTitleOffset(1.15);
	hPt_CohJpsi_Sim ->GetXaxis()->SetTitleSize(0.05);
	hPt_CohJpsi_Sim ->GetXaxis()->SetTitleOffset(1.10);
	hPt_CohJpsi_Sim ->SetMarkerStyle(24);
	hPt_CohJpsi_Sim ->SetMarkerSize(0.8);
	hPt_CohJpsi_Sim ->SetMarkerColor(1);
	hPt_CohJpsi_Sim ->SetLineColor(1);
	hPt_CohJpsi_Sim ->SetLineWidth(2);
	hPt_CohJpsi_Sim ->SetAxisRange(x1,    x2,    "x");
	hPt_CohJpsi_Sim ->SetAxisRange(3.e-7, 0.535,     "y");
	hPt_CohJpsi_Sim ->SetYTitle("");
	//hPt_CohJpsi_Sim ->Draw("pe");

	//1  a            2.74758e+00   9.98959e+04   3.79448e-05  -1.56031e+00
	//   2  b            1.30268e+02   5.12593e+04   5.90896e-05  -1.64300e+00
	//   3  mu           4.51289e-04   7.97396e-01   9.73847e-04  -1.52831e+00
	//   4  nNormData    9.36259e-01   5.88699e+02   3.67113e-05  -1.50959e+00
	//   5  sigma        1.81133e-01   5.33221e-01   2.94832e-04  -6.91552e-01
	//   6  sigma2       4.51223e-02   9.16081e-01   9.56864e-05  -1.14269e+00
	//   ERR DEF= 0.5


	RooRealVar mPt("mPt", "p_{T}^{#mu#mu} (GeV)", x1, x2);
	RooRealVar nNorm1(      "nNorm1",      "nNorm1",      1.e-2, 1.e-5, 1.e3);
	RooRealVar mu1(         "mu1",         "mu1",         0.0004, 0.,    0.1 );
	RooRealVar sigma1(      "sigma1",      "sigma1",      0.04,  0,     1. );
	RooRealVar ratioSigma1( "ratioSigma1", "ratioSigma1", 5.0,   0.,    10.0 );
	RooRealVar a1( "a1", "a1",  0.95, 0., 1.);

	RooGenericPdf *ptPdf4Sim  = new RooGenericPdf("ptPdf4Sim",  "ptPdf4Sim",  "mPt*(a1*TMath::Exp(-pow(mPt-mu1,2)/(2.*pow(sigma1,2))) + (1.-a1)*TMath::Exp( -pow(mPt-mu1,2)/(2.*pow(sigma1*ratioSigma1,2)) ))", RooArgSet(mPt, mu1,sigma1,ratioSigma1, a1) ); //2 gaussion share same mean but diff. width
	//RooGenericPdf *ptPdf = new RooGenericPdf("ptPdf", "ptPdf", "mPt*(a*TMath::Exp(-pow(mPt-mu,2)/(2.*pow(sigma,2))) + b*TMath::Exp( -pow(mPt-mu,2)/(2.*pow(sigma*ratioSigma,2)) ))", RooArgSet(mPt, mu,sigma,ratioSigma, a,b) );
	//RooGenericPdf *ptPdf = new RooGenericPdf("ptPdf", "ptPdf", "mPt*(a*TMath::Exp(-pow(mPt-mu,2)/(2.*pow(sigma,2)))+b*TMath::Exp( -pow(mPt-mu2,2)/(2.*pow(sigma2,2)) ))", RooArgSet(mPt, mu,sigma, mu2, sigma2, a,b) );
	
	RooAddPdf totPtPdf4Sim( "totPtPdf4Sim",  "totPtPdf4Sim",  RooArgList(*ptPdf4Sim),  RooArgList(nNorm1)  );
	//
	//------------------------------------------------------------------------------------------------------------
	RooDataHist dataPt_inSim("dataPt_inSim", "dataPt_inSim", mPt, hPt_CohJpsi_Sim );

	totPtPdf4Sim.fitTo( dataPt_inSim, Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE), Save());
	//------------------------------------------------------------------------------------------------------------
	
	c1->cd();
	setPad(0.13, 0.05, 0.05, 0.13);
	gPad->SetLogy(1);
	int nFrameMBins = 10; //(x2-x1)/hPt_CohJpsi_Sim->GetBinWidth(1);

	RooPlot *framePt_inSim = mPt.frame(Range(x1, x2), Title(""), Bins(nFrameMBins));
	//framePt_inSim ->GetYaxis()->SetRangeUser(1.e-4, hPt_CohJpsi_Sim->GetMaximum()*1.5);
	framePt_inSim ->SetMaximum(hPt_CohJpsi_Sim->GetMaximum()*1.0);
	//framePt_inSim ->SetMinimum(1.e-3);
	framePt_inSim ->GetYaxis()->SetTitleOffset(0.80);

	framePt_inSim ->GetYaxis()->SetNdivisions(8);
	framePt_inSim ->GetYaxis()->SetTitleSize(0.05);
	framePt_inSim ->GetYaxis()->SetLabelSize(0.05);
	framePt_inSim ->GetYaxis()->SetTitleOffset(1.25);

	framePt_inSim ->GetXaxis()->SetNdivisions(8);
	framePt_inSim ->GetXaxis()->SetTitleSize(0.05);
	framePt_inSim ->GetXaxis()->SetLabelSize(0.05);
	framePt_inSim ->GetXaxis()->SetTitleOffset(1.15);
	framePt_inSim ->GetXaxis()->SetLabelFont(40);
	framePt_inSim ->SetTickLength(0.02);

	dataPt_inSim.plotOn( framePt_inSim, MarkerStyle(24), MarkerSize(0.6), MarkerColor(1), LineColor(1), LineWidth(2), DrawOption("pz"));
	totPtPdf4Sim.plotOn( framePt_inSim, LineColor(2), LineStyle(1), LineWidth(2));
	//dataPt_DataFit.plotOn( framePt_inSim, MarkerStyle(20), MarkerSize(0.8), MarkerColor(4), LineColor(4), LineStyle(1), LineWidth(2), DrawOption("pz"));
	//totPtPdfData  .plotOn( framePt_inSim, LineColor(2), LineStyle(1), LineWidth(2));

	//	cout<<endl;
	//	cout<<"******** Print frame ********"<<endl;
	//	framePt_inSim->Print();
	//	cout<<"******** End ********"<<endl;
	//	cout<<endl;

	framePt_inSim->Draw() ;
	//	RooRealVar nNorm1(      "nNorm1",      "nNorm1",      1.e-2, 1.e-5, 1.e3);
	//	RooRealVar mu1(         "mu1",         "mu1",         0.0004, 0.,    0.1 );
	//	RooRealVar sigma1(      "sigma1",      "sigma1",      0.04,  0,     1. );
	//	RooRealVar ratioSigma1( "ratioSigma1", "ratioSigma1", 5.0,   0.,    10.0 );
	//	RooRealVar a1( "a1", "a1",  0.95, 0., 1.);
	
	drawLatex(0.15, 0.90, "fit by pt*(N_{1}*g1(#mu,#sigma_{1}, pt)+N2*g2(#mu,#sigma_{1}*ratio, pt))",  mTextFont, mTextSize, mTextColor);
	drawLatex(0.55, 0.85, Form("#mu = %.5f #pm %.5f",        mu1.getVal(),    mu1.getError() ),        mTextFont, mTextSize, mTextColor);
	drawLatex(0.55, 0.80, Form("#sigma_{1} = %.5f #pm %.5f", sigma1.getVal(), sigma1.getError() ),     mTextFont, mTextSize, mTextColor);
	drawLatex(0.55, 0.72, Form("#frac{#sigma_{2}}{#sigma_{1}} = %.2f #pm %.3f", ratioSigma1.getVal(), ratioSigma1.getError() ), mTextFont, mTextSize, mTextColor);
	drawLatex(0.55, 0.62, Form("#frac{N_{1}}{N_{1}+N_{2}} = %.2f #pm %.3f",   a1.getVal(),    a1.getError() ),        mTextFont, mTextSize, mTextColor);

	//Double_t chi2ndf = framePt_inSim->chiSquare("totPtPdf_Norm[mPt]", "h_dataPt", 6); //h_dataPt,RooCurve::totPtPdf_Norm[mPt]
	//drawLatex(0.18, 0.72, Form("#chi^{2}/ndf = %1.1f", chi2ndf),                                 mTextFont, mTextSize, mTextColor);

	//------------------------------------------------------------------------------------------------------------
	c1->SaveAs("./plots_4fittemp/fitCohJpsi_inSim.pdf");
	//------------------------------------------------------------------------------------------------------------

	//------------------------------------------------------------------------------------------------------------
	RooRealVar mPt2("mPt2", "p_{T}^{#mu#mu} (GeV)", x1, 0.10);
	RooRealVar nNorm2( "nNorm2", "nNorm2",  1.000, 1.e-5, 1.e3);
	RooRealVar mu2(    "mu2",    "mu2",     0.004, 0.,    0.1 );
	RooRealVar sigma2( "sigma2", "sigma2",  0.04,   0,     0.1 );
	//RooRealVar ratioSigma2( "ratioSigma2", "ratioSigma2", 5.2,  0.,     10.0 );
	RooConstVar ratioSigma2("ratioSigma2", "ratioSigma2", ratioSigma1.getVal());
	RooRealVar a2( "a2", "a2",  0.900, 0., 1.);
	RooGenericPdf *ptPdf4Data = new RooGenericPdf("ptPdf4Data", "ptPdf4Data", "mPt2*(a2*TMath::Exp(-pow(mPt2-mu2,2)/(2.*pow(sigma2,2))) + (1.-a2)*TMath::Exp( -pow(mPt2-mu2,2)/(2.*pow(sigma2*ratioSigma2,2)) ))", RooArgSet(mPt2, mu2,sigma2,ratioSigma2, a2) );
	
	RooAddPdf totPtPdf4Data("totPtPdf4Data", "totPtPdf4Data", RooArgList(*ptPdf4Data), RooArgList(nNorm2) );
	//------------------------------------------------------------------------------------------------------------
	hPt_CohJpsi_fromDataFit -> Scale( 1./hPt_CohJpsi_fromDataFit->Integral() );
	RooDataHist dataPt_inData("dataPt_inData", "dataPt_inData", mPt2, hPt_CohJpsi_fromDataFit);
	//------------------------------------------------------------------------------------------------------------
	totPtPdf4Data.fitTo( dataPt_inData, Extended(kTRUE), SumW2Error(kTRUE), Hesse(kTRUE), Minos(kFALSE), Save());
	//------------------------------------------------------------------------------------------------------------
	
	c1->cd();
	setPad(0.13, 0.05, 0.05, 0.13);
	gPad->SetLogy(1);

	RooPlot *framePt_inData = mPt2.frame(Range(x1, x2), Title(""), Bins(nFrameMBins));
	//framePt_inData ->GetYaxis()->SetRangeUser(1.e-4, hPt_CohJpsi_Sim->GetMaximum()*1.5);
	framePt_inData ->SetMaximum(hPt_CohJpsi_Sim->GetMaximum()*1.0);
	//framePt_inData ->SetMinimum(1.e-3);
	framePt_inData ->GetYaxis()->SetTitleOffset(0.80);
	framePt_inData ->GetYaxis()->SetNdivisions(8);
	framePt_inData ->GetYaxis()->SetTitleSize(0.05);
	framePt_inData ->GetYaxis()->SetLabelSize(0.05);
	framePt_inData ->GetYaxis()->SetTitleOffset(1.25);

	framePt_inData ->GetXaxis()->SetNdivisions(8);
	framePt_inData ->GetXaxis()->SetTitleSize(0.05);
	framePt_inData ->GetXaxis()->SetLabelSize(0.05);
	framePt_inData ->GetXaxis()->SetTitleOffset(1.15);
	framePt_inData ->GetXaxis()->SetLabelFont(40);
	framePt_inData ->SetTickLength(0.02);

	dataPt_inData.plotOn( framePt_inData, MarkerStyle(24), MarkerSize(0.6), MarkerColor(1), LineColor(1), LineWidth(2), DrawOption("pz"));
	totPtPdf4Data.plotOn( framePt_inData, LineColor(2), LineStyle(1), LineWidth(2));
	//dataPt_DataFit.plotOn( framePt_inData, MarkerStyle(20), MarkerSize(0.8), MarkerColor(4), LineColor(4), LineStyle(1), LineWidth(2), DrawOption("pz"));
	//totPtPdfData  .plotOn( framePt_inData, LineColor(2), LineStyle(1), LineWidth(2));

	//	cout<<endl;
	//	cout<<"******** Print frame ********"<<endl;
	//	framePt_inData->Print();
	//	cout<<"******** End ********"<<endl;
	//	cout<<endl;


	//Double_t chi2ndf = framePt_inData->chiSquare("totPtPdf_Norm[mPt2]", "h_dataPt", 6); 
	//drawLatex(0.18, 0.72, Form("#chi^{2}/ndf = %1.1f", chi2ndf),                                 mTextFont, mTextSize, mTextColor);
	//
	
	//	RooConstVar nNorm3( "nNorm3", "nNorm3",  nNorm2.getVal()*101 );
	//	RooConstVar mu3(    "mu3",    "mu3",     mu2.getVal()    );
	//	RooConstVar sigma3( "sigma3", "sigma3",  sigma2.getVal() );
	//	//RooConstVar ratioSigma2( "ratioSigma2", "ratioSigma2", 5.2,  0.,     10.0 );
	//	RooConstVar ratioSigma("ratioSigma", "ratioSigma", ratioSigma1.getVal() );
	//	RooConstVar a3( "a3", "a3",  a2.getVal() );
	//	RooGenericPdf *ptPdf4Data_fix = new RooGenericPdf("ptPdf4Data_fix", "ptPdf4Data_fix", "mPt2*(a3*TMath::Exp(-pow(mPt2-mu3,2)/(2.*pow(sigma3,2))) + (1.-a3)*TMath::Exp( -pow(mPt2-mu3,2)/(2.*pow(sigma3*ratioSigma,2)) ))", RooArgSet(mPt2, mu3,sigma3,ratioSigma, a3) );
	//	
	//	ptPdf4Data_fix->plotOn( framePt_inData, LineColor(5), LineStyle(5), LineWidth(2));


	framePt_inData->Draw() ;
	
	drawLatex(0.15, 0.90, "fit by pt*(N_{1}*g1(#mu,#sigma_{1}, pt)+N2*g2(#mu,#sigma_{1}*FixRatio, pt))",  mTextFont, mTextSize, mTextColor);
	drawLatex(0.55, 0.85, Form("#mu = %.5f #pm %.5f",        mu2.getVal(),    mu2.getError() ),        mTextFont, mTextSize, mTextColor);
	drawLatex(0.55, 0.80, Form("#sigma_{1} = %.5f #pm %.5f", sigma2.getVal(), sigma2.getError() ),     mTextFont, mTextSize, mTextColor);
	drawLatex(0.55, 0.72, Form("#frac{#sigma_{2}}{#sigma_{1}} = %.2f", ratioSigma2.getVal() ), mTextFont, mTextSize, mTextColor);
	drawLatex(0.55, 0.62, Form("#frac{N_{1}}{N_{1}+N_{2}} = %.2f #pm %.3f",   a2.getVal(),    a2.getError() ),        mTextFont, mTextSize, mTextColor);
	//------------------------------------------------------------------------------------------------------------
	c1->SaveAs("./plots_4fittemp/fitCohJpsi_inData.pdf");
	//------------------------------------------------------------------------------------------------------------

	//overlay the extrapolate function with the simulated templates and its fit and see what is the difference
	c1->cd();
	
	RooConstVar nNorm3( "nNorm3", "nNorm3",  nNorm2.getVal()*101 );
	RooConstVar mu3(    "mu3",    "mu3",     mu2.getVal()    );
	RooConstVar sigma3( "sigma3", "sigma3",  sigma2.getVal() );
	//RooConstVar ratioSigma2( "ratioSigma2", "ratioSigma2", 5.2,  0.,     10.0 );
	RooConstVar ratioSigma("ratioSigma", "ratioSigma", ratioSigma1.getVal() );
	RooConstVar a3( "a3", "a3",  a2.getVal() );
	RooGenericPdf *ptPdf4Data_fix = new RooGenericPdf("ptPdf4Data_fix", "ptPdf4Data_fix", "mPt*(a3*TMath::Exp(-pow(mPt-mu3,2)/(2.*pow(sigma3,2))) + (1.-a3)*TMath::Exp( -pow(mPt-mu3,2)/(2.*pow(sigma3*ratioSigma,2)) ))", RooArgSet(mPt, mu3,sigma3,ratioSigma, a3) );
	//RooDataHist dataPt_inData_cp("dataPt_inData_cp", "dataPt_inData_cp", mPt, hPt_CohJpsi_fromDataFit);
	//dataPt_inData_cp.plotOn(framePt_inSim, MarkerStyle(24), MarkerSize(0.6), MarkerColor(4), LineColor(4), LineWidth(2), DrawOption("pz") );
	ptPdf4Data_fix->plotOn( framePt_inSim, LineColor(4), LineStyle(1), LineWidth(2));
	framePt_inSim->Draw() ;
	

	cout<<endl;
	cout<<"******** Print frame ********"<<endl;
	framePt_inSim->Print();
	cout<<"******** End ********"<<endl;
	cout<<endl;

	//(RooHist::h_dataPt_inSim,RooCurve::totPtPdf4Sim_Norm[mPt],RooCurve::ptPdf4Data_fix_Norm[mPt])

	TLegend  *leg =  new TLegend(0.44, 0.65, 0.84, 0.88);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextFont(mTextFont);
	leg->SetTextSize(0.035);
	leg->AddEntry( framePt_inSim->findObject("h_dataPt_inSim"),            "Pt Shape from STARLight",   "p");
	leg->AddEntry( framePt_inSim->findObject("totPtPdf4Sim_Norm[mPt]"),    "fitted by pt*(a*g1+b*g2)",  "l");
	leg->AddEntry( framePt_inSim->findObject("ptPdf4Data_fix_Norm[mPt]"),  "curve fitted to Real Data", "l");
	leg->Draw("same");

	
	c1->SaveAs("./plots_4fittemp/fitCohJpsi_Overlay.pdf");
//
}

//------------------------------------------------------
Double_t CrystalBall2(Double_t *x, Double_t *par)
{
	Double_t N      = par[0];
	Double_t mu     = par[1];
	Double_t s      = par[2];
	Double_t n1     = par[3];
	Double_t alpha1 = par[4];
	Double_t n2     = par[5];
	Double_t alpha2 = par[6];

	Double_t A = TMath::Power(n1/fabs(alpha1), n1) * TMath::Exp(-alpha1*alpha1/2.);
	Double_t B = n1/fabs(alpha1) - fabs(alpha1);

	Double_t C = TMath::Power(n2/fabs(alpha2), n2) * TMath::Exp(-alpha2*alpha2/2.);
	Double_t D = n2/fabs(alpha2) - fabs(alpha2);

	Double_t norm = (x[0]-mu)/s;

	if(norm < -alpha1) 
	{
		return N * A * TMath::Power(B-norm, -n1);
	}
	else if(norm < alpha2) 
	{
		return N * TMath::Exp(-0.5*norm*norm);
	}
	else 
	{
		return N * C * TMath::Power(D+norm, -n2);
	}
}

