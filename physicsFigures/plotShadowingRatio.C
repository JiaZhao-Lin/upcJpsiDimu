#include "fit2D.C"
#include "../common/LoadSignal.C"
#include "readTheory.C"
#include "Map_IO.C"

const double Sqrt_s     = 5020;
const double Gamma_beam = 2672.9;
const double Mass_N     = (0.93827+0.93957) / 2;
// the effective Lorentz factor in the rest frame of the collision partner A. It is related to that in the c.m. system of the two colliding heavy ions
// Eq.4.2 Prog. Parr. Nucl. Phys.. Vol. 39, pp. 503-564. 1997
const double Gamma_eff  = 2 * Gamma_beam * Gamma_beam -1;

struct Converter
{
	static double y2x(const double y)
	{
		return (mJpsi_PDG / Sqrt_s) * exp(-y);
	}
	
	static double y2W(const double y)
	{
		return sqrt( (2 * Gamma_beam * Mass_N * mJpsi_PDG) * exp(y) );
	}

	static double W2y(const double W)
	{
		return log( W * W / (2 * Gamma_beam * Mass_N * mJpsi_PDG) );
	}

	static double y2w(const double y)
	{
		//Rap to photon energy
		return mJpsi_PDG / 2 * exp(y);
	}

	static double y2w_eff(const double y)
	{
		//Rap to photon energy in A rest frame
		return w2w_eff( y2w(y) );
	}

	static double w2y(const double w)
	{
		//Photon energy to Rap
		return log(2 * w /mJpsi_PDG);
	}

	static double w2w_eff(const double w)
	{
		//convert photon energy in gamma-A cm frame to A rest frame
		return w * Gamma_eff;
	}

	static std::vector<double> W2w_eff_v( std::vector<double> W)
	{
		std::vector<double> v_;
		//convert photon energy in gamma-A cm frame to A rest frame
		for (int i = 0; i < W.size(); i++)
		{
			v_.push_back( y2w_eff(W2y( W[i])) );
		}
		return v_;
	}
};


void Cal_R_Error(std::map<TString, std::vector<double>> &ShadowRatio_ParamsMap)
{
	for (int i = 0; i < ShadowRatio_ParamsMap.at("Raps").size(); ++i)
	{
		ShadowRatio_ParamsMap.at("R_Err").push_back( 0.5 * ShadowRatio_ParamsMap.at("R")[i] * ShadowRatio_ParamsMap.at("Sigmas_Err")[i]/ShadowRatio_ParamsMap.at("Sigmas")[i]);
	}
}

void plotSigmaVsW( std::map<TString, std::vector<double>> &ShadowRatio_ParamsMap, const std::map<TString, std::vector<double>> &TotalSysUncer_Map, const int flag4Axis = 0 ) //0: logY only, 1: logX and LogY
{
	for (int i = 0; i < ALICE_Run2_FwdRap_y.size(); ++i)
	{
		//Calculating the values
		ALICE_Run2_FwdRap_W        		.push_back( Converter::y2W(ALICE_Run2_FwdRap_y[i]) );
		// ALICE_Run2_FwdRap_w_eff        	.push_back( Converter::y2w_eff(  ALICE_Run2_FwdRap_y[i]  ) );
		ALICE_Run2_FwdRap_Sigma    		.push_back( ALICE_Run2_FwdRap_dSigmady[i] / ALICE_Run2_FwdRap_Flux[i] );
		ALICE_Run2_FwdRap_Sigma_StatErr	.push_back( ALICE_Run2_FwdRap_Sigma[i] * ALICE_Run2_FwdRap_dSigmady_StatErr[i] / ALICE_Run2_FwdRap_dSigmady[i]		);
		ALICE_Run2_FwdRap_Sigma_SysErrLow.push_back( ALICE_Run2_FwdRap_Sigma[i] * TMath::Hypot(ALICE_Run2_FwdRap_dSigmady_SysErrLow[i] / ALICE_Run2_FwdRap_dSigmady[i],
																							 ALICE_Run2_FwdRap_Flux_Err[i] / ALICE_Run2_FwdRap_Flux[i])		);
		ALICE_Run2_FwdRap_Sigma_SysErrHig.push_back( ALICE_Run2_FwdRap_Sigma[i] * TMath::Hypot(ALICE_Run2_FwdRap_dSigmady_SysErrHig[i] / ALICE_Run2_FwdRap_dSigmady[i],
																							 ALICE_Run2_FwdRap_Flux_Err[i] / ALICE_Run2_FwdRap_Flux[i])		);
		
		ALICE_Run2_FwdRap_x        		.push_back( Converter::y2x(ALICE_Run2_FwdRap_y[i]) );
		ALICE_Run2_FwdRap_R 			.push_back( sqrt(ALICE_Run2_FwdRap_Sigma[i]/ALICE_Run2_FwdRap_Sigma_IA[i])	);
		ALICE_Run2_FwdRap_R_StatErr		.push_back( 0.5 * ALICE_Run2_FwdRap_R[i] * ALICE_Run2_FwdRap_Sigma_StatErr[i]/ALICE_Run2_FwdRap_Sigma[i]	);
		ALICE_Run2_FwdRap_R_SysErrLow	.push_back( 0.5 * ALICE_Run2_FwdRap_R[i] * TMath::Hypot( ALICE_Run2_FwdRap_Sigma_SysErrLow[i]/ALICE_Run2_FwdRap_Sigma[i],
																								ALICE_Run2_FwdRap_Sigma_IA_Err[i]/ALICE_Run2_FwdRap_Sigma_IA[i])	);
		ALICE_Run2_FwdRap_R_SysErrHig	.push_back( 0.5 * ALICE_Run2_FwdRap_R[i] * TMath::Hypot( ALICE_Run2_FwdRap_Sigma_SysErrHig[i]/ALICE_Run2_FwdRap_Sigma[i],
																								ALICE_Run2_FwdRap_Sigma_IA_Err[i]/ALICE_Run2_FwdRap_Sigma_IA[i])	);
	}
	for (int i = 0; i < LHCb_Run2_FwdRap_y.size(); ++i)
	{
		//Calculating the values
		LHCb_Run2_FwdRap_W        		.push_back( Converter::y2W(LHCb_Run2_FwdRap_y[i]) );
		// LHCb_Run2_FwdRap_w_eff     		.push_back( Converter::y2w_eff(  LHCb_Run2_FwdRap_y[i]  ) ) ;
		LHCb_Run2_FwdRap_Sigma    		.push_back( LHCb_Run2_FwdRap_dSigmady[i] / LHCb_Run2_FwdRap_Flux[i] );
		LHCb_Run2_FwdRap_Sigma_StatErr	.push_back( LHCb_Run2_FwdRap_Sigma[i] * LHCb_Run2_FwdRap_dSigmady_StatErr[i] / LHCb_Run2_FwdRap_dSigmady[i]		);
		LHCb_Run2_FwdRap_Sigma_SysErrLow.push_back( LHCb_Run2_FwdRap_Sigma[i] * TMath::Hypot(LHCb_Run2_FwdRap_dSigmady_SysErrLow[i] / LHCb_Run2_FwdRap_dSigmady[i],
																							 LHCb_Run2_FwdRap_Flux_Err[i] / LHCb_Run2_FwdRap_Flux[i])		);
		LHCb_Run2_FwdRap_Sigma_SysErrHig.push_back( LHCb_Run2_FwdRap_Sigma[i] * TMath::Hypot(LHCb_Run2_FwdRap_dSigmady_SysErrHig[i] / LHCb_Run2_FwdRap_dSigmady[i],
																							 LHCb_Run2_FwdRap_Flux_Err[i] / LHCb_Run2_FwdRap_Flux[i])		);
		
		LHCb_Run2_FwdRap_x        		.push_back( Converter::y2x(LHCb_Run2_FwdRap_y[i]) );
		LHCb_Run2_FwdRap_R 			.push_back( sqrt(LHCb_Run2_FwdRap_Sigma[i]/LHCb_Run2_FwdRap_Sigma_IA[i])	);
		LHCb_Run2_FwdRap_R_StatErr		.push_back( 0.5 * LHCb_Run2_FwdRap_R[i] * LHCb_Run2_FwdRap_Sigma_StatErr[i]/LHCb_Run2_FwdRap_Sigma[i]	);
		LHCb_Run2_FwdRap_R_SysErrLow	.push_back( 0.5 * LHCb_Run2_FwdRap_R[i] * TMath::Hypot( LHCb_Run2_FwdRap_Sigma_SysErrLow[i]/LHCb_Run2_FwdRap_Sigma[i],
																								LHCb_Run2_FwdRap_Sigma_IA_Err[i]/LHCb_Run2_FwdRap_Sigma_IA[i])	);
		LHCb_Run2_FwdRap_R_SysErrHig	.push_back( 0.5 * LHCb_Run2_FwdRap_R[i] * TMath::Hypot( LHCb_Run2_FwdRap_Sigma_SysErrHig[i]/LHCb_Run2_FwdRap_Sigma[i],
																								LHCb_Run2_FwdRap_Sigma_IA_Err[i]/LHCb_Run2_FwdRap_Sigma_IA[i])	);
	}
	// std::vector<double> ALICE_Run2_MidRap_w_eff;
	// ALICE_Run2_MidRap_w_eff        	.push_back( Converter::y2w_eff(  0  ) );

	auto c = new TCanvas();
	// c->SetLogx();

	//gPad->SetPad(0.0,0.25,1.0,0.96);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.12);
    //gPad->SetLeftMargin(0.05);
    gPad->SetRightMargin(0.05);
	
	TH2D* SigmaVsW;
	auto X_AXIS_ERR 			= std::vector<double>(ShadowRatio_ParamsMap.at("Ws").size(), 3.8);
	auto X_AXIS_ERR_ALICE_Fwd 	= std::vector<double>(ShadowRatio_ParamsMap.at("Ws").size(), 3.8);
	auto X_AXIS_ERR_ALICE_Mid 	= std::vector<double>(ShadowRatio_ParamsMap.at("Ws").size(), 3.8);
	auto X_AXIS_ERR_LHCb	 	= std::vector<double>(ShadowRatio_ParamsMap.at("Ws").size(), 3.8);
	std::map<TString, std::vector<double>> IA_Map = {{"Ws",		{}}};	// For plotting IA
	for (int i = 5; i < 500; ++i){	IA_Map.at("Ws").push_back(i);	}
	getImpulseApprox(IA_Map);

	if(flag4Axis==0) //logY only
	{
		c->SetLogy();
		SigmaVsW = new TH2D("SigmaVsW", ";W_{#gammaN}^{Pb} (GeV);#sigma(#gamma Pb #rightarrow J/#psi Pb) (mb);", 10,0,420, 10, 0.004, 0.19);
		//SigmaVsW = new TH2D("SigmaVsW", ";W_{#gammaN}^{Pb} (GeV);#sigma(#gamma Pb #rightarrow J/#psi Pb) (mb);", 10,0,420, 10, 0.004, 0.13);
		// SigmaVsW = new TH2D("SigmaVsW", ";#omega_{#gamma}(in target Pb rest frame) (GeV);#sigma(#gamma Pb #rightarrow J/#psi Pb) (mb);", 10,0,3e8, 10, 0.004, 0.13);
	}
	else //logX and logY as default one
	{
		c->SetLogx();
		c->SetLogy();
		SigmaVsW = new TH2D("SigmaVsW", ";W_{#gammaN}^{Pb} (GeV);#sigma(#gamma A #rightarrow J/#psi A) (mb);", 10,14,520, 10, 0.004, 0.13);
		// SigmaVsW = new TH2D("SigmaVsW", ";#omega_{#gamma}(in target Pb rest frame) (GeV);#sigma(#gamma Pb #rightarrow J/#psi Pb) (mb);", 10,1e5,3e8, 10, 0.004, 0.13);
		X_AXIS_ERR 			= {11,2,12,1.7,14,1.5};
		X_AXIS_ERR_ALICE_Fwd= {0.65,0.8};
		X_AXIS_ERR_ALICE_Mid= {5};
		X_AXIS_ERR_LHCb		= {0.55,0.6};
	}

	SigmaVsW->GetYaxis()->CenterTitle();
    SigmaVsW->GetXaxis()->CenterTitle();
	SigmaVsW->GetYaxis()->SetTitleSize(0.05);
	SigmaVsW->GetYaxis()->SetTitleSize(0.05);
	SigmaVsW->GetYaxis()->SetTitleOffset(0.99);
	SigmaVsW->GetYaxis()->SetLabelSize(0.04);
	SigmaVsW->GetXaxis()->SetTitleSize(0.05);
	SigmaVsW->GetXaxis()->SetTitleOffset(0.98);
	SigmaVsW->GetXaxis()->SetLabelSize(0.04);
	SigmaVsW->SetTickLength(0.04);
	SigmaVsW->Draw("");

	auto Sigmas_TotalSysErr = TotalSysUncer_Map.at("Sigmas_TotalSysUncer");
	for (int i = 0; i < Sigmas_TotalSysErr.size(); ++i)
	{
		Sigmas_TotalSysErr[i] *= 0.01 * TotalSysUncer_Map.at("Sigmas")[i];
	}

	// std::vector<double> CMS_w_eff = {};
	// for (int i = 0; i < ShadowRatio_ParamsMap.at("Raps").size(); i++)
	// {
	// 	CMS_w_eff.push_back(Converter::y2w_eff(ShadowRatio_ParamsMap.at("Raps")[i]));
	// }

	TGraphErrors* ge_CMS        			 = new TGraphErrors(ShadowRatio_ParamsMap.at("Ws").size(),	&ShadowRatio_ParamsMap.at("Ws")[0],	&ShadowRatio_ParamsMap.at("Sigmas")[0],	0,	&ShadowRatio_ParamsMap.at("Sigmas_Err")[0]	);
	TGraphAsymmErrors* gae_CMS 	= new TGraphAsymmErrors(
												ShadowRatio_ParamsMap.at("Ws").size(),
												&ShadowRatio_ParamsMap.at("Ws")[0],        &ShadowRatio_ParamsMap.at("Sigmas")[0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&Sigmas_TotalSysErr[0], &Sigmas_TotalSysErr[0]);
	
	TGraphErrors* ge_ALICE_Run2_MidRap	= new TGraphErrors(ALICE_Run2_MidRap_W.size(),	&ALICE_Run2_MidRap_W[0],	&ALICE_Run2_MidRap_Sigma[0],	0,	&ALICE_Run2_MidRap_Sigma_StatErr[0]	);
	TGraphErrors* ge_ALICE_Run2_FwdRap	= new TGraphErrors(ALICE_Run2_FwdRap_W.size(),	&ALICE_Run2_FwdRap_W[0],	&ALICE_Run2_FwdRap_Sigma[0],	0,	&ALICE_Run2_FwdRap_Sigma_StatErr[0]	);
	TGraphAsymmErrors* gae_ALICE_Run2_MidRap 	= new TGraphAsymmErrors(
												ALICE_Run2_MidRap_W.size(),
												&ALICE_Run2_MidRap_W[0],        &ALICE_Run2_MidRap_Sigma[0],
												&X_AXIS_ERR_ALICE_Mid[0],  &X_AXIS_ERR_ALICE_Mid[0],
												&ALICE_Run2_MidRap_Sigma_SysErrLow[0], &ALICE_Run2_MidRap_Sigma_SysErrHig[0]);
	TGraphAsymmErrors* gae_ALICE_Run2_FwdRap 	= new TGraphAsymmErrors(
												ALICE_Run2_FwdRap_W.size(),
												&ALICE_Run2_FwdRap_W[0],        &ALICE_Run2_FwdRap_Sigma[0],
												&X_AXIS_ERR_ALICE_Fwd[0],  &X_AXIS_ERR_ALICE_Fwd[0],
												&ALICE_Run2_FwdRap_Sigma_SysErrLow[0], &ALICE_Run2_FwdRap_Sigma_SysErrHig[0]);
	TGraphErrors* ge_LHCb_Run2_FwdRap	= new TGraphErrors(LHCb_Run2_FwdRap_W.size(),	&LHCb_Run2_FwdRap_W[0],	&LHCb_Run2_FwdRap_Sigma[0],	0,	&LHCb_Run2_FwdRap_Sigma_StatErr[0]	);
	TGraphAsymmErrors* gae_LHCb_Run2_FwdRap 	= new TGraphAsymmErrors(
												LHCb_Run2_FwdRap_W.size(),
												&LHCb_Run2_FwdRap_W[0],        &LHCb_Run2_FwdRap_Sigma[0],
												&X_AXIS_ERR_LHCb[0],  &X_AXIS_ERR_LHCb[0],
												&LHCb_Run2_FwdRap_Sigma_SysErrLow[0], &LHCb_Run2_FwdRap_Sigma_SysErrHig[0]);

	TGraphErrors* ge_CGCnoFluct = new TGraphErrors(CGC_JpsiNoFluct_W.size(),	&CGC_JpsiNoFluct_W[0],	&CGC_JpsiNoFluct_CohXsec[0],	0,	0);
	TGraphErrors* ge_IA        	= new TGraphErrors(IA_Map.at("Ws").size(),	&IA_Map.at("Ws")[0],	&IA_Map["Sigmas_IA"][0],	0,	0);

	

	//drawLatex(0.15, 0.86, "Pb+Pb #rightarrow Pb+Pb+J/#psi #sqrt{s_{NN}} = 5.02 TeV",  42,        0.05,      1 );

	if(flag4Axis == 1)
	{
		TLegend  *legData =  new TLegend(0.14, 0.66, 0.37, 0.82);
		legData->SetFillStyle(0);
		legData->SetFillColor(0);
		legData->SetTextSize(0.04);
		legData->AddEntry(ge_CMS,               "CMS",               "p");
		legData->AddEntry(ge_ALICE_Run2_FwdRap, "ALICE* (-4 < y < -3.5)", "p");
		legData->AddEntry(ge_ALICE_Run2_MidRap, "ALICE* (|y| < 0.15)",  "p");
		legData->AddEntry(ge_LHCb_Run2_FwdRap,  "LHCb* (-4.5 < y < -3.5)", "p");
		legData->Draw("same");

		TLegend  *legTheory =  new TLegend(0.68, 0.15, 0.98, 0.50);
		legTheory->SetFillStyle(0);
		legTheory->SetFillColor(0);
		legTheory->SetTextSize(0.035);
		legTheory->AddEntry(ge_IA,                "Impulse Approx.", "l");
		legTheory->AddEntry(ge_CGCnoFluct,        "CGC IPsat",                   "l");

		drawGG("Sigmas",legTheory);

		drawLTA_Sigmas_R("Sigmas", legTheory);

		drawbBK("Sigmas",legTheory);

		drawCD("Sigmas",legTheory);

		legTheory->Draw("same");

		// TF1 *f1 	= new TF1("f1","log10(x)",3.5469263e-5,0.048933106);
		// TGaxis *A1 	= new TGaxis(14, 0.10,520, 0.10,"f1");
		// A1->SetTitleOffset(1);
		// A1->SetTitle("X");
		// A1->SetTicks("+");
		// A1->SetTitleSize(0.05);
		// A1->SetLabelSize(0.01);
		// A1->Draw("sames");

		// TH2D *hAxis = new TH2D("hAxis","hAxis",10,14,520, 10, 0.006, 0.10);
		// hAxis->GetYaxis()->SetLabelSize(0.2);
		// hAxis->SetLineColor(kRed);
		// hAxis->GetYaxis()->SetTickSize(0.5);
		// hAxis->Draw("Same X+");

        drawLatex(0.15, 0.84, "Pb+Pb #rightarrow Pb'+Pb'+J/#psi",  42,        0.05,      1 );
        drawLatex(0.56, 0.94, "PbPb 1.52 nb^{-1} (5.02 TeV)",      42,        0.05,      1 );
        drawLatex(0.15, 0.94, "#bf{CMS} #it{Preliminary}",                          42,        0.05,      1 );
	}
	else
	{
		TLegend  *legData =  new TLegend(0.20, 0.17, 0.50, 0.37);
		legData->SetFillStyle(0);
		legData->SetFillColor(0);
		legData->SetTextSize(0.04);
		legData->AddEntry(ge_CMS,               "CMS",                    "p");
		legData->AddEntry(ge_ALICE_Run2_FwdRap, "ALICE* (-4 < y < -3.5)", "p");
		legData->AddEntry(ge_ALICE_Run2_MidRap, "ALICE* (|y| < 0.15)",    "p");
		legData->AddEntry(ge_LHCb_Run2_FwdRap,  "LHCb* (-4.5 < y < -3.5)","p");
		legData->Draw("same");

		TLegend  *legTheory =  new TLegend(0.58, 0.31, 0.85, 0.46);
		legTheory->SetFillStyle(0);
		legTheory->SetFillColor(0);
		legTheory->SetTextSize(0.033);
		legTheory->AddEntry(ge_IA,                "Impulse Approx.", "l");
		legTheory->AddEntry(ge_CGCnoFluct,        "CGC IPsat",       "l");
		drawGG("Sigmas",legTheory);
		legTheory->Draw("same");
		

		TLegend  *leg_LTA =  new TLegend(0.58, 0.26, 0.95, 0.31);
		leg_LTA->SetNColumns(2);
		leg_LTA->SetFillStyle(0);
		leg_LTA->SetFillColor(0);
		leg_LTA->SetTextSize(0.030);
		drawLTA_Sigmas_R("Sigmas", leg_LTA);
		leg_LTA->Draw("same");
	
		TLegend  *leg_bBK =  new TLegend(0.58, 0.21, 0.91, 0.26);
		leg_bBK->SetNColumns(2);
		leg_bBK->SetFillStyle(0);
		leg_bBK->SetFillColor(0);
		leg_bBK->SetTextSize(0.030);
		drawbBK("Sigmas", leg_bBK);
		leg_bBK->Draw("same");
		
		TLegend  *leg_CD =  new TLegend(0.58, 0.16, 0.94, 0.21);
		leg_CD->SetNColumns(3);
		leg_CD->SetFillStyle(0);
		leg_CD->SetFillColor(0);
		leg_CD->SetTextSize(0.030);
		drawCD("Sigmas", leg_CD);
		leg_CD->Draw("same");	

		drawLatex(0.18,0.86,"x~3.8#times10^{-3}",42, 0.025, 1);
		drawLatex(0.29,0.86,"x~9.6#times10^{-4}",42, 0.025, 1);
		drawLatex(0.48,0.86,"x~2.4#times10^{-4}",42, 0.025, 1);
		drawLatex(0.67,0.86,"x~1.1#times10^{-4}",42, 0.025, 1);
		drawLatex(0.87,0.86,"x~6.0#times10^{-5}",42, 0.025, 1);

		drawLatex(0.15, 0.80, "Pb+Pb #rightarrow Pb'+Pb'+J/#psi",  42,        0.05,      1 );
        drawLatex(0.56,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
        drawLatex(0.15,0.94,"#bf{CMS} #it{Preliminary}",42, 0.05, 1);
	}
	
	ge_CGCnoFluct->SetMarkerColor(1);
	ge_CGCnoFluct->SetLineColor(11);
	ge_CGCnoFluct->SetLineStyle(5);
	ge_CGCnoFluct->SetLineWidth(2);
	ge_CGCnoFluct->Draw("lsame");

	ge_IA->SetMarkerColor(2);
	ge_IA->SetLineStyle(5);
	ge_IA->SetLineColor(1);
	ge_IA->SetLineWidth(2);
	ge_IA->Draw("lsame");

	gae_CMS ->SetMarkerStyle(24);
	gae_CMS ->SetFillColorAlpha(16, 0.7);
	gae_CMS ->SetFillStyle(1001);
	gae_CMS ->Draw("2same");
	ge_CMS->SetMarkerStyle(20);
	ge_CMS->SetMarkerSize(1.6);
	ge_CMS->SetMarkerColor(2);
	ge_CMS->SetLineColor(2);
	ge_CMS->SetLineWidth(2);
	ge_CMS->Draw("pezsame");

	gae_ALICE_Run2_MidRap ->SetMarkerStyle(24);
	gae_ALICE_Run2_MidRap ->SetFillColorAlpha(16, 0.7);
	gae_ALICE_Run2_MidRap ->SetFillStyle(1001);
	gae_ALICE_Run2_MidRap ->Draw("2same");
	ge_ALICE_Run2_MidRap->SetMarkerStyle(24);
	ge_ALICE_Run2_MidRap->SetMarkerColor(4);
    ge_ALICE_Run2_MidRap->SetMarkerSize(1.6);
	ge_ALICE_Run2_MidRap->SetLineColor(4);
	ge_ALICE_Run2_MidRap->SetLineWidth(2);
	ge_ALICE_Run2_MidRap->Draw("pezsame");

	gae_ALICE_Run2_FwdRap ->SetMarkerStyle(24);
	gae_ALICE_Run2_FwdRap ->SetFillColorAlpha(16, 0.7);
	gae_ALICE_Run2_FwdRap ->SetFillStyle(1001);
	gae_ALICE_Run2_FwdRap ->Draw("2same");
	ge_ALICE_Run2_FwdRap->SetMarkerStyle(25);
	ge_ALICE_Run2_FwdRap->SetMarkerColor(4);
    ge_ALICE_Run2_FwdRap->SetMarkerSize(1.5);
	ge_ALICE_Run2_FwdRap->SetLineColor(4);
	ge_ALICE_Run2_FwdRap->SetLineWidth(2);
	ge_ALICE_Run2_FwdRap->Draw("pezsame");

	gae_LHCb_Run2_FwdRap ->SetMarkerStyle(24);
	gae_LHCb_Run2_FwdRap ->SetFillColorAlpha(16, 0.7);
	gae_LHCb_Run2_FwdRap ->SetFillStyle(1001);
	gae_LHCb_Run2_FwdRap ->Draw("2same");
	ge_LHCb_Run2_FwdRap->SetMarkerStyle(26);
	ge_LHCb_Run2_FwdRap->SetMarkerColor(4);
    ge_LHCb_Run2_FwdRap->SetMarkerSize(1.5);
	ge_LHCb_Run2_FwdRap->SetLineColor(4);
	ge_LHCb_Run2_FwdRap->SetLineWidth(2);
	ge_LHCb_Run2_FwdRap->Draw("pezsame");

	if (flag4Axis)
	{
		c->SaveAs("outplots/SigmaVsW_logXY.png");
		c->SaveAs("outplots/SigmaVsW_logXY.pdf");
	}
	else
	{
		c->SaveAs("outplots/SigmaVsW_logY.png");
		c->SaveAs("outplots/SigmaVsW_logY.pdf");
	}

//	delete c;
}

void plotRvsX( std::map<TString, std::vector<double>> &ShadowRatio_ParamsMap, const std::map<TString, std::vector<double>> &TotalSysUncer_Map )
{	
	auto c = new TCanvas();

	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.12);
	//gPad->SetLeftMargin(0.05);
	gPad->SetRightMargin(0.05);


	c->SetLogx();

	//TH2D* htem2d = new TH2D("htem2d", "", 10,3.0e-5,5.0e-2, 10, 0, 1);
	TH2D* htem2d = new TH2D("htem2d", "", 10,1.0e-5,5e-2, 10, 0.2, 1.05);
	htem2d->SetYTitle("R^{Pb}_{g}(x, #mu^{2}=2.4 GeV^{2})");
	htem2d->SetXTitle("x");

	htem2d->GetYaxis()->SetTitleSize(0.06);
	htem2d->GetYaxis()->SetTitleOffset(0.85);
	htem2d->GetYaxis()->SetLabelSize(0.04);
	htem2d->GetYaxis()->CenterTitle();
    htem2d->GetXaxis()->CenterTitle();
	htem2d->GetXaxis()->SetTitleSize(0.07);
	htem2d->GetXaxis()->SetTitleOffset(0.69);
	htem2d->GetXaxis()->SetLabelSize(0.04);
	htem2d->SetTickLength(0.04);
	htem2d->Draw();

	// std::vector<double> X_AXIS_ERR = {5e-6,1.5e-4,4e-6,2e-4,3e-6,2.5e-4};
	std::vector<double> X_AXIS_ERR = {8e-6,3e-4,7e-6,4e-4,5e-6,4.5e-4};
	auto R_TotalSysErr = TotalSysUncer_Map.at("R_TotalSysUncer");
	for (int i = 0; i < R_TotalSysErr.size(); ++i)
	{
		R_TotalSysErr[i] *= 0.01 * TotalSysUncer_Map.at("R")[i];
	}
	TGraphErrors* ge_CMS        	= new TGraphErrors(ShadowRatio_ParamsMap.at("Xs").size(),    &ShadowRatio_ParamsMap.at("Xs")[0], 
					&ShadowRatio_ParamsMap.at("R")[0],	&ShadowRatio_ParamsMap["Xs_Err"][0],	&ShadowRatio_ParamsMap.at("R_Err")[0]);
	TGraphAsymmErrors* gae_CMS 	= new TGraphAsymmErrors(
												ShadowRatio_ParamsMap.at("Xs").size(),
												&ShadowRatio_ParamsMap.at("Xs")[0],        &ShadowRatio_ParamsMap.at("R")[0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&R_TotalSysErr[0], &R_TotalSysErr[0]);

	TGraphErrors* ge_ALICE_Run2_MidRap	= new TGraphErrors(ALICE_Run2_MidRap_x.size(),	&ALICE_Run2_MidRap_x[0],	&ALICE_Run2_MidRap_R[0],	0,	&ALICE_Run2_MidRap_R_StatErr[0]	);
	TGraphErrors* ge_ALICE_Run2_FwdRap	= new TGraphErrors(ALICE_Run2_FwdRap_x.size(),	&ALICE_Run2_FwdRap_x[0],	&ALICE_Run2_FwdRap_R[0],	0,	&ALICE_Run2_FwdRap_R_StatErr[0]	);
	TGraphErrors* ge_LHCb_Run2_FwdRap	= new TGraphErrors(LHCb_Run2_FwdRap_x.size(),	&LHCb_Run2_FwdRap_x[0],	&LHCb_Run2_FwdRap_R[0],	0,	&LHCb_Run2_FwdRap_R_StatErr[0]	);
	X_AXIS_ERR = {5e-5};
	TGraphAsymmErrors* gae_ALICE_Run2_MidRap 	= new TGraphAsymmErrors(
												ALICE_Run2_MidRap_x.size(),
												&ALICE_Run2_MidRap_x[0],        &ALICE_Run2_MidRap_R[0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&ALICE_Run2_MidRap_R_SysErrLow[0], &ALICE_Run2_MidRap_R_SysErrHig[0]);
	X_AXIS_ERR = {2.5e-3,2e-3};
	TGraphAsymmErrors* gae_ALICE_Run2_FwdRap 	= new TGraphAsymmErrors(
												ALICE_Run2_FwdRap_x.size(),
												&ALICE_Run2_FwdRap_x[0],        &ALICE_Run2_FwdRap_R[0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&ALICE_Run2_FwdRap_R_SysErrLow[0], &ALICE_Run2_FwdRap_R_SysErrHig[0]);
	TGraphAsymmErrors* gae_LHCb_Run2_FwdRap 	= new TGraphAsymmErrors(
												LHCb_Run2_FwdRap_x.size(),
												&LHCb_Run2_FwdRap_x[0],        &LHCb_Run2_FwdRap_R[0],
												&X_AXIS_ERR[0],  &X_AXIS_ERR[0],
												&LHCb_Run2_FwdRap_R_SysErrLow[0], &LHCb_Run2_FwdRap_R_SysErrHig[0]);
	gae_CMS ->SetMarkerStyle(24);
	gae_CMS ->SetFillColorAlpha(16, 0.7);
	gae_CMS ->SetFillStyle(1001);
	ge_CMS->SetMarkerStyle(20);
	ge_CMS->SetMarkerColor(2);
    ge_CMS->SetLineColor(2);
    ge_CMS->SetMarkerSize(1.6);
	ge_CMS->SetLineWidth(2);

	gae_ALICE_Run2_MidRap ->SetMarkerStyle(24);
	gae_ALICE_Run2_MidRap ->SetFillColorAlpha(16, 0.7);
	gae_ALICE_Run2_MidRap ->SetFillStyle(1001);
	ge_ALICE_Run2_MidRap->SetMarkerStyle(24);
	ge_ALICE_Run2_MidRap->SetMarkerColor(4);
    ge_ALICE_Run2_MidRap->SetMarkerSize(1.6);
	ge_ALICE_Run2_MidRap->SetLineColor(4);
	ge_ALICE_Run2_MidRap->SetLineWidth(2);

	gae_ALICE_Run2_FwdRap ->SetMarkerStyle(24);
	gae_ALICE_Run2_FwdRap ->SetFillColorAlpha(16, 0.7);
	gae_ALICE_Run2_FwdRap ->SetFillStyle(1001);
	ge_ALICE_Run2_FwdRap->SetMarkerStyle(25);
	ge_ALICE_Run2_FwdRap->SetMarkerColor(4);
    ge_ALICE_Run2_FwdRap->SetMarkerSize(1.5);
	ge_ALICE_Run2_FwdRap->SetLineColor(4);
	ge_ALICE_Run2_FwdRap->SetLineWidth(2);

	gae_LHCb_Run2_FwdRap ->SetMarkerStyle(24);
	gae_LHCb_Run2_FwdRap ->SetFillColorAlpha(16, 0.7);
	gae_LHCb_Run2_FwdRap ->SetFillStyle(1001);
	ge_LHCb_Run2_FwdRap->SetMarkerStyle(26);
	ge_LHCb_Run2_FwdRap->SetMarkerColor(4);
    ge_LHCb_Run2_FwdRap->SetMarkerSize(1.5);
	ge_LHCb_Run2_FwdRap->SetLineColor(4);
	ge_LHCb_Run2_FwdRap->SetLineWidth(2);

	TLegend  *leg =  new TLegend(0.13, 0.62, 0.40, 0.8);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.04);
	leg->AddEntry(ge_CMS,               "CMS",                     "p");
	leg->AddEntry(ge_ALICE_Run2_FwdRap, "ALICE* (-4 < y < -3.5)",  "p");
	leg->AddEntry(ge_ALICE_Run2_MidRap, "ALICE* (|y| < 0.15)",     "p");
	leg->AddEntry(ge_LHCb_Run2_FwdRap,  "LHCb* (-4.5 < y < -3.5)", "p");
    leg->Draw("same");

    TLegend  *leg1 =  new TLegend(0.53, 0.35, 0.95, 0.40);
	leg1->SetNColumns(2);
    leg1->SetFillStyle(0);
    leg1->SetFillColor(0);
    leg1->SetTextSize(0.033);
	drawGG("R",leg1);
	leg1->Draw("same");
	//drawLTA_Sigmas_R("R", leg1);
	//drawbBK("R", leg1);
	//drawCD("R", leg1);

	TLegend  *leg_LTA =  new TLegend(0.53, 0.30, 0.95, 0.35);
	leg_LTA->SetNColumns(2);
	leg_LTA->SetFillStyle(0);
	leg_LTA->SetFillColor(0);
	leg_LTA->SetTextSize(0.033);
	drawLTA_Sigmas_R("R", leg_LTA);
	leg_LTA->Draw("same");

	TLegend  *leg_bBK =  new TLegend(0.53, 0.25, 0.91, 0.30);
	leg_bBK->SetNColumns(2);
	leg_bBK->SetFillStyle(0);
	leg_bBK->SetFillColor(0);
	leg_bBK->SetTextSize(0.033);
	drawbBK("R", leg_bBK);
	leg_bBK->Draw("same");

	TLegend  *leg_CD =  new TLegend(0.53, 0.20, 0.94, 0.25);
	leg_CD->SetNColumns(3);
	leg_CD->SetFillStyle(0);
	leg_CD->SetFillColor(0);
	leg_CD->SetTextSize(0.033);
	drawCD("R", leg_CD);
	leg_CD->Draw("same");
	

	gae_CMS ->Draw("2same");
	ge_CMS->Draw("pezsame");
	gae_ALICE_Run2_MidRap ->Draw("2same");
	ge_ALICE_Run2_MidRap->Draw("pezsame");
	gae_ALICE_Run2_FwdRap ->Draw("2same");
	ge_ALICE_Run2_FwdRap->Draw("pezsame");
	gae_LHCb_Run2_FwdRap ->Draw("2same");
	ge_LHCb_Run2_FwdRap->Draw("pezsame");

	drawLatex(0.15, 0.84, "Pb+Pb #rightarrow Pb'+Pb'+J/#psi",  42,        0.05,      1 );
	//drawLatex(0.48, 0.17, "Nuclear suppression factor", 42, 0.055, 1);

    drawLatex(0.56,0.94,"PbPb 1.52 nb^{-1} (5.02 TeV)",42, 0.05, 1);
    drawLatex(0.15,0.94,"#bf{CMS} #it{Preliminary}",42, 0.05, 1);

	c->SaveAs("outplots/ShadowingRatiovsX.png");
	c->SaveAs("outplots/ShadowingRatiovsX.pdf");

//	delete c;

	return ShadowRatio_ParamsMap;
}

std::map<TString, std::vector<double>> getParamsMap(const TString infile = "../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.appliedTnP.root", const double flux_uncer = 1.0, const TString fluxSubCase = "_SigNN68p3R6p67a0p56")
{
	//---------------------------------Read Map------------------------------------------
	std::map<TString, std::vector<double>> ShadowRatio_ParamsMap = {
		{"Raps",		{}},	{"Raps_Err",		{}},
		{"Sigmas",		{}},	{"Sigmas_Err",		{}},
		{"Sigmas_IA",	{}},	{"Sigmas_IA_Err",	{}},
		{"R",			{}},	{"R_Err",			{}},
		{"Xs",			{}},	{"Xs_Err",			{}},
		{"Ws",			{}},	{"Ws_Err",			{}},

		// {"Rap",			{}},	{"RapErr",				{}},
		// {"Xsec_AnAn",	{}},	{"XsecErr_AnAn",		{}},
		// {"Xsec_0n0n",	{}},	{"XsecErr_0n0n",		{}},
		// {"Xsec_0nXnSum",{}},	{"XsecErr_0nXnSum",		{}},
		// {"Xsec_XnXn",	{}},	{"XsecErr_XnXn",		{}},
	};

	LoadJpsiXsec loadJpsiXsec(infile);

	{
		auto temp_map = loadJpsiXsec.GetMap();
		ShadowRatio_ParamsMap.insert(temp_map.begin(),	temp_map.end());
	}

	for (int i = 0; i < ShadowRatio_ParamsMap.at("Rap").size(); ++i)
	{
		ShadowRatio_ParamsMap.at("Rap")[i] *= -1;
	}
	//-------------------------------------------------------------------------------------

	//---------------------------------Manual Test Map------------------------------------------
	//WS
	// std::map<TString, std::vector<double>> ShadowRatio_ParamsMap = {
	// 	{"Raps",		{}},	{"Raps_Err",		{}},
	// 	{"Sigmas",		{}},	{"Sigmas_Err",		{}},
	// 	{"Sigmas_IA",	{}},	{"Sigmas_IA_Err",	{}},
	// 	{"R",			{}},	{"R_Err",			{}},
	// 	{"Xs",			{}},	{"Xs_Err",			{}},
	// 	{"Ws",			{}},	{"Ws_Err",			{}},
	// 	{"Rap",			{-1.7000, -1.8000, -1.9000, -2.0000, -2.1000, -2.2000, -2.3000, -2.4000}},	{"RapErr",				{0,0,0,0,0,0,0,0}},
	// 	{"Xsec_AnAn",	{3.6673,3.5814,3.4918,3.4005,3.2986,3.2020,3.0970,2.9921}},	{"XsecErr_AnAn",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_0n0n",	{2.6450, 2.5837, 2.5214, 2.4603, 2.3931, 2.3332, 2.2687, 2.2067}},	{"XsecErr_0n0n",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_0nXnSum",{.7501,.7296,.7069,.6823,.6546,.6256,.5940,.5610}},	{"XsecErr_0nXnSum",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_XnXn",	{.2722,.2682,.2634,.2578,.2509,.2433,.2343,.2244}},	{"XsecErr_XnXn",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// };
	// std::map<TString, std::vector<double>> tem = {
	// 	{"Rap",			{1.7000, 1.8000, 1.9000, 2.0000, 2.1000, 2.2000, 2.3000, 2.4000}},	{"RapErr",				{0,0,0,0,0,0,0,0}},
	// 	{"Xsec_AnAn",	{3.6673,3.5814,3.4918,3.4005,3.2986,3.2020,3.0970,2.9921}},	{"XsecErr_AnAn",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_0n0n",	{2.6450, 2.5837, 2.5214, 2.4603, 2.3931, 2.3332, 2.2687, 2.2067}},	{"XsecErr_0n0n",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_0nXnSum",{.7501,.7296,.7069,.6823,.6546,.6256,.5940,.5610}},	{"XsecErr_0nXnSum",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_XnXn",	{.2722,.2682,.2634,.2578,.2509,.2433,.2343,.2244}},	{"XsecErr_XnXn",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// };

	//SS
	// std::map<TString, std::vector<double>> ShadowRatio_ParamsMap = {
	// 	{"Raps",		{}},	{"Raps_Err",		{}},
	// 	{"Sigmas",		{}},	{"Sigmas_Err",		{}},
	// 	{"Sigmas_IA",	{}},	{"Sigmas_IA_Err",	{}},
	// 	{"R",			{}},	{"R_Err",			{}},
	// 	{"Xs",			{}},	{"Xs_Err",			{}},
	// 	{"Ws",			{}},	{"Ws_Err",			{}},
	// 	{"Rap",			{-1.7000, -1.8000, -1.9000, -2.0000, -2.1000, -2.2000, -2.3000, -2.4000  }},	{"RapErr",				{0.,0.,0,0,0,0,0,0}},
	// 	{"Xsec_AnAn",	{2.5587,  2.5201,  2.4798,  2.4395,  2.3908,  2.3473,  2.2968,  2.2466}},	{"XsecErr_AnAn",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_0n0n",	{1.9093,	1.8860,	1.8621,	1.8395,	1.8112,	1.7887,	1.7614,	1.7356}},	{"XsecErr_0n0n",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_0nXnSum",{.4788,.4662,.4528,.4386,.4224,.4059,.3879,.3692}},	{"XsecErr_0nXnSum",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_XnXn",	{.1706,.1679,.1648,.1614,.1572,.1527,.1475,.1418}},	{"XsecErr_XnXn",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// };

	// std::map<TString, std::vector<double>> tem = {
	// 	{"Rap",			{1.7000, 1.8000, 1.9000, 2.0000, 2.1000, 2.2000, 2.3000, 2.4000  }},	{"RapErr",				{0.,0.,0,0,0,0,0,0}},
	// 	{"Xsec_AnAn",	{2.5587,  2.5201,  2.4798,  2.4395,  2.3908,  2.3473,  2.2968,  2.2466}},	{"XsecErr_AnAn",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_0n0n",	{1.9093,	1.8860,	1.8621,	1.8395,	1.8112,	1.7887,	1.7614,	1.7356}},	{"XsecErr_0n0n",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_0nXnSum",{.4788,.4662,.4528,.4386,.4224,.4059,.3879,.3692}},	{"XsecErr_0nXnSum",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// 	{"Xsec_XnXn",	{.1706,.1679,.1648,.1614,.1572,.1527,.1475,.1418}},	{"XsecErr_XnXn",		{0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}},
	// };
	//-------------------------------------------------------------------------------------

	{auto Temp = fit2D(loadJpsiXsec.GetMap(), flux_uncer, fluxSubCase); ShadowRatio_ParamsMap.at("Raps") = Temp[0]; ShadowRatio_ParamsMap["Raps_Err"] = Temp[1]; ShadowRatio_ParamsMap.at("Sigmas") = Temp[2]; ShadowRatio_ParamsMap.at("Sigmas_Err") = Temp[3];}
	for (int i = 0; i < ShadowRatio_ParamsMap.at("Raps").size(); ++i)
	{
		ShadowRatio_ParamsMap.at("Xs") 		.push_back( Converter::y2x(ShadowRatio_ParamsMap.at("Raps")[i]) 	);
		ShadowRatio_ParamsMap["Xs_Err"]		.push_back( 0 );
		ShadowRatio_ParamsMap.at("Ws") 		.push_back( Converter::y2W(ShadowRatio_ParamsMap.at("Raps")[i]) 	);
		ShadowRatio_ParamsMap["Ws_Err"] 	.push_back( 0 );

		cout<<Form("y: %f,	x: %f,	W: %f", ShadowRatio_ParamsMap.at("Raps")[i], ShadowRatio_ParamsMap.at("Xs")[i], ShadowRatio_ParamsMap.at("Ws")[i])<<endl;
		cout<<"x from W: "<<pow(mJpsi_PDG,2)/pow(ShadowRatio_ParamsMap.at("Ws")[i],2)<<endl;
	}
	
	{getImpulseApprox(ShadowRatio_ParamsMap);}

	for (int i = 0; i < ShadowRatio_ParamsMap.at("Raps").size(); ++i)
	{
		ShadowRatio_ParamsMap.at("R").push_back( sqrt( ShadowRatio_ParamsMap.at("Sigmas")[i]/ShadowRatio_ParamsMap.at("Sigmas_IA")[i] ) );
	}
	
	Cal_R_Error(ShadowRatio_ParamsMap);

	return ShadowRatio_ParamsMap;
}

void plotShadowingRatio()
{
	const int template_option  = 1; //0:OldCohJpsi; 1:NewCohJpsi (w R+1fm);
	const TString template_Name[2] = {"", "_NewCohJpsi"};

	const int   RunTnPcase         	= 1;	//Default 1
	const TString TnPcases[4]		= {"", ".appliedTnP", ".appliedTnP_Low", ".appliedTnP_Hig"};

	const int   RunHFcase         	= 0;	//Default 0
	const TString HFcases[4]		= {"", ".looseHF", ".tightHF", ".removeHF"};

	//---------------------------------Remake Map------------------------------------------
	// auto ShadowRatio_ParamsMap = getParamsMap(Form("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins%s%s%s.root",template_Name[template_option].Data(), TnPcases[RunTnPcase].Data(), HFcases[RunHFcase].Data()));
	// auto TotalSysUncer_Map 		= readMap(Form("rootfiles/TotalSysUncer%s_Map.root", template_Name[template_option].Data()));
	// plotSigmaVsW(ShadowRatio_ParamsMap,TotalSysUncer_Map);
	// plotRvsX(ShadowRatio_ParamsMap,TotalSysUncer_Map);
	// saveMap(ShadowRatio_ParamsMap, Form("rootfiles/Results%s_Map.root", template_Name[template_option].Data()));
	//-------------------------------------------------------------------------------------
	
	//----------------------------Read Map From Root File----------------------------------
	auto ShadowRatio_ParamsMap 	= readMap(Form("rootfiles/Results%s_Map.root", template_Name[template_option].Data()));
	auto TotalSysUncer_Map 		= readMap(Form("rootfiles/TotalSysUncer%s_Map.root", template_Name[template_option].Data()));
	plotSigmaVsW(ShadowRatio_ParamsMap,TotalSysUncer_Map);
	plotRvsX(ShadowRatio_ParamsMap,TotalSysUncer_Map);
	//-------------------------------------------------------------------------------------
}
