/*
# Author: 
# Kong Tu, BNL
# Time: 
# March 30 2022

# Purpose:
# This macros is intended for calculating the nuclear Form Factor
# and its relation to connect the heavy nucleus gammaA->VM A cross
# section to gammaP->VM P cross section. This approach, without
# considering any nuclear effect, is known as the ``Impulse Approximation".
# When comparing with the data, the ``leading twist approximation of nuclear
# shadowing effect" can be experimentally extracted. 

# Useful references 
  (Guzey, Strikman, Zhalov, et al):
# Phys. Lett. B 726 (2013) 290-295, https://arxiv.org/abs/1305.1724
# JHEP 10 (2013) 207, https://arxiv.org/abs/1305.1724
# Eur.Phys.J.C 74 (2014) 7, 2942, https://arxiv.org/abs/1312.6486
# Nucleus distribution and properties, https://arxiv.org/pdf/1408.2549.pdf
*/

#include "RiceStyle.h"
using namespace std;
static const double hbarc    = 0.1973269718;
static const double hbarcmev = hbarc*1000.;
static const double pi       = 3.141592654;
static const double twoPi    = 2 * pi;
static const double alpha    = 1/137.035999074;

double formFactor(Double_t *t, Double_t *par)
{
  double A_nucleus = par[0];
  double _Radius = par[1];
  double _rho0 = 0.169551;
  // nuclear form factor, from Klein Nystrand PRC 60 (1999) 014903, Eq. 14
  const double R    = _Radius;
  const double q    = sqrt(t[0]);
  const double arg1 = q * R / hbarc;
  const double arg2 = hbarc / (q * R);
  const double sph  = (sin(arg1) - arg1 * cos(arg1)) * 3. * arg2 * arg2 * arg2;
  const double a0   = par[2];
  double ff = sph / (1. + (a0 * a0 * t[0]) / (hbarc * hbarc));
  ff = ff * A_nucleus;
  return ff*ff;
        
}
std::vector<std::vector<double>> runUPC_AAModel(const std::vector<double> Ws, TString nucleus="Pb"){

  //constants
  double Mvm=3.09;//mass jpsi
  double M_N= (0.93827+0.93957) / 2;

  //coherent (elastic) proton d\sigma/dt | t=0 
  //as a function of W
  TF1* giveGammaP = new TF1("giveGammaP","[0]*TMath::Power(1-TMath::Power(([1]+[2])/x,2),1.5)*TMath::Power((x*x)/(100*100),[3])",4,506.);
  giveGammaP->SetTitle("#frac{d#sigma}{dt}(W_{#gamma p}, t=0)");
  giveGammaP->SetParameter(0, 342.);
  giveGammaP->SetParameter(1, Mvm);
  giveGammaP->SetParameter(2, M_N);
  giveGammaP->SetParameter(3, 0.4);

  TCanvas* c1 = new TCanvas("c1","c1",1,1,900,600);
  giveGammaP->GetYaxis()->SetTitle("#frac{d#sigma}{dt}(W_{#gamma p}, t=0) (nb/GeV^{2})");
  giveGammaP->GetXaxis()->SetTitle("W_{#gamma p} (GeV)");
  giveGammaP->Draw();

  //nuclear Form Factor squared 
  TF1* nuclearFF = new TF1("nuclearFF", formFactor, 1e-8, 0.2, 2);
  if( nucleus == "Pb" ){
    nuclearFF->SetParameter(0,208);
    nuclearFF->SetParameter(1,6.62);
    nuclearFF->SetParameter(2,0.535);
  }
  else if( nucleus == "Au"){
    nuclearFF->SetParameter(0,197); 
    nuclearFF->SetParameter(1,6.38);
    nuclearFF->SetParameter(2,0.546);
  }
  else{
    throw std::runtime_error("runUPC_AAModel No nucleus Setting!!!");
  }
  
  //IA histogram as a function of W
  TH1D* h_IA = new TH1D("h_IA",";W_{#gamma p} (GeV); #sigma^{#gamma+A #rightarrow VM+A} (W_{#gamma p}) [nb]",500,4.5,505.5);
  TH1D* h_IA_= new TH1D(*h_IA);
  for(int ibin=0;ibin<h_IA->GetNbinsX();ibin++){
    double w = h_IA->GetBinCenter(ibin+1);
    double gammaP_xs = giveGammaP->Eval(w);
    double t_min = (TMath::Power(Mvm,4)*M_N*M_N) / TMath::Power(w,4);
    double Phi_A = nuclearFF->Integral(t_min,0.2);//integrate from t_min to infinity.

    //Eq. 12 in Eur.Phys.J.C 74 (2014) 7, 2942
    double IA = gammaP_xs * Phi_A;
    h_IA->SetBinContent(ibin+1, IA);
    h_IA->SetBinError(ibin+1, IA*0.02);//2% error for now
  }

  TCanvas* c2 = new TCanvas("c2","c2",1,1,900,600);
  h_IA->SetStats(kFALSE);
  h_IA->SetTitle("Impulse Approximation");
  h_IA->SetMarkerStyle(24);
  h_IA->Draw("PE");

  if(nucleus=="Au"){
    TBox* box1 = new TBox(17,10e3,36,30e3);
    box1->SetFillColorAlpha(kRed+2,0.2);
    box1->SetFillStyle(1001);
    box1->SetLineWidth(0);
    box1->Draw("same");
    box1->Draw("same");

    TLatex *latex1 = new TLatex(0.16, 0.15, "STAR (2016)");
    latex1->SetNDC();
    latex1->SetTextSize(23);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    latex1->Draw("same");

    TBox* box2 = new TBox(86,43e3,183,95e3);
    box2->SetFillColorAlpha(kBlue+2,0.2);
    box2->SetFillStyle(1001);
    box2->SetLineWidth(0);
    box2->Draw("same");
    box2->Draw("same");

    TLatex *latex2 = new TLatex(0.4, 0.3, "STAR (2022+)");
    latex2->SetNDC();
    latex2->SetTextSize(23);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    latex2->Draw("same");
  }
  else if(nucleus=="Pb"){
    TBox* box1 = new TBox(37,20e3,414,180e3);
    box1->SetFillColorAlpha(kRed+2,0.2);
    box1->SetFillStyle(1001);
    box1->SetLineWidth(0);
    box1->Draw("same");
    box1->Draw("same");

    TLatex *latex1 = new TLatex(0.16, 0.75, "CMS |y_{VM}| < 2.4 (Run 2 & 3)");
    latex1->SetNDC();
    latex1->SetTextSize(23);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    latex1->Draw("same");
  }
  



  // const std::vector<double> Ws = {41.503042, 50.691930, 306.668307, 374.565516};
  std::vector<double> IAs = {}, IAs_Err={};
  h_IA_->SetMarkerColor(kRed);
  for (int i = 0; i < Ws.size(); ++i)
  {
    int bin   = h_IA->FindBin( Ws[i] );
    double s  = h_IA->GetBinContent( bin );
    double s_err  = h_IA->GetBinError( bin );
    h_IA_->SetBinContent(bin, s);
    IAs.push_back(s/1e6);
    IAs_Err.push_back(s_err/1e6);
  }
  TCanvas* c3 = new TCanvas("c3","c3",1,1,900,600);
  THStack *hs = new THStack("hs", ";W_{#gamma p} (GeV); #sigma^{#gamma+A #rightarrow VM+A} (W_{#gamma p}) [nb]");
  hs->Add(h_IA); hs->Add(h_IA_);
  hs->Draw("NOSTACK P");

  for (int i = 0; i < Ws.size(); ++i)
  {
    TLatex *latex1 = new TLatex(0.33, 0.35-0.04*i, Form("W_{#gamma p} = %.2f GeV, #sigma^{#gamma+A #rightarrow VM+A} (W_{#gamma p}) = %f #pm %f mb",Ws[i],IAs[i],IAs_Err[i]));
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    latex1->Draw("same");
  }
  {
    TBox* box1 = new TBox(37,20e3,414,180e3);
    box1->SetFillColorAlpha(kRed+2,0.1);
    box1->SetFillStyle(1001);
    box1->SetLineWidth(0);
    box1->Draw("same");
    box1->Draw("same");

    TLatex *latex1 = new TLatex(0.16, 0.80, "CMS |y_{VM}| < 2.4 (Run 2 & 3)");
    latex1->SetNDC();
    latex1->SetTextSize(23);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    latex1->Draw("same");
  }

  c1->SaveAs("./ImpulseApproximation/runUPC_AAModel_0.png");
  c2->SaveAs("./ImpulseApproximation/runUPC_AAModel_1.png");
  c3->SaveAs("./ImpulseApproximation/runUPC_AAModel_2.png");
  delete c1;
  delete c2;
  delete c3;
  return {IAs, IAs_Err};
}