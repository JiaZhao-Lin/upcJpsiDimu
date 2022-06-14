double NuclearFormFactor(Double_t *t, Double_t *par)
{
  const double hbarc    = 0.1973269718;

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

void getImpulseApprox(std::map<TString, std::vector<double>> &ShadowRatio_ParamsMap)
{
  //constants
  const double JpsiMass   = 3.096916;
  const double M_N= (0.93827+0.93957) / 2;

  //coherent (elastic) proton d\sigma/dt | t=0 
  //as a function of W
  TF1* DSigmaDt = new TF1("DSigmaDt","[0]*TMath::Power(1-TMath::Power(([1]+[2])/x,2),1.5)*TMath::Power((x*x)/(100*100),[3])",4,4000); //original:4,506.
  DSigmaDt->SetTitle("#frac{d#sigma}{dt}(W_{#gamma p}, t=0)");
  DSigmaDt->SetParameter(0, 342.);
  DSigmaDt->SetParameter(1, JpsiMass);
  DSigmaDt->SetParameter(2, M_N);
  DSigmaDt->SetParameter(3, 0.4);

  DSigmaDt->GetYaxis()->SetTitle("#frac{d#sigma}{dt}(W_{#gamma p}, t=0) (nb/GeV^{2})");
  DSigmaDt->GetXaxis()->SetTitle("W_{#gamma p} (GeV)");

  //nuclear Form Factor squared 
  TF1* NuclearFF = new TF1("NuclearFF", NuclearFormFactor, 1e-8, 0.2, 3);

  NuclearFF->SetParameter(0,208);
  NuclearFF->SetParameter(1,6.62);
  NuclearFF->SetParameter(2,0.535);

  std::vector<double> IAs = {}, IAs_Err={};
  for (int i = 0; i < ShadowRatio_ParamsMap["Ws"].size(); ++i)
  {
    const double w = ShadowRatio_ParamsMap["Ws"][i];
    const double gammaP_xs = DSigmaDt->Eval(w);
    const double t_min = (TMath::Power(JpsiMass,4)*M_N*M_N) / TMath::Power(w,4);
    const double Phi_A = NuclearFF->Integral(t_min,0.2);//integrate from t_min to infinity. Original: 0.2

    //Eq. 12 in Eur.Phys.J.C 74 (2014) 7, 2942
    double IA = gammaP_xs * Phi_A;
    IA *= 1e-6;  // nb to mb 
    double IA_Err = IA * 0.054;

    IAs.push_back(IA);
    IAs_Err.push_back(IA_Err);
  }

  ShadowRatio_ParamsMap["Sigmas_IA"]        = IAs;
  ShadowRatio_ParamsMap["Sigmas_IA_Err"]    = IAs_Err;
}

void ImpulseApprox()
{
  //testing
  std::map<TString, std::vector<double>> TestMap = 
  {
    {"Ws",    {17.961384, 20.352915,124.67672,  200,300,400,500}}
  };

  getImpulseApprox(TestMap);
  for (int i = 0; i < TestMap["Ws"].size(); ++i)
  {
    cout<< "Ws:" << TestMap["Ws"][i] <<"  IA:"<<TestMap["Sigmas_IA"][i]<<"  IA_Err:"<<TestMap["Sigmas_IA_Err"][i]<<endl;
  }

}