//follow equation 12 of PLR 87.192301

#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"

#include <Riostream.h>
#include "TLegend.h"
#include "TLegendEntry.h"

#include "Math/IFunction.h"
#include <cmath>
#include "TSystem.h"
#include "TAxis.h"
#include "TPaveLabel.h"


const double R = 6.67; //6.67 +/- 0.03 fm
const double a = 0.56; // 0.56 +/-0.03 fm
const double pi= 3.1415926;
const double Q2= 0.0; //for coherent VM
const double Alpha = 1/137.;
const double c = 299792458.; //m/s
const double unit_fm2tob = 100.; //1b = 100 fm^2;
double Fv2ee(TString VMname="Jpsi");
double Mv(TString VMname="Jpsi");
double fun_term3(double *x, double *par);

void calBDLxsec(const TString VMname="Jpsi")
{
	const double fv    = Fv2ee(VMname);
	const double mv    = Mv(VMname);
	const double term1 = pow((2.*pi*pow(R,2)), 2)/(16.*pi);
	const double term2 = (3*fv*pow(mv,3))/(Alpha*pow(pow(mv,2)+Q2,2));

	//TF1* f4term3 = new TF1("f4term3", "(4./(x*pow(6.67,2)))*pow(TMath::BesselJ1(sqrt(x)*6.67),2)", 0., 10.);
	
	TF1* f4term3 = new TF1("f4term3", fun_term3, 0,1000, 0);
	
	const double term3 = f4term3->Integral(0., 1000);

	const double xsec = term1*term2*term3; //in unit of fm^2;

	const double xsec_mb = (xsec/unit_fm2tob)*1000.; //1b = 1000mb;

	cout<<"xsec from BDL: "<<xsec<<endl;

}


double Fv2ee(TString VMname="Jpsi")
{
	const double unit_kev2gev = 1.e-6;
	double Gamma_VM2ee = -99.;

	if(VMname=="Jpsi")
	{
		Gamma_VM2ee = 5.55; //Fee = 5.55 +/-0.14+/-0.02 kev PDG
	}
	else if(VMname=="Rho")
	{
		Gamma_VM2ee = 7.04; //Fee = 7.04 +/-0.06 kev PDG
	}
	
	return Gamma_VM2ee*unit_kev2gev;
}

double Mv(TString VMname="Jpsi")
{
	double mass = -99.;

	if(VMname=="Jpsi")
	{
		mass = 3.0969; //
	}
	else if(VMname=="Rho")
	{
		mass = 0.77526; //
	}
	
	return mass;
}
double fReject4QED(double *x, double *par)
{
	if((x[0]>2.80&& x[0]<3.35) || (x[0]>3.50 && x[0]<3.85))
	{
		TF1::RejectPoint();
		return 0;
	}

	return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3);
}

double fun_term3(double *x, double *par)
{
	double t    = x[0]; //positive t here

	//please put sqrt(t) and R in the same unit
	const double unit_fm2GeV = (1.e3/197.3); //1fm = 1./(197.3 MeV) = 1000./(197.3 GeV) 
	double temf = (4./(t*pow(R,2)))*pow(TMath::BesselJ1(sqrt(t)*R*unit_fm2GeV),2);

	return temf;
}
