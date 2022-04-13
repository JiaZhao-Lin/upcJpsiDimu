#include "./ImpulseApproximation/runUPC_AAModel.C"
#include "fit2D.C"

const double JpsiMass = 3.096916;
const double Sqrt_s = 5020;
const double Gamma_beam = 2672.9;
const double Mass_N = (0.93827+0.93957) / 2;

// const std::vector<double>	Raps			= { -2.2,	-1.8,	1.8,	2.2};
// const std::vector<double>	Raps_Err		= { 0.2,	0.2,	0.2,	0.2};
// const std::vector<double> 	Sigmas 		= {	0.0180,	0.0200,	0.0695,	0.0604},  
// 							Sigmas_Err 	= {	0.0003,	0.0009,	0.0027,	0.0018};

							// Sigmas_IA	= { 0.028276, 0.033986, 0.144246, 0.169381},
							// Sigmas_IA_Err={ 0.000566, 0.000680, 0.002885, 0.003388};

std::vector<double> Raps, Raps_Err;
std::vector<double> Sigmas, Sigmas_Err;
std::vector<double> Sigmas_IA, Sigmas_IA_Err;
std::vector<double> R, R_Err;

std::vector<double> 	Xs = {},  Xs_Err = {};
std::vector<double> 	Ws = {},  Ws_Err = {};

void Cal_R_Error()
{
	for (int i = 0; i < Raps.size(); ++i)
	{
		double tem = R[i]*R[i] * sqrt(pow(Sigmas_Err[i]/Sigmas[i],2) + pow(Sigmas_IA_Err[i]/Sigmas_IA[i],2) );
		R_Err.push_back( tem/R[i]/2);
	}
}

double y2x(const double y)
{
	return (JpsiMass / Sqrt_s) * exp(-y);
}

double y2W(const double y)
{
	return sqrt( (2 * Gamma_beam * Mass_N * JpsiMass) * exp(y) );
}

void plotResults()
{
	{auto Temp = fit2D(); Raps = Temp[0]; Raps_Err = Temp[1]; Sigmas = Temp[2]; Sigmas_Err = Temp[3];}
	for (int i = 0; i < Raps.size(); ++i)
	{
		Xs 		.push_back( y2x(Raps[i]) 	);
		Xs_Err	.push_back( 0 );
		Ws 		.push_back( y2W(Raps[i]) 	);
		Ws_Err 	.push_back( 0 );

		cout<<Form("y: %f,	x: %f,	W: %f", Raps[i], Xs[i], Ws[i])<<endl;
		cout<<"x from W: "<<pow(JpsiMass,2)/pow(Ws[i],2)<<endl;
	}
	{auto Temp = runUPC_AAModel(Ws); Sigmas_IA = Temp[0]; Sigmas_IA_Err = Temp[1];}

	for (int i = 0; i < Raps.size(); ++i){
		R.push_back( sqrt( Sigmas[i]/Sigmas_IA[i] ) );
	}
	Cal_R_Error();

	auto c = new TCanvas();
	c->SetLogx();

	TGraphErrors* gr = new TGraphErrors(Xs.size(), &Xs[0], &R[0], &Xs_Err[0], &R_Err[0]);
	gr->Draw("AP");
	gr->SetTitle(";x;R");
	gr->GetXaxis()->SetLimits(1e-5,1e-1);
	gr->GetYaxis()->SetRangeUser(0,1);

	c->SaveAs("outplots/RvsX.png");
}