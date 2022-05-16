#include "plotShadowingRatio.C"


void Relative_Uncer( std::map<TString, std::vector<double>> &m1,	std::map<TString, std::vector<double>> m2, const TString param)
{
	std::vector<double> out;
	for (int i = 0; i < m1[param].size(); ++i)
	{
		out.push_back( 100 * abs(m1[param][i] - m2[param][i]) / m2[param][i] );
	}
	m1[param+"_Uncer"] = out;
}

void plot_Uncer(std::map<TString, std::vector<double>> &Map, const std::map<TString, std::vector<double>> &Map_Default,
				TString paramX, TString paramY, int style, int color, TLegend *legend, TString legendName)
{
	Relative_Uncer(Map,	Map_Default,	paramY);

	auto gr	= new TGraph(Map[paramX].size(),		&Map[paramX][0],		&Map[paramY+"_Uncer"][0]);
	gr->SetMarkerSize(1.1);
	gr->SetMarkerColor(color);
	gr->Draw("psame");

	legend->AddEntry(gr, legendName, "p");
}

void getFitUncertainty()
{
	std::map<TString, std::vector<double>> CB_Poly3_Map 		= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly3_PUShuai_6RapBins.root");
	std::map<TString, std::vector<double>> CB_Poly4_Map 		= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_Poly4_PUShuai_6RapBins.root");
	std::map<TString, std::vector<double>> CB_FixCBA_Poly3_Map 	= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_FixCBA_Poly3_PUShuai_6RapBins.root");
	std::map<TString, std::vector<double>> CB_FixCBN_Poly3_Map 	= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_FixCBN_Poly3_PUShuai_6RapBins.root");
	std::map<TString, std::vector<double>> CB_FixCBAN_Poly3_Map = getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CB_FixCBAN_Poly3_PUShuai_6RapBins.root");
	std::map<TString, std::vector<double>> CBG_Poly3_Map 		= getParamsMap("../signalExt/JpsiXsecValues/JpsiXsec_CBG_Poly3_PUShuai_6RapBins.root");

	const std::map<TString, std::vector<double>> JpsiXsec_Default		= CB_Poly3_Map;

	// for (int i = 0; i < CB_Poly3_Map["Sigmas"].size(); ++i)
	// {
	// 	cout<<abs(CB_Poly3_Map["Sigmas"][i]-CBG_Poly3_Map["Sigmas"][i])/CB_Poly3_Map["Sigmas"][i]<<endl;
	// }

	//-----------------------------------W Unicertainty----------------------------------
	auto c1 = new TCanvas();
	auto Sigma_FitUncer	= new TH2D("Sigma_FitUncertainty", ";W_{#gammaPb} (GeV);#sigma Uncer.(%);", 10,14,520, 10,0,40);
	Sigma_FitUncer->Draw();

	TLegend  *legendW =  new TLegend(0.15, 0.70, 0.30, 0.90);
	legendW->SetFillStyle(0);
	legendW->SetTextSize(0.040);
	plot_Uncer(CB_Poly4_Map,			JpsiXsec_Default,	"Ws",	"Sigmas",	1,	1,	legendW,	"CB_Poly4");
	plot_Uncer(CB_FixCBA_Poly3_Map,		JpsiXsec_Default,	"Ws",	"Sigmas",	1,	2,	legendW,	"CB_FixCBA_Poly3");
	plot_Uncer(CB_FixCBN_Poly3_Map,		JpsiXsec_Default,	"Ws",	"Sigmas",	1,	3,	legendW,	"CB_FixCBN_Poly3");
	plot_Uncer(CB_FixCBAN_Poly3_Map,	JpsiXsec_Default,	"Ws",	"Sigmas",	1,	4,	legendW,	"CB_FixCBAN_Poly3");
	plot_Uncer(CBG_Poly3_Map,			JpsiXsec_Default,	"Ws",	"Sigmas",	1,	6,	legendW,	"CBG_Poly3");
	legendW->Draw("same");
	//-----------------------------------------------------------------------------------

	//-----------------------------------R Unicertainty----------------------------------
	auto c2 = new TCanvas();
	c2->SetLogx();
	auto R_FitUncer	= new TH2D("R_FitUncertainty", ";x;R^{Pb}_{g} Uncer.(%);", 10,4.0e-5,1.1e-2, 10, 0, 30);
	R_FitUncer->Draw();

	TLegend  *legendR =  new TLegend(0.15, 0.70, 0.30, 0.90);
	legendR->SetFillStyle(0);
	legendR->SetTextSize(0.040);
	plot_Uncer(CB_Poly4_Map,			JpsiXsec_Default,	"Xs",	"R",	1,	1,	legendR,	"CB_Poly4");
	plot_Uncer(CB_FixCBA_Poly3_Map,		JpsiXsec_Default,	"Xs",	"R",	1,	2,	legendR,	"CB_FixCBA_Poly3");
	plot_Uncer(CB_FixCBN_Poly3_Map,		JpsiXsec_Default,	"Xs",	"R",	1,	3,	legendR,	"CB_FixCBN_Poly3");
	plot_Uncer(CB_FixCBAN_Poly3_Map,	JpsiXsec_Default,	"Xs",	"R",	1,	4,	legendR,	"CB_FixCBAN_Poly3");
	plot_Uncer(CBG_Poly3_Map,			JpsiXsec_Default,	"Xs",	"R",	1,	6,	legendR,	"CBG_Poly3");
	legendR->Draw("same");
	//-----------------------------------------------------------------------------------
}