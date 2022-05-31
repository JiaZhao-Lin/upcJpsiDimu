#include "../common/function.C"
#include "PhysParameters.h"


const TVector2 CenterPoint{0.02,	0.05};
const double SearchX_Range 	= 0.02;
const double SearchY_Range 	= 0.05;
const double StepSizeX		= 0.0001;
const double StepSizeY		= 0.0001;
const double mTinyNum = 1e-5;

TH2D * h_Chi2Contour0 = new TH2D("h_Chi2Contour0","Chi^2;#sigma(y1);#sigma(y2)",
												2.*SearchX_Range/StepSizeX+1, CenterPoint.X() - SearchX_Range - mTinyNum, CenterPoint.X() + SearchX_Range + mTinyNum,
												2.*SearchY_Range/StepSizeY+1, CenterPoint.Y() - SearchY_Range - mTinyNum, CenterPoint.Y() + SearchY_Range + mTinyNum);
TH2D * h_Chi2Contour1 = new TH2D("h_Chi2Contour1","Chi^2;#sigma(y1);#sigma(y2)",
												2.*SearchX_Range/StepSizeX+1, CenterPoint.X() - SearchX_Range - mTinyNum, CenterPoint.X() + SearchX_Range + mTinyNum,
												2.*SearchY_Range/StepSizeY+1, CenterPoint.Y() - SearchY_Range - mTinyNum, CenterPoint.Y() + SearchY_Range + mTinyNum);
TH2D * h_Chi2Contour2 = new TH2D("h_Chi2Contour2","Chi^2;#sigma(y1);#sigma(y2)",
												2.*SearchX_Range/StepSizeX+1, CenterPoint.X() - SearchX_Range - mTinyNum, CenterPoint.X() + SearchX_Range + mTinyNum,
												2.*SearchY_Range/StepSizeY+1, CenterPoint.Y() - SearchY_Range - mTinyNum, CenterPoint.Y() + SearchY_Range + mTinyNum);


std::vector<TH2D *> h_Chi2Contours = {h_Chi2Contour0, h_Chi2Contour1, h_Chi2Contour2};

double Cal_Chi2(const TVector2 point, const int RapBin)
{
	double chi2 = 0;

	chi2 += pow( Xsec_0n0n[RapBin] - point * flux_0n0n[RapBin], 2) / pow( XsecErr_0n0n[RapBin], 2);
	chi2 += pow( Xsec_0nXnSum[RapBin] - point * flux_0nXnSum[RapBin], 2) / pow( XsecErr_0nXnSum[RapBin], 2);
	chi2 += pow( Xsec_XnXn[RapBin] - point * flux_XnXn[RapBin], 2) / pow( XsecErr_XnXn[RapBin], 2);

	chi2 /= 6;
	return chi2;
}

TVector2 SearchPoint(const int RapBin)
{
	const int NStepX			= 2*SearchX_Range/StepSizeX;
	const int NStepY			= 2*SearchY_Range/StepSizeY;

	double LeastChi2 = Cal_Chi2(CenterPoint,	RapBin);
	TVector2 LeastChi2Point	(CenterPoint);
	TVector2 StartPoint		(CenterPoint.X() - SearchX_Range, CenterPoint.Y() - SearchY_Range);

	for (int i = 0; i < NStepX + 1; ++i)
	{
		for (int j = 0; j < NStepY + 1; ++j)
		{
			TVector2 CurrentPoint(StartPoint.X() + StepSizeX * i,	StartPoint.Y() + StepSizeY * j);
			double CurrentChi2 = Cal_Chi2(CurrentPoint, RapBin);

			h_Chi2Contours[RapBin]->Fill(CurrentPoint.X(), CurrentPoint.Y(), CurrentChi2);
			
			if ( CurrentChi2 < LeastChi2)
			{
				LeastChi2 = CurrentChi2;
				LeastChi2Point = TVector2(CurrentPoint);
			}
		}
	}
	cout << Form("The LeastChi2 Point (\u03C3(y = %.1f), \u03C3(y = %.1f)) is (%f, %f) with chi2: %.2f", Rap[RapBin], -Rap[RapBin],
																				LeastChi2Point.X(), LeastChi2Point.Y(),LeastChi2)<<endl;

	return LeastChi2Point;
}

void findLeastChi2()
{
	std::vector<TVector2> LeastChi2s;

	for (int i = 0; i < Rap.size(); ++i)
	{
		LeastChi2s.push_back(SearchPoint(i));
	}
	// Double_t level0[6] = {0,Cal_Chi2(LeastChi2s[0],0) + 1 ,Cal_Chi2(LeastChi2s[0],0) + 2,Cal_Chi2(LeastChi2s[0],0) + 3,Cal_Chi2(LeastChi2s[0],0) + 4,Cal_Chi2(LeastChi2s[0],0) + 5};
	// TCanvas * c0 = new TCanvas();
	// h_Chi2Contours[0]->Draw("cont1");
	// h_Chi2Contours[0]->SetContour(6,level0);
	// drawLatex(0.15, 0.85, Form("(#sigma(y = %.1f), #sigma(y = %.1f)) = (%.4f, %.4f), #chi^{2}: %.2f", Rap[0], -Rap[0],
	// 									LeastChi2s[0].X(), LeastChi2s[0].Y(),Cal_Chi2(LeastChi2s[0],0)),      42,       0.05,      1);
	// c0->SaveAs("outplots/Chi2Contours_0.png");

	// Double_t level1[6] = {0,Cal_Chi2(LeastChi2s[1],1) + 1 ,Cal_Chi2(LeastChi2s[1],1) + 2,Cal_Chi2(LeastChi2s[1],1) + 3,Cal_Chi2(LeastChi2s[1],1) + 4,Cal_Chi2(LeastChi2s[1],1) + 5};
	// TCanvas * c1 = new TCanvas();
	// h_Chi2Contours[1]->Draw("cont1");
	// h_Chi2Contours[1]->SetContour(6,level1);
	// drawLatex(0.15, 0.85, Form("(#sigma(y = %.1f), #sigma(y = %.1f)) = (%.4f, %.4f), #chi^{2}: %.2f", Rap[1], -Rap[1],
	// 									LeastChi2s[1].X(), LeastChi2s[1].Y(),Cal_Chi2(LeastChi2s[1],1)),      42,       0.05,      1);
	// c1->SaveAs("outplots/Chi2Contours_1.png");

	for (int i = 0; i < Rap.size(); ++i)
	{
		Double_t level0[6] = {0,Cal_Chi2(LeastChi2s[i],i) + 1 ,Cal_Chi2(LeastChi2s[i],i) + 2,Cal_Chi2(LeastChi2s[i],i) + 3,Cal_Chi2(LeastChi2s[i],i) + 4,Cal_Chi2(LeastChi2s[i],i) + 5};
		TCanvas * c = new TCanvas();
		h_Chi2Contours[i]->Draw("cont1");
		h_Chi2Contours[i]->SetContour(6,level0);
		drawLatex(0.15, 0.85, Form("(#sigma(y = %.2f), #sigma(y = %.2f)) = (%.4f, %.4f), #chi^{2}: %.2f", Rap[i], -Rap[i],
											LeastChi2s[i].X(), LeastChi2s[i].Y(),Cal_Chi2(LeastChi2s[i],i)),      42,       0.045,      1);
		c->SaveAs(Form("outplots/Chi2Contours_%d.png",i));

	}
}