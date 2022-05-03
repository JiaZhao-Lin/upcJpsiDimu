
const int nDmatrixs = 3; 
TMatrixD MigrationMatrix(nDmatrixs,nDmatrixs);
TMatrixD MigrationMatrixInvert(nDmatrixs,nDmatrixs);

const double f_00	=	0.889779, f_0X	=	(0.053063600+0.050845800)/2,	f_XX	=	0.0014546+0.00140865+0.00150574+0.0019427;
std::vector<double> f_v{f_00, f_0X,	f_XX};

double Cal_Chi2(const double Ps, const double Pm)
{
	std::vector<double> f_Err_v;
	TVector3	f_ZB	=	{f_00,	f_0X,	f_XX};
	double nDen = 1e6; double nDenErr = sqrt(nDen);
	for (int i = 0; i < 3; ++i)
	{
		double ratioErr = TMath::Abs( ((1.-2.*f_v[i])*pow(sqrt(nDen*f_v[i]),2)+pow(f_v[i]*nDenErr,2))/pow(nDen,2) );
		ratioErr = sqrt(ratioErr);
		f_Err_v.push_back(ratioErr);
	}

	double A = pow( (1 - Ps), 2) * (1-Pm);
	double B = Ps * (1-Pm) * (1-Ps);
	double C = Pm + (1-Pm) * Ps*Ps;

	TVector3	tem(A,	B,	C);
	double chi2 = 0;
	chi2 += pow((f_ZB - tem).X(),2) / pow(f_Err_v[0],2);
	chi2 += pow((f_ZB - tem).Y(),2) / pow(f_Err_v[1],2);
	chi2 += pow((f_ZB - tem).Z(),2) / pow(f_Err_v[2],2);
	return chi2;
}

double Cal_Squares(const double Ps, const double Pm)
{
	TVector3	f_ZB	=	{f_00,	f_0X,	f_XX};

	double A = pow( (1 - Ps), 2) * (1-Pm);
	double B = Ps * (1-Pm) * (1-Ps);
	double C = Pm + (1-Pm) * Ps*Ps;

	TVector3	tem(A,	B,	C);

	return (f_ZB - tem).Mag2();
}

void getInteractionProbability()
{
	const TVector2 CenterPoint{0.0559,	0.00174};
	const double SearchPs_Range 	= 0.01;
	const double SearchPm_Range 	= 0.01;
	const double StepSizePs		= 0.0001;
	const double StepSizePm		= 0.00001;
	const int NStepPs			= 2*SearchPs_Range/StepSizePs;
	const int NStepPm			= 2*SearchPm_Range/StepSizePm;

	double LeastSquares = Cal_Squares(CenterPoint.X(),	CenterPoint.Y());
	TVector2 LeastSquaresPoint	(CenterPoint);
	TVector2 StartPoint		(CenterPoint.X() - SearchPs_Range, CenterPoint.Y() - SearchPm_Range);


	for (int i = 0; i < NStepPs + 1; ++i)
	{
		for (int j = 0; j < NStepPm + 1; ++j)
		{
			TVector2 CurrentPoint(StartPoint.X() + StepSizePs * i,	StartPoint.Y() + StepSizePm * j);
			double CurrentSquares = Cal_Squares(CurrentPoint.X(),	CurrentPoint.Y());

			if ( CurrentSquares < LeastSquares)
			{
				LeastSquares = CurrentSquares;
				LeastSquaresPoint = TVector2(CurrentPoint);
			}
		}
	}
	cout<<"getInteractionProbability-->"<<endl;
	cout<<"Ps: "<<LeastSquaresPoint.X()<<" , Pm: "<<LeastSquaresPoint.Y()<<endl<<endl;
	cout<<"us: "<<-log(1-LeastSquaresPoint.X())<<" , um: "<<-log(1-LeastSquaresPoint.Y())<<endl<<endl;
}


TMatrixD getPileUp_CorrFactor()
{
	// const double Ps = 0.0559; const double Pm = 0.00174;	//ATLAST Number
	const double Ps = 0.0552; const double Pm = 0.00325;	//From CMS ZB Data with Least Square Method

	const double MigrationMatrixElement[nDmatrixs*nDmatrixs] = {	(1-Ps)*(1-Ps)*(1-Pm)	,		0 				,	0,
											   						2*Ps*(1-Ps)*(1-Pm)		,		(1-Ps)*(1-Pm)	,	0,
											   						Pm+Ps*Ps*(1-Pm)			,		Pm+Ps*(1-Pm)	,	1};

	MigrationMatrix       = TMatrixD(nDmatrixs,	nDmatrixs,	MigrationMatrixElement);
	MigrationMatrixInvert = MigrationMatrix; MigrationMatrixInvert.Invert();

	cout<<"getPileUp_CorrFactor--->The MigrationMatrix is:"<<endl<<endl;
	for (int i = 0; i < nDmatrixs; ++i)
	{
		for (int j = 0; j < nDmatrixs; ++j)
		{
			cout<<MigrationMatrix[i][j]<<",	";
		}
		cout<<endl;
	}
	cout<<endl;

	return MigrationMatrixInvert;
}

void PileUp_Corr( double NJpsi_inMFit[][7], double NerrJpsi_inMFit[][7], const double nDiffRapBins )
{
	getPileUp_CorrFactor();

	for (int iy = 0; iy < nDiffRapBins; ++iy)
	{
		double tem[nDmatrixs] = 
		{
			NJpsi_inMFit[1][iy],
			NJpsi_inMFit[4][iy],
			NJpsi_inMFit[5][iy]
		};

		TMatrixD N_Meastured_iy(nDmatrixs,	1, tem); 
		TMatrixD N_True_iy = MigrationMatrixInvert * N_Meastured_iy;

		cout<<"N_True_iy.GetNrows(): "<<N_True_iy.GetNrows()<<endl;
		cout<<"N_True_iy.GetNcols(): "<<N_True_iy.GetNcols()<<endl;
		
		NJpsi_inMFit[1][iy] = N_True_iy[0][0];
		NJpsi_inMFit[4][iy] = N_True_iy[1][0];
		NJpsi_inMFit[5][iy] = N_True_iy[2][0];

		auto delta0 = sqrt( pow(MigrationMatrixInvert[0][0] * NerrJpsi_inMFit[1][iy], 2) +
							pow(MigrationMatrixInvert[0][1] * NerrJpsi_inMFit[4][iy], 2) +
							pow(MigrationMatrixInvert[0][2] * NerrJpsi_inMFit[5][iy], 2) );

		auto delta1 = sqrt( pow(MigrationMatrixInvert[1][0] * NerrJpsi_inMFit[1][iy], 2) +
							pow(MigrationMatrixInvert[1][1] * NerrJpsi_inMFit[4][iy], 2) +
							pow(MigrationMatrixInvert[1][2] * NerrJpsi_inMFit[5][iy], 2) );
		
		auto delta2 = sqrt( pow(MigrationMatrixInvert[2][0] * NerrJpsi_inMFit[1][iy], 2) +
							pow(MigrationMatrixInvert[2][1] * NerrJpsi_inMFit[4][iy], 2) +
							pow(MigrationMatrixInvert[2][2] * NerrJpsi_inMFit[5][iy], 2) );
		
		NerrJpsi_inMFit[1][iy] = delta0;
		NerrJpsi_inMFit[4][iy] = delta1;
		NerrJpsi_inMFit[5][iy] = delta2;
		
	}
}

void PileUp_Corr()
{
	getInteractionProbability();
}
