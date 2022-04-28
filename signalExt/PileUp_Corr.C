
const int nDmatrixs = 4; 
TMatrixD MigrationMatrix(nDmatrixs,nDmatrixs);
TMatrixD MigrationMatrixInvert(nDmatrixs,nDmatrixs);

TMatrixD getPileUp_CorrFactor()
{
	const double Pzb00 = 0.889779;
	const double Pzb0X = 0.0275305+0.0255331;
	const double PzbX0 = 0.0256458+0.0252;
	const double PzbXX = 0.00150574+0.0019427+0.0014546+0.00140865;

	const double P00_00 = Pzb00;
	const double P0X_00 = Pzb0X; const double P0X_0X = Pzb00 + Pzb0X;
	const double PX0_00 = PzbX0; const double PX0_X0 = Pzb00 + PzbX0;
	const double PXX_00 = PzbXX; const double PXX_0X = PzbX0 + PzbXX; const double PXX_X0 = Pzb0X + PzbXX;
	const double PXX_XX = 1.0;

	const double MigrationMatrixElement[16] = {	P00_00,	0,		0,		0,
											   	P0X_00,	P0X_0X,	0,		0,
											   	PX0_00,	0,		PX0_X0,	0,
											   	PXX_00,	PXX_0X,	PXX_X0,	PXX_XX };

	MigrationMatrix       = TMatrixD(nDmatrixs,	nDmatrixs,	MigrationMatrixElement);
	MigrationMatrixInvert = MigrationMatrix; MigrationMatrixInvert.Invert();

	cout<<"getPileUp_CorrFactor--->The MigrationMatrix is:"<<endl<<endl;
	for (int i = 0; i < nDmatrixs; ++i)
	{
		for (int j = 0; j < nDmatrixs; ++j)
		{
			cout<<MigrationMatrixInvert[i][j]<<",	";
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
			NJpsi_inMFit[2][iy],
			NJpsi_inMFit[3][iy],
			NJpsi_inMFit[5][iy]
		};

		TMatrixD N_Meastured_iy(nDmatrixs,	1, tem); 
		TMatrixD N_True_iy = MigrationMatrixInvert * N_Meastured_iy;

		cout<<"N_True_iy.GetNrows(): "<<N_True_iy.GetNrows()<<endl;
		cout<<"N_True_iy.GetNcols(): "<<N_True_iy.GetNcols()<<endl;
		
		NJpsi_inMFit[1][iy] = N_True_iy[0][0];
		NJpsi_inMFit[2][iy] = N_True_iy[1][0];
		NJpsi_inMFit[3][iy] = N_True_iy[2][0];
		NJpsi_inMFit[5][iy] = N_True_iy[3][0];

		NJpsi_inMFit[4][iy] = N_True_iy[1][0] + N_True_iy[2][0]; //0nXnSum =  0nXn + Xn0n

		auto delta0 = sqrt( pow(MigrationMatrixInvert[0][0] * NerrJpsi_inMFit[1][iy], 2) +
							pow(MigrationMatrixInvert[0][1] * NerrJpsi_inMFit[2][iy], 2) +
							pow(MigrationMatrixInvert[0][2] * NerrJpsi_inMFit[3][iy], 2) +
							pow(MigrationMatrixInvert[0][3] * NerrJpsi_inMFit[5][iy], 2) );

		auto delta1 = sqrt( pow(MigrationMatrixInvert[1][0] * NerrJpsi_inMFit[1][iy], 2) +
							pow(MigrationMatrixInvert[1][1] * NerrJpsi_inMFit[2][iy], 2) +
							pow(MigrationMatrixInvert[1][2] * NerrJpsi_inMFit[3][iy], 2) +
							pow(MigrationMatrixInvert[1][3] * NerrJpsi_inMFit[5][iy], 2) );
		
		auto delta2 = sqrt( pow(MigrationMatrixInvert[2][0] * NerrJpsi_inMFit[1][iy], 2) +
							pow(MigrationMatrixInvert[2][1] * NerrJpsi_inMFit[2][iy], 2) +
							pow(MigrationMatrixInvert[2][2] * NerrJpsi_inMFit[3][iy], 2) +
							pow(MigrationMatrixInvert[2][3] * NerrJpsi_inMFit[5][iy], 2) );
		
		auto delta3 = sqrt( pow(MigrationMatrixInvert[3][0] * NerrJpsi_inMFit[1][iy], 2) +
							pow(MigrationMatrixInvert[3][1] * NerrJpsi_inMFit[2][iy], 2) +
							pow(MigrationMatrixInvert[3][2] * NerrJpsi_inMFit[3][iy], 2) +
							pow(MigrationMatrixInvert[3][3] * NerrJpsi_inMFit[5][iy], 2) );
		
		NerrJpsi_inMFit[1][iy] = delta0;
		NerrJpsi_inMFit[2][iy] = delta1;
		NerrJpsi_inMFit[3][iy] = delta2;
		NerrJpsi_inMFit[5][iy] = delta3;
		
		NerrJpsi_inMFit[4][iy] = sqrt( pow(NerrJpsi_inMFit[2][iy],2) + pow(NerrJpsi_inMFit[3][iy],2) ); //0nXnSumErr = 
	}
}
