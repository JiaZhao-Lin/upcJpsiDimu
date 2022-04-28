
TMatrixD MigrationMatrix(4,4);
TMatrixD MigrationMatrixInvert(4,4);

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

	MigrationMatrix = TMatrixD(4,	4,	MigrationMatrixElement);
	MigrationMatrixInvert = MigrationMatrix; MigrationMatrixInvert.Invert();

	cout<<"getPileUp_CorrFactor--->The MigrationMatrix is:"<<endl<<endl;
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
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
	for (int i = 0; i < nDiffRapBins; ++i)
	{
		double tem[4] = {NJpsi_inMFit[1][i],	NJpsi_inMFit[2][i],	NJpsi_inMFit[3][i],	NJpsi_inMFit[5][i]};
		TMatrixD N_Meastured_iy(4,	1, tem); 
		TMatrixD N_True_iy = MigrationMatrixInvert * N_Meastured_iy;
		NJpsi_inMFit[1][i] = N_True_iy[0][0];
		NJpsi_inMFit[2][i] = N_True_iy[1][0];
		NJpsi_inMFit[3][i] = N_True_iy[2][0];
		NJpsi_inMFit[5][i] = N_True_iy[3][0];

		NJpsi_inMFit[4][i] = N_True_iy[1][0] + N_True_iy[2][0];

		auto delta0 = sqrt( pow(MigrationMatrixInvert[0][0] * NerrJpsi_inMFit[1][i],2) + 
							pow(MigrationMatrixInvert[0][1] * NerrJpsi_inMFit[2][i],2) + 
							pow(MigrationMatrixInvert[0][2] * NerrJpsi_inMFit[3][i],2) + 
							pow(MigrationMatrixInvert[0][3] * NerrJpsi_inMFit[5][i],2) );
		auto delta1 = sqrt( pow(MigrationMatrixInvert[1][0] * NerrJpsi_inMFit[1][i],2) + 
							pow(MigrationMatrixInvert[1][1] * NerrJpsi_inMFit[2][i],2) + 
							pow(MigrationMatrixInvert[1][2] * NerrJpsi_inMFit[3][i],2) + 
							pow(MigrationMatrixInvert[1][3] * NerrJpsi_inMFit[5][i],2) );
		auto delta2 = sqrt( pow(MigrationMatrixInvert[2][0] * NerrJpsi_inMFit[1][i],2) + 
							pow(MigrationMatrixInvert[2][1] * NerrJpsi_inMFit[2][i],2) + 
							pow(MigrationMatrixInvert[2][2] * NerrJpsi_inMFit[3][i],2) + 
							pow(MigrationMatrixInvert[2][3] * NerrJpsi_inMFit[5][i],2) );
		auto delta3 = sqrt( pow(MigrationMatrixInvert[3][0] * NerrJpsi_inMFit[1][i],2) + 
							pow(MigrationMatrixInvert[3][1] * NerrJpsi_inMFit[2][i],2) + 
							pow(MigrationMatrixInvert[3][2] * NerrJpsi_inMFit[3][i],2) + 
							pow(MigrationMatrixInvert[3][3] * NerrJpsi_inMFit[5][i],2) );
		NerrJpsi_inMFit[1][i] = delta0;
		NerrJpsi_inMFit[2][i] = delta1;
		NerrJpsi_inMFit[3][i] = delta2;
		NerrJpsi_inMFit[5][i] = delta3;
		
		NerrJpsi_inMFit[4][i] = sqrt( pow(NerrJpsi_inMFit[2][i],2) + pow(NerrJpsi_inMFit[3][i],2) );

	}

}
