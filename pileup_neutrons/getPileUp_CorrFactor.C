
void getPileUp_CorrFactor()
{
	const double Pzb00 = 0.889779;
	const double Pzb0X = 0.0275305+0.0255331;
	const double PzbX0 = 0.0256458+0.0252;
	const double PzbXX = 0.00150574+0.0019427+0.0014546+0.00140865;

	const double P00_00 = Pzb00;
	const double P0X_00 = Pzb0X; const double P0X_0X = Pzb00 + Pzb0X;
	const double PX0_00 = PzbX0; const double PX0_X0 = Pzb00 + PzbX0;
	const double PXX_00 = PzbXX; const double PXX_0X = Pzb0X + PzbXX; const double PXX_X0 = PzbX0 + PzbXX;
	const double PXX_XX = 1.0;
}
