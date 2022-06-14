#ifndef PHYSPARAETERS_H
#define PHYSPARAETERS_H

//---------------------------------------Systematic Uncertainties---------------------------------------------
const double Lumi_Uncer = 0.015;
const double Jpsi2MuMuBR_Uncer = 0.033/5.961;

//------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------
//-------------------------------------------Run2 ALICE MidRap------------------------------------------------
const std::vector<double> ALICE_Run2_MidRap_W  				= 	{ 124.67672		},
						ALICE_Run2_MidRap_Sigma	    		= 	{ 4.07/2/86.483 },
						ALICE_Run2_MidRap_Sigma_StatErr		=	{ ALICE_Run2_MidRap_Sigma[0] * 0.11/4.07 },
						ALICE_Run2_MidRap_Sigma_SysErrLow	=	{ ALICE_Run2_MidRap_Sigma[0] * TMath::Hypot( 0.23/4.07,	 0.05)	}, // 5% error for flux
						ALICE_Run2_MidRap_Sigma_SysErrHig	=	{ ALICE_Run2_MidRap_Sigma[0] * TMath::Hypot( 0.23/4.07,	 0.05)	}, // 5% error for flux
						ALICE_Run2_MidRap_x	    			= 	{ 0.00061691554 },
						ALICE_Run2_MidRap_Sigma_IA	    	= 	{ 0.0658072 	},
						ALICE_Run2_MidRap_Sigma_IA_Err		= 	{ 0.00355359 	},
						ALICE_Run2_MidRap_R	    			= 	{ sqrt(ALICE_Run2_MidRap_Sigma[0] / ALICE_Run2_MidRap_Sigma_IA[0]) },
						ALICE_Run2_MidRap_R_StatErr	    	= 	{ 0.5 * ALICE_Run2_MidRap_R[0] * ALICE_Run2_MidRap_Sigma_StatErr[0]/ALICE_Run2_MidRap_Sigma[0] },
						ALICE_Run2_MidRap_R_SysErrLow  		= 	{ 0.5 * ALICE_Run2_MidRap_R[0] * TMath::Hypot( ALICE_Run2_MidRap_Sigma_SysErrLow[0]/ALICE_Run2_MidRap_Sigma[0],
																								ALICE_Run2_MidRap_Sigma_IA_Err[0]/ALICE_Run2_MidRap_Sigma_IA[0]) },
						ALICE_Run2_MidRap_R_SysErrHig  		= 	{ 0.5 * ALICE_Run2_MidRap_R[0] * TMath::Hypot( ALICE_Run2_MidRap_Sigma_SysErrHig[0]/ALICE_Run2_MidRap_Sigma[0],
																								ALICE_Run2_MidRap_Sigma_IA_Err[0]/ALICE_Run2_MidRap_Sigma_IA[0]) };
//-------------------------------------------Run2 ALICE FwdRap------------------------------------------------
const std::vector<double> ALICE_Run2_FwdRap_y  				= 	{ -3.875,		-3.625 	},
						ALICE_Run2_FwdRap_dSigmady			=	{ 1.615,		1.938 	},
						ALICE_Run2_FwdRap_dSigmady_StatErr	=	{ 0.060,		0.0420 	},
						ALICE_Run2_FwdRap_dSigmady_SysErrLow=	{ 0.147,		0.190 	},
						ALICE_Run2_FwdRap_dSigmady_SysErrHig=	{ 0.135,		0.166 	},
						ALICE_Run2_FwdRap_Flux 				=	{ 206.795,		199.037 },
						ALICE_Run2_FwdRap_Flux_Err 			=	{ 206.795*0.05,	199.037*0.05}, // 5% error for flux
						ALICE_Run2_FwdRap_Sigma_IA 			=	{ 0.0104829,	0.0128178		},
						ALICE_Run2_FwdRap_Sigma_IA_Err 		=	{ 0.000566075,	0.000692159	};

std::vector<double> 	ALICE_Run2_FwdRap_W  				= 	{},	
						ALICE_Run2_FwdRap_Sigma				=	{},	
						ALICE_Run2_FwdRap_Sigma_StatErr		=	{},
						ALICE_Run2_FwdRap_Sigma_SysErrLow	=	{},
						ALICE_Run2_FwdRap_Sigma_SysErrHig	=	{},
						ALICE_Run2_FwdRap_x 				=	{},
						ALICE_Run2_FwdRap_R 				=	{},
						ALICE_Run2_FwdRap_R_StatErr 		=	{},
						ALICE_Run2_FwdRap_R_SysErrLow 		=	{},
						ALICE_Run2_FwdRap_R_SysErrHig 		=	{};
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------
//-------------------------------------------2022 CMS 2 Rap Bin-----------------------------------------------
// const std::vector<double>  Rap 				= { -2.2,  -1.8 };
// const std::vector<double>  Rap_Err  		= { 0.20,  0.20 };
// //----------------------------------------Rap:  -2.2,		2.2 		  -1.8,		1.8 
// const std::vector<TVector2>  flux_0n0n    = { { 124.474,	7.355 },	{ 112.075,	13.158	}	}; //N(y)
// const std::vector<TVector2>  flux_0nXnSum = { { 18.291,		7.517 },	{ 18.262,	10.386	}	}; //N(y)
// const std::vector<TVector2>  flux_XnXn    = { { 6.521,		3.856 },	{ 6.520,	4.812	}	}; //N(y)
// // const std::vector<TVector2>  flux_0n0n    = { { 725.484, 0.526 }, { 437.863, 1.405 }   }; //n(E)
// // const std::vector<TVector2>  flux_0nXnSum = { { 106.609, 0.538 }, { 71.349, 1.109 } }; //n(E)
// // const std::vector<TVector2>  flux_XnXn    = { { 38.009, 0.276 },  { 25.471, 0.514 } }; //n(E)

// //-------------------------------------Rap: |-2.2|,		|-1.8| 
// const std::vector<double>  Xsec_0n0n      = { 2.65,		3.099 },     XsecErr_0n0n      = { 0.029,	0.08  };
// const std::vector<double>  Xsec_0nXnSum   = { 0.979,	1.203 },     XsecErr_0nXnSum   = { 0.016,	0.023 };
// const std::vector<double>  Xsec_XnXn      = { 0.259,	0.352 },     XsecErr_XnXn      = { 0.008,	0.016 };
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------
//-------------------------------------------2022 CMS 3 Rap Bin-----------------------------------------------
// const std::vector<double>  Rap 				= { -2.25,  -2.0, -1.75};
// const std::vector<double>  Rap_Err  		= { 0.15,   0.10, 0.15 };
//----------------------------------------Rap:  -2.25,		2.25 		  -2.0,		2.0 		  -1.75,	1.75
// const std::vector<TVector2>  flux_0n0n    = { { 126.025,	6.806 },	{ 118.272,	10.012	},	{ 110.526,	14.027	}	}; //N(y)
// const std::vector<TVector2>  flux_0nXnSum = { { 18.294,		7.151 },	{ 18.279,	8.980	},	{ 18.257,	10.721	}	}; //N(y)
// const std::vector<TVector2>  flux_XnXn    = { { 6.521,		3.716 },	{ 6.521,	4.367	},	{ 6.519,	4.911	}	}; //N(y)
//----------------------------------------Rap:  1.75,	-1.75 		  2.0,		-2.0 		  2.25,		-2.25
// const std::vector<TVector2>  flux_0n0n    = { { 17.281,	116.004 },	{ 12.755,	123.755	},	{ 9.011,	131.512	}	}; //N(y)
// const std::vector<TVector2>  flux_0nXnSum = { { 11.857,	18.274 },	{ 10.231,	18.290	},	{ 8.448,	18.301	}	}; //N(y)
// const std::vector<TVector2>  flux_XnXn    = { { 5.234,	6.520 },	{ 4.766,	6.521	},	{ 4.184,	6.522	}	}; //N(y)


//----------------CrystalBall:			  Rap: |-2.25|,	|-2.0|,	|-1.75| --------------------------------------
// const std::vector<double>  Xsec_0n0n      = { 2.512,	2.928,	3.008 },     XsecErr_0n0n      = { 0.035,	0.032,	0.057 };
// const std::vector<double>  Xsec_0nXnSum   = { 0.929,	1.092,	1.267 },     XsecErr_0nXnSum   = { 0.018,	0.018,	0.035 };
// const std::vector<double>  Xsec_XnXn      = { 0.259,	0.310,	0.324 },     XsecErr_XnXn      = { 0.009,	0.009,	0.017 };
//------------------------------------------------------------------------------------------------------------

//----------------CrystalBall PU ATLAST:  Rap: |-2.25|,	|-2.0|,	|-1.75| --------------------------------------
// const std::vector<double>  Xsec_0n0n      = { 2.824,	3.307,	3.371 },     XsecErr_0n0n      = { 0.039,	0.042,	0.064 };
// const std::vector<double>  Xsec_0nXnSum   = { 0.700,	0.824,	0.998 },     XsecErr_0nXnSum   = { 0.019,	0.020,	0.037 };
// const std::vector<double>  Xsec_XnXn      = { 0.201,	0.243,	0.247 },     XsecErr_XnXn      = { 0.009,	0.014,	0.017 };
//------------------------------------------------------------------------------------------------------------

////----------------CrystalBall PU CMS:	  Rap: |-2.25|,	|-2.0|,	|-1.75| --------------------------------------
//const std::vector<double>  Xsec_0n0n      = { 2.824,	3.307,	3.371 },     XsecErr_0n0n      = { 0.039,	0.042,	0.064 };
//const std::vector<double>  Xsec_0nXnSum   = { 0.700,	0.824,	0.998 },     XsecErr_0nXnSum   = { 0.019,	0.020,	0.037 };
//const std::vector<double>  Xsec_XnXn      = { 0.201,	0.243,	0.247 },     XsecErr_XnXn      = { 0.009,	0.014,	0.017 };
////------------------------------------------------------------------------------------------------------------

//----------------CrystalBall PU Shuai:	  Rap: |-2.25|,	|-2.0|,	|-1.75| --------------------------------------
 // const std::vector<double>  Xsec_0n0n      = { 2.823,	3.307,	3.371 },     XsecErr_0n0n      = { 0.039,	0.042,	0.064 };
 // const std::vector<double>  Xsec_0nXnSum   = { 0.702,	0.819,	1.006 },     XsecErr_0nXnSum   = { 0.019,	0.028,	0.037 };
 // const std::vector<double>  Xsec_XnXn      = { 0.201,	0.243,	0.246 },     XsecErr_XnXn      = { 0.009,	0.014,	0.017 };
//------------------------------------------------------------------------------------------------------------

//----------------CrystalBall + Gaussian: Rap: |-2.25|,	|-2.0|,	|-1.75| --------------------------------------
// const std::vector<double>  Xsec_0n0n      = { 2.549,	2.957,	3.019 },     XsecErr_0n0n      = { 0.040,	0.039,	0.060 };
// const std::vector<double>  Xsec_0nXnSum   = { 0.929,	1.096,	1.264 },     XsecErr_0nXnSum   = { 0.018,	0.018,	0.035 };
// const std::vector<double>  Xsec_XnXn      = { 0.258,	0.309,	0.327 },     XsecErr_XnXn      = { 0.017,	0.009,	0.018 };

//----------------DCB:					  Rap: |-2.25|,	|-2.0|,	|-1.75| --------------------------------------
// const std::vector<double>  Xsec_0n0n      = { 2.538,	2.947,	3.512 },     XsecErr_0n0n      = { 0.029,	0.030,	0.062 };
// const std::vector<double>  Xsec_0nXnSum   = { 0.929,	1.120,	1.369 },     XsecErr_0nXnSum   = { 0.016,	0.018,	0.035 };
// const std::vector<double>  Xsec_XnXn      = { 0.259,	0.331,	0.368 },     XsecErr_XnXn      = { 0.009,	0.009,	0.018 };
//------------------------------------------------------------------------------------------------------------
//Zaochen DCB
// const std::vector<double>  Xsec_0n0n      = { 2.551,	2.935,	3.031 },     XsecErr_0n0n      = { 0.030,	0.030,	0.055 };
// const std::vector<double>  Xsec_0nXnSum   = { 0.936,	1.096,	1.367 },     XsecErr_0nXnSum   = { 0.017,	0.017,	0.035 };
// const std::vector<double>  Xsec_XnXn      = { 0.257,	0.323,	0.383 },     XsecErr_XnXn      = { 0.010,	0.009,	0.018 };
//------------------------------------------------------------------------------------------------------------

////----------------AsyDCB:				  Rap: |-2.25|,	|-2.0|,	|-1.75| --------------------------------------
//const std::vector<double>  Xsec_0n0n      = { 2.544,	2.940,	3.500 },     XsecErr_0n0n      = { 0.029,	0.030,	0.061 };
//const std::vector<double>  Xsec_0nXnSum   = { 0.931,	1.130,	1.262 },     XsecErr_0nXnSum   = { 0.016,	0.018,	0.033 };
//const std::vector<double>  Xsec_XnXn      = { 0.257,	0.331,	0.325 },     XsecErr_XnXn      = { 0.009,	0.009,	0.016 };
////------------------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
//-----------------------------CGC Prediction 5.02TeV gamma + Pb -> JPsi + Pb --------------------------------
const std::vector<double>  CGC_JpsiNoFluct_W	= {37, 47, 59,75, 95, 121, 153, 194, 246, 311, 394, 599  },
						   CGC_JpsiNoFluct_CohXsec 		= {0.01759987558438179, 0.020489530220598236, 0.023324292435666014, 0.02714433695284619, 0.03137306416640146, 0.036389263389293006, 0.041829124040095714, 0.04793443973914736, 0.05473231107873076, 0.06212140529983476, 0.07032820925862308, 0.0874426817784331},
						   CGC_JpsiNoFluct_InCohXsec	= {0.005314082049208781, 0.006081113144220922, 0.00684657722834769, 0.008057272329685066, 0.009032600408812817, 0.010363958837189501, 0.01173848797350525, 0.013306012198701258, 0.014398401192654236, 0.016650908354141573, 0.01846420961715733, 0.022516649102744427};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
#endif


