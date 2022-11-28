#include "./ImpulseApproximation/ImpulseApprox.C"
#include "../common/DataReader.C"


void drawLTA(TString Case, TLegend * leg)
{
	std::vector<TString> Name = {"LTA SS", "LTA WS", "EPS09 central"};
	std::vector<int> colors = {2,8,4};

	MultiDataReader MDR({"../physicsFigures/inputfiles/LTA_Jpsi_strong_shadowing.dat", "../physicsFigures/inputfiles/LTA_Jpsi_weak_shadowing.dat", "../physicsFigures/inputfiles/EPS09_central_Jpsi.dat"},
						{"y", "AnAn", "0n0n", "0nXnSum", "XnXn"});

	for (int i = 0; i < MDR.SDRs.size()-1; ++i)
	{
		auto gr_AnAn 	= new TGraph(MDR.GetVec(i, "y").size(),	MDR.GetVec(i, "y").data(),	MDR.GetVec(i, "AnAn").data());
		auto gr_0n0n 	= new TGraph(MDR.GetVec(i, "y").size(),	MDR.GetVec(i, "y").data(),	MDR.GetVec(i, "0n0n").data());
		auto gr_0nXnSum = new TGraph(MDR.GetVec(i, "y").size(),	MDR.GetVec(i, "y").data(),	MDR.GetVec(i, "0nXnSum").data());
		auto gr_XnXn 	= new TGraph(MDR.GetVec(i, "y").size(),	MDR.GetVec(i, "y").data(),	MDR.GetVec(i, "XnXn").data());
	
		if (Case == "AnAn"){
			gr_AnAn->SetLineColor(colors[i]);
            gr_AnAn->SetLineStyle(9);
            gr_AnAn->SetLineWidth(3);
            gr_AnAn->Draw("same");
			leg->AddEntry(gr_AnAn,	Name[i].Data(),		"l" );
		}
		else if (Case == "0n0n"){
			gr_0n0n->SetLineColor(colors[i]);
            gr_0n0n->SetLineStyle(9);
            gr_0n0n->SetLineWidth(3);
            gr_0n0n->Draw("same");
			leg->AddEntry(gr_0n0n,	Name[i].Data(),		"l" );
		}
		else if (Case == "0nXnSum"){
			gr_0nXnSum->SetLineColor(colors[i]);
            gr_0nXnSum->SetLineStyle(9);
            gr_0nXnSum->SetLineWidth(3);
            gr_0nXnSum->Draw("same");
			leg->AddEntry(gr_0nXnSum,	Name[i].Data(),	"l" );
		}
		else if (Case == "XnXn"){
			gr_XnXn->SetLineColor(colors[i]);
            gr_XnXn->SetLineStyle(9);
            gr_XnXn->SetLineWidth(3);
            gr_XnXn->Draw("same");
			leg->AddEntry(gr_XnXn,	Name[i].Data(),		"l" );
		}
	}
}


void drawLTA_Sigmas_R(TString Case, TLegend * leg)
{
	std::vector<TString> Name = {"LTA_SS","LTA_WS"};
	std::vector<int> colors = {2,8,4};

	SingleDataReader SDR1("./inputfiles/LTA_Jpsi_Sigma_R_1.dat", 
						{"y", "Xs_Left", "R_StrSuppr_Left", "R_WeakSuppr_Left", "Ws_Left", "Sigmas_StrSuppr_Left", "Sigmas_WeakSuppr_Left", "Xs_Right", "R_StrSuppr_Right", "R_WeakSuppr_Right", "Ws_Right", "Sigmas_StrSuppr_Right", "Sigmas_WeakSuppr_Right"});
	SingleDataReader SDR2("./inputfiles/LTA_Jpsi_Sigma_R_2.dat",
						{"y", "Ws_Left", "Sigmas_StrSuppr_Left", "Sigmas_WeakSuppr_Left", "Ws_Right", "Sigmas_StrSuppr_Right", "Sigmas_WeakSuppr_Right"});

	std::map<TString, std::vector<double>> Map ={
		{"Xs",				{}},
		{"R_WeakSuppr",		{}},
		{"R_StrSuppr",		{}},
		{"Ws",				{}},
		{"Sigmas_WeakSuppr",{}},
		{"Sigmas_StrSuppr",	{}}
	};

	/*Using R with IA calculated from our code for now. There is 5-10% difference in calculated form factor value compares to Vadim*/
	Map.at("Xs")				=	SDR1.GetVec("Xs_Right");
	// Map.at("R_WeakSuppr")		=	SDR1.GetVec("R_WeakSuppr_Right");
	// Map.at("R_StrSuppr")		=	SDR1.GetVec("R_StrSuppr_Right");
	Map.at("Ws")				=	SDR2.GetVec("Ws_Right");
	Map.at("Sigmas_WeakSuppr")	=	SDR2.GetVec("Sigmas_WeakSuppr_Right");
	Map.at("Sigmas_StrSuppr")	=	SDR2.GetVec("Sigmas_StrSuppr_Right");
	auto Xs_temp 			= SDR1.GetVec("Xs_Left");			std::reverse(Xs_temp.begin(),Xs_temp.end());
	// auto R_WeakSuppr_temp 	= SDR1.GetVec("R_WeakSuppr_Left");	std::reverse(R_WeakSuppr_temp.begin(),R_WeakSuppr_temp.end());
	// auto R_StrSuppr_temp 	= SDR1.GetVec("R_StrSuppr_Left");	std::reverse(R_StrSuppr_temp.begin(),R_StrSuppr_temp.end());
	auto Ws_temp 			= SDR2.GetVec("Ws_Left");			std::reverse(Ws_temp.begin(),Ws_temp.end());
	auto Sigmas_WeakSuppr_temp 	= SDR2.GetVec("Sigmas_WeakSuppr_Left");	std::reverse(Sigmas_WeakSuppr_temp.begin(),Sigmas_WeakSuppr_temp.end());
	auto Sigmas_StrSuppr_temp 	= SDR2.GetVec("Sigmas_StrSuppr_Left");	std::reverse(Sigmas_StrSuppr_temp.begin(),Sigmas_StrSuppr_temp.end());

	Map.at("Xs")			.insert(Map.at("Xs").begin(),			Xs_temp.begin(),	Xs_temp.end());
	// Map.at("R_WeakSuppr")	.insert(Map.at("R_WeakSuppr").begin(),	R_WeakSuppr_temp.begin(),	R_WeakSuppr_temp.end());
	// Map.at("R_StrSuppr")	.insert(Map.at("R_StrSuppr").begin(),	R_StrSuppr_temp.begin(),	R_StrSuppr_temp.end());
	Map.at("Ws")			.insert(Map.at("Ws").begin(),			Ws_temp.begin(),	Ws_temp.end());
	Map.at("Sigmas_WeakSuppr")	.insert(Map.at("Sigmas_WeakSuppr").begin(),	Sigmas_WeakSuppr_temp.begin(),	Sigmas_WeakSuppr_temp.end());
	Map.at("Sigmas_StrSuppr")	.insert(Map.at("Sigmas_StrSuppr").begin(),	Sigmas_StrSuppr_temp.begin(),	Sigmas_StrSuppr_temp.end());
	
	getImpulseApprox(Map);
	for (int i = 0; i < Map.at("Ws").size(); ++i)
	{
		Map.at("R_WeakSuppr")		.push_back(sqrt(Map.at("Sigmas_WeakSuppr")[i]/Map.at("Sigmas_IA")[i]));
		Map.at("R_StrSuppr")		.push_back(sqrt(Map.at("Sigmas_StrSuppr")[i]/Map.at("Sigmas_IA")[i]));
	}

	auto gr_SigmasVsW_StrSuppr 	= new TGraph(Map.at("Ws").size(),	
												Map.at("Ws").data(),			Map.at("Sigmas_StrSuppr").data());
	auto gr_SigmasVsW_WeakSuppr = new TGraph(Map.at("Ws").size(),	
												Map.at("Ws").data(),			Map.at("Sigmas_WeakSuppr").data());
	auto gr_RvsX_StrSuppr	 	= new TGraph(Map.at("Xs").size(),	
												Map.at("Xs").data(),			Map.at("R_StrSuppr").data());
	auto gr_RvsX_WeakSuppr	 	= new TGraph(Map.at("Xs").size(),	
												Map.at("Xs").data(),			Map.at("R_WeakSuppr").data());

	if (Case == "Sigmas")
	{
		// gr_SigmasVsW_StrSuppr->SetFillColor(colors[0]);
   		// gr_SigmasVsW_StrSuppr->SetFillStyle(3005);
   		gr_SigmasVsW_StrSuppr->SetLineColor(colors[0]);
   		gr_SigmasVsW_StrSuppr->SetLineStyle(9);
   		gr_SigmasVsW_StrSuppr->SetLineWidth(3);
   		gr_SigmasVsW_StrSuppr->Draw("lsame");
		leg->AddEntry(gr_SigmasVsW_StrSuppr,	Name[0].Data(),		"l" );

		// gr_SigmasVsW_WeakSuppr->SetFillColor(colors[1]);
   		// gr_SigmasVsW_WeakSuppr->SetFillStyle(3005);
   		gr_SigmasVsW_WeakSuppr->SetLineColor(colors[1]);
   		gr_SigmasVsW_WeakSuppr->SetLineStyle(9);
   		gr_SigmasVsW_WeakSuppr->SetLineWidth(3);
   		gr_SigmasVsW_WeakSuppr->Draw("lsame");
		leg->AddEntry(gr_SigmasVsW_WeakSuppr,	Name[1].Data(),		"l" );
	}
	else if (Case == "R")
	{
		// gr_RvsX_StrSuppr->SetFillColor(colors[0]);
   		// gr_RvsX_StrSuppr->SetFillStyle(3005);
		gr_RvsX_StrSuppr->SetLineColor(colors[0]);
		gr_RvsX_StrSuppr->SetLineStyle(9);
		gr_RvsX_StrSuppr->SetLineWidth(3);
   		gr_RvsX_StrSuppr->Draw("lsame");
		leg->AddEntry(gr_RvsX_StrSuppr,	Name[0].Data(),		"l" );

		// gr_RvsX_WeakSuppr->SetFillColor(colors[1]);
   		// gr_RvsX_WeakSuppr->SetFillStyle(3005);
   		gr_RvsX_WeakSuppr->SetLineColor(colors[1]);
   		gr_RvsX_WeakSuppr->SetLineStyle(9);
   		gr_RvsX_WeakSuppr->SetLineWidth(3);
   		gr_RvsX_WeakSuppr->Draw("lsame");
		leg->AddEntry(gr_RvsX_WeakSuppr,	Name[1].Data(),		"l" );
	}
}


void drawGG(TString Case, TLegend * leg)
{
	std::vector<TString> Name = {"GG-hs"};
	std::vector<int> colors = {2,8,37};

	std::map<TString, std::vector<double>> Map ={
		{"Xs",			{}},
		{"Ws",			{}},
		{"Sigmas",		{}},
		{"R",			{}},
	};

	SingleDataReader SDR1("./inputfiles/GG-hs_Jpsi_Sigma.txt", 
						{"X", "Sigmas"});

	Map.at("Xs")	 = SDR1.GetVec("X");	
	Map.at("Sigmas") = SDR1.GetVec("Sigmas");
	int N = Map.at("Xs").size()-4;	
	Map.at("Xs").resize(N);	Map.at("Sigmas").resize(N);
	
	for(int j=0; j<Map.at("Xs").size(); j++ )
	{
		Map.at("Ws")			.push_back(sqrt(pow(mJpsi_PDG,2)/(Map.at("Xs")[j])));
		Map.at("Sigmas")[j] *=	0.001;
	}

	getImpulseApprox(Map);
	for (int i = 0; i < N; ++i)
	{
		Map.at("R")		.push_back(sqrt(Map.at("Sigmas")[i]/Map.at("Sigmas_IA")[i]));
	}

	auto gr_SigmasVsW 	= new TGraph(Map.at("Ws").size(),	
								Map.at("Ws").data(),			Map.at("Sigmas").data());
	auto gr_RvsW 		= new TGraph(Map.at("Xs").size(),	
								Map.at("Xs").data(),			Map.at("R").data());

	if (Case == "Sigmas")
	{
		// gr_SigmasVsW->SetFillColor(colors[0]);
   		// gr_SigmasVsW->SetFillStyle(3005);
   		gr_SigmasVsW->SetLineColor(colors[2]);
   		gr_SigmasVsW->SetLineStyle(5);
   		gr_SigmasVsW->SetLineWidth(3);
   		gr_SigmasVsW->Draw("lsame");
		leg->AddEntry(gr_SigmasVsW,	Name[0].Data(),		"l" );

	}
	else if (Case == "R")
	{
		// gr_RvsW->SetFillColor(colors[0]);
   		// gr_RvsW->SetFillStyle(3005);
		gr_RvsW->SetLineColor(colors[2]);
		gr_RvsW->SetLineStyle(5);
		gr_RvsW->SetLineWidth(3);
   		gr_RvsW->Draw("lsame");
		leg->AddEntry(gr_RvsW,	Name[0].Data(),		"l" );
	}
}


void drawbBK(TString Case, TLegend * leg)
{
	std::vector<TString> Name	 = {"bBK_GG",	"bBK_A"};

	std::vector<int> colors = {2,8,4};

	MultiDataReader MDR1({"./inputfiles/bBK_GG_Jpsi_Sigma.txt", "./inputfiles/bBK_A_Jpsi_Sigma.txt"},
						{"Xs", "Ws", "Sigmas"});
	MultiDataReader MDR2({"./inputfiles/bBK_GG_Jpsi_Xsec.txt", "./inputfiles/bBK_A_Jpsi_Xsec.txt"},
						{"Raps", "Xsec"});

	for (int i = 0; i < MDR1.SDRs.size(); ++i)
	{
		std::map<TString, std::vector<double>> Map ={
			{"Ws",			{MDR1.GetVec(i, "Ws")}},
			{"Sigmas",		{MDR1.GetVec(i, "Sigmas")}},
			{"Raps",		{MDR2.GetVec(i, "Raps")}},
			{"Xsec",		{MDR2.GetVec(i, "Xsec")}},
			{"Xs",			{MDR1.GetVec(i, "Xs")}},
			{"R",			{}}
		};

		getImpulseApprox(Map);
		for (int j = 0; j < Map.at("Ws").size(); ++j)
		{
			Map.at("Sigmas")[j] *= 0.001;
			Map.at("R")		.push_back(sqrt(Map.at("Sigmas")[j]/Map.at("Sigmas_IA")[j]));
		}

		auto gr_AnAn 	= new TGraph(Map.at("Raps").size(),	Map.at("Raps").data(),	Map.at("Xsec").data());
		auto gr_Simgas 	= new TGraph(Map.at("Ws").size(),	Map.at("Ws").data(),	Map.at("Sigmas").data());
		auto gr_R 		= new TGraph(Map.at("Xs").size(),	Map.at("Xs").data(),	Map.at("R").data());
	
		if (Case == "Xsec"){
			gr_AnAn->SetLineColor(colors[i]);
            gr_AnAn->SetLineStyle(5);
            gr_AnAn->SetLineWidth(3);
            gr_AnAn->Draw("same");
			leg->AddEntry(gr_AnAn,	Name[i].Data(),		"l" );
		}
		else if (Case == "Sigmas"){
			gr_Simgas->SetLineColor(colors[i]);
            gr_Simgas->SetLineStyle(5);
            gr_Simgas->SetLineWidth(3);
            gr_Simgas->Draw("same");
			leg->AddEntry(gr_Simgas,	Name[i].Data(),		"l" );
		}
		else if (Case == "R"){
			gr_R->SetLineColor(colors[i]);
            gr_R->SetLineStyle(5);
            gr_R->SetLineWidth(3);
            gr_R->Draw("same");
			leg->AddEntry(gr_R,	Name[i].Data(),		"l" );
		}
	}
}


void drawCD(TString Case, TLegend * leg)
{
	std::vector<TString> Name	 = {"CD_BGK",	"CD_GBW",	"CD_IIM"};
	std::vector<int> colors = {2,8,4};

	MultiDataReader MDR1({"./inputfiles/CD_BGK_Jpsi_Sigma.dat", "./inputfiles/CD_GBW_Jpsi_Sigma.dat", "./inputfiles/CD_IIM_Jpsi_Sigma.dat"},
						{"Ws", "Sigmas"});
	MultiDataReader MDR2({"./inputfiles/CD_BGK_Jpsi_Xsec.dat", "./inputfiles/CD_GBW_Jpsi_Xsec.dat", "./inputfiles/CD_IIM_Jpsi_Xsec.dat"},
						{"Raps", "Xsec"});

	for (int i = 0; i < MDR1.SDRs.size(); ++i)
	{
		std::map<TString, std::vector<double>> Map ={
			{"Ws",			{MDR1.GetVec(i, "Ws")}},
			{"Sigmas",		{MDR1.GetVec(i, "Sigmas")}},
			{"Raps",		{MDR2.GetVec(i, "Raps")}},
			{"Xsec",		{MDR2.GetVec(i, "Xsec")}},
			{"Xs",			{}},
			{"R",			{}}
		};
		getImpulseApprox(Map);
		for (int j = 0; j < Map.at("Ws").size(); ++j)
		{
			Map.at("Sigmas")[j] *= 0.001;
			Map.at("Xs")	.push_back(pow(mJpsi_PDG,2)/pow(Map.at("Ws")[j],2));
			Map.at("R")		.push_back(sqrt(Map.at("Sigmas")[j]/Map.at("Sigmas_IA")[j]));
		}
		for (int j = 0; j < Map.at("Raps").size(); ++j)
		{
			Map.at("Xsec")[j] *= 0.001;
		}

		auto gr_AnAn	= new TGraph(Map.at("Raps").size(),	Map.at("Raps").data(),	Map.at("Xsec").data());
		auto gr_Simgas 	= new TGraph(Map.at("Ws").size(),	Map.at("Ws").data(),	Map.at("Sigmas").data());
		auto gr_R 		= new TGraph(Map.at("Xs").size(),	Map.at("Xs").data(),	Map.at("R").data());
	
		if (Case == "Xsec"){
			gr_AnAn->SetLineColor(colors[i]);
            gr_AnAn->SetLineStyle(6);
            gr_AnAn->SetLineWidth(3);
            gr_AnAn->Draw("same");
			leg->AddEntry(gr_AnAn,	Name[i].Data(),		"l" );
		}
		else if (Case == "Sigmas"){
			gr_Simgas->SetLineColor(colors[i]);
            gr_Simgas->SetLineStyle(6);
            gr_Simgas->SetLineWidth(3);
            gr_Simgas->Draw("same");
			leg->AddEntry(gr_Simgas,	Name[i].Data(),		"l" );
		}
		else if (Case == "R"){
			gr_R->SetLineColor(colors[i]);
            gr_R->SetLineStyle(6);
            gr_R->SetLineWidth(3);
            gr_R->Draw("same");
			leg->AddEntry(gr_R,	Name[i].Data(),		"l" );
		}
	}
}