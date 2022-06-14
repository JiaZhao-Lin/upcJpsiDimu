#include "./ImpulseApproximation/ImpulseApprox.C"

std::vector< std::map<TString, std::vector<double>> > readLTAtheory() 
{
	std::vector<TString> fileName = {"LTA_Jpsi_weak_shadowing", "LTA_Jpsi_strong_shadowing", "EPS09_central_Jpsi"};
	std::vector< std::map<TString, std::vector<double>> > Maps;
	
	for (int i = 0; i < fileName.size(); ++i)
	{
		std::map<TString, std::vector<double>> Map ={
			{"Raps",			{}},
			{"Xsec_AnAn",	{}},
			{"Xsec_0n0n",	{}},
			{"Xsec_0nXnSum",{}},
			{"Xsec_XnXn",	{}},
		};

		TTree *tree = new TTree(); 
		tree->ReadFile(Form("./inputfiles/%s.dat",fileName[i].Data()), "y:AnAn:0n0n:0nXnSum:XnXn");

		float fRaps,fXsec_AnAn,fXsec_0n0n,fXsec_0nXnSum,fXsec_XnXn;

		tree->SetBranchAddress("y", 		&fRaps);
		tree->SetBranchAddress("AnAn",		&fXsec_AnAn);
		tree->SetBranchAddress("0n0n",		&fXsec_0n0n);
		tree->SetBranchAddress("0nXnSum",	&fXsec_0nXnSum);
		tree->SetBranchAddress("XnXn",		&fXsec_XnXn);

		int N = tree->GetEntries();

		for(int j=0; j<N; j++ )
		{
			tree->GetEntry(j);
			Map.at("Raps")			.push_back(double(fRaps));
			Map.at("Xsec_AnAn")		.push_back(double(fXsec_AnAn));
			Map.at("Xsec_0n0n")		.push_back(double(fXsec_0n0n));
			Map.at("Xsec_0nXnSum")	.push_back(double(fXsec_0nXnSum));
			Map.at("Xsec_XnXn")		.push_back(double(fXsec_XnXn));
		}

		Maps.push_back(Map);
	}
	return Maps;
}

void drawLTA(TString Case, TLegend * leg)
{
	std::vector<TString> Name = {"LTA WS", "LTA SS", "EPS09 central"};
	std::vector<int> colors = {2,8,4};
	auto LTAs = readLTAtheory();

	for (int i = 0; i < LTAs.size()-1; ++i)
	{
		auto Map = LTAs[i];
		auto gr_AnAn = new TGraph(Map.at("Raps").size(),	Map.at("Raps").data(),	Map.at("Xsec_AnAn").data());
		auto gr_0n0n = new TGraph(Map.at("Raps").size(),	Map.at("Raps").data(),	Map.at("Xsec_0n0n").data());
		auto gr_0nXnSum = new TGraph(Map.at("Raps").size(),	Map.at("Raps").data(),	Map.at("Xsec_0nXnSum").data());
		auto gr_XnXn = new TGraph(Map.at("Raps").size(),	Map.at("Raps").data(),	Map.at("Xsec_XnXn").data());
	
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

std::map<TString, std::vector<double>> readLTAtheory_Sigmas_R() 
{
	std::vector<TString> fileName = {"LTA_Jpsi_Sigma_R"};
	std::vector< std::map<TString, std::vector<double>> > Maps;
	

	std::map<TString, std::vector<double>> Map ={
		{"Xs",				{}},
		{"R_WeakSuppr",		{}},
		{"R_StrSuppr",		{}},
		{"Ws",				{}},
		{"Sigmas_WeakSuppr",{}},
		{"Sigmas_StrSuppr",	{}}
	};


	TTree *tree = new TTree(); 
	tree->ReadFile(Form("./inputfiles/%s_1.dat",fileName[0].Data()), "y:Xs_Left:R_StrSuppr_Left:R_WeakSuppr_Left:Ws_Left:Sigmas_StrSuppr_Left:Sigmas_WeakSuppr_Left:Xs_Right:R_StrSuppr_Right:R_WeakSuppr_Right:Ws_Right:Sigmas_StrSuppr_Right:Sigmas_WeakSuppr_Right");
	float fy, fXs_Left, fR_StrSuppr_Left, fR_WeakSuppr_Left, fWs_Left, fSigmas_StrSuppr_Left, fSigmas_WeakSuppr_Left, fXs_Right, fR_StrSuppr_Right, fR_WeakSuppr_Right, fWs_Right, fSigmas_StrSuppr_Right, fSigmas_WeakSuppr_Right;

	tree->SetBranchAddress("Xs_Left", 				&fXs_Left);
	tree->SetBranchAddress("R_StrSuppr_Left",		&fR_StrSuppr_Left);
	tree->SetBranchAddress("R_WeakSuppr_Left",		&fR_WeakSuppr_Left);
	tree->SetBranchAddress("Xs_Right",				&fXs_Right);
	tree->SetBranchAddress("R_StrSuppr_Right",		&fR_StrSuppr_Right);
	tree->SetBranchAddress("R_WeakSuppr_Right",		&fR_WeakSuppr_Right);

	int N = tree->GetEntries();
	for(int j=0; j<N; j++ )
	{
		tree->GetEntry(j);
		Map.at("Xs")			.push_back(double(fXs_Right));
		Map.at("R_StrSuppr")	.push_back(double(fR_StrSuppr_Right));
		Map.at("R_WeakSuppr")	.push_back(double(fR_WeakSuppr_Right));
	}
	for(int j=0; j<N; j++ )
	{
		tree->GetEntry(j);
		Map.at("Xs")			.insert(Map.at("Xs").begin(),			double(fXs_Left));
		Map.at("R_StrSuppr")	.insert(Map.at("R_StrSuppr").begin(),	double(fR_StrSuppr_Left));
		Map.at("R_WeakSuppr")	.insert(Map.at("R_WeakSuppr").begin(),	double(fR_WeakSuppr_Left));
	}


	TTree *tree2 = new TTree(); 
	tree2->ReadFile(Form("./inputfiles/%s_2.dat",fileName[0].Data()), "y:Ws_Left:Sigmas_StrSuppr_Left:Sigmas_WeakSuppr_Left:Ws_Right:Sigmas_StrSuppr_Right:Sigmas_WeakSuppr_Right");

	// float fy,fXs,fR_WeakSuppr_Left,fR_WeakSuppr_Right,
	// 	fWs, fSigmas_WeakSuppr_Left,fSigma_WeakSupprs_Right;

	tree2->SetBranchAddress("Ws_Left", 					&fWs_Left);
	tree2->SetBranchAddress("Sigmas_StrSuppr_Left",		&fSigmas_StrSuppr_Left);
	tree2->SetBranchAddress("Sigmas_WeakSuppr_Left",	&fSigmas_WeakSuppr_Left);
	tree2->SetBranchAddress("Ws_Right",					&fWs_Right);
	tree2->SetBranchAddress("Sigmas_StrSuppr_Right",	&fSigmas_StrSuppr_Right);
	tree2->SetBranchAddress("Sigmas_WeakSuppr_Right",	&fSigmas_WeakSuppr_Right);

	N = tree2->GetEntries();
	for(int j=0; j<N; j++ )
	{
		tree2->GetEntry(j);
		Map.at("Ws")			.push_back(double(fWs_Right));
		Map.at("Sigmas_StrSuppr")	.push_back(double(fSigmas_StrSuppr_Right));
		Map.at("Sigmas_WeakSuppr")	.push_back(double(fSigmas_WeakSuppr_Right));
	}
	for(int j=0; j<N; j++ )
	{
		tree2->GetEntry(j);
		Map.at("Ws")				.insert(Map.at("Ws").begin(),			double(fWs_Left));
		Map.at("Sigmas_StrSuppr")	.insert(Map.at("Sigmas_StrSuppr").begin(),	double(fSigmas_StrSuppr_Left));
		Map.at("Sigmas_WeakSuppr")	.insert(Map.at("Sigmas_WeakSuppr").begin(),	double(fSigmas_WeakSuppr_Left));
	}

	return Map;
}

void drawLTA_Sigmas_R(TString Case, TLegend * leg)
{
	std::vector<TString> Name = {"LTA_SS","LTA_WS"};
	std::vector<int> colors = {2,8,4};
	auto LTAs = readLTAtheory_Sigmas_R();

	auto Map = LTAs;
	auto gr_SigmasVsW_StrSuppr = new TGraph(Map.at("Ws").size(),	
												Map.at("Ws").data(),			Map.at("Sigmas_StrSuppr").data());
	auto gr_SigmasVsW_WeakSuppr = new TGraph(Map.at("Ws").size(),	
												Map.at("Ws").data(),			Map.at("Sigmas_WeakSuppr").data());
	auto gr_RvsX_StrSuppr	 = new TGraph(Map.at("Xs").size(),	
												Map.at("Xs").data(),			Map.at("R_StrSuppr").data());
	auto gr_RvsX_WeakSuppr	 = new TGraph(Map.at("Xs").size(),	
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

std::map<TString, std::vector<double>> readGGtheory() 
{
	const double JpsiMass   = 3.096916;

	std::vector<TString> fileName = {"GG-hs_Jpsi_Sigma"};

	std::map<TString, std::vector<double>> Map ={
		{"Xs",			{}},
		{"Ws",			{}},
		{"Sigmas",		{}},
		{"R",			{}},
	};

	TTree *tree = new TTree(); 
	tree->ReadFile(Form("./inputfiles/%s.txt",fileName[0].Data()), "X:Sigmas");

	float fXs, fSimgas;

	tree->SetBranchAddress("X", 		&fXs);
	tree->SetBranchAddress("Sigmas",	&fSimgas);

	int N = tree->GetEntries()-4;

	for(int j=0; j<N; j++ )
	{
		tree->GetEntry(j);
		Map.at("Xs")			.push_back(double(fXs));
		Map.at("Ws")			.push_back(sqrt(pow(JpsiMass,2)/(double(fXs))));
		Map.at("Sigmas")		.push_back(double(fSimgas)*0.001);
	}

	std::map<TString, std::vector<double>> Temp_Map;
	getImpulseApprox(Map);
	for (int i = 0; i < N; ++i)
	{
		Map.at("R")		.push_back(sqrt(Map.at("Sigmas")[i]/Map.at("Sigmas_IA")[i]));
	}

	return Map;
}

void drawGG(TString Case, TLegend * leg)
{
	std::vector<TString> Name = {"GG-hs"};
	std::vector<int> colors = {2,8,37};
	auto LTAs = readGGtheory();

	auto Map = LTAs;
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

std::vector< std::map<TString, std::vector<double>> > readbBKtheory() 
{
	std::vector<TString> fileName_Sigma	 = {"bBK_GG_Jpsi_Sigma",	"bBK_A_Jpsi_Sigma"};
	std::vector<TString> fileName_Xsec	 = {"bBK_GG_Jpsi_Xsec",		"bBK_A_Jpsi_Xsec"};

	std::vector< std::map<TString, std::vector<double>> > Maps;
	
	for (int i = 0; i < fileName_Sigma.size(); ++i)
	{
		std::map<TString, std::vector<double>> Map ={
			{"Raps",			{}},
			{"Xsec_AnAn",	{}},
			{"Xs",			{}},
			{"Ws",			{}},
			{"Sigmas",		{}},
			{"R",			{}}
		};

		TTree *tree = new TTree(); 
		tree->ReadFile(Form("./inputfiles/%s.txt",fileName_Sigma[i].Data()), "X:W:Sigma");

		TTree *tree1 = new TTree(); 
		tree1->ReadFile(Form("./inputfiles/%s.txt",fileName_Xsec[i].Data()), "Raps:Xsec");

		float fXs,fWs,fSigmas,fRaps,fXsec;

		tree->SetBranchAddress("X", 		&fXs);
		tree->SetBranchAddress("W",			&fWs);
		tree->SetBranchAddress("Sigma",		&fSigmas);

		tree1->SetBranchAddress("Raps", 	&fRaps);
		tree1->SetBranchAddress("Xsec",		&fXsec);

		int N = tree->GetEntries();
		for(int j=0; j<N; j++ )
		{
			tree->GetEntry(j);
			Map.at("Xs")			.push_back(double(fXs));
			Map.at("Ws")			.push_back(double(fWs));
			Map.at("Sigmas")		.push_back(double(fSigmas)*0.001);
		}

		getImpulseApprox(Map);
		for (int i = 0; i < N; ++i)
		{
			Map.at("R")		.push_back(sqrt(Map.at("Sigmas")[i]/Map.at("Sigmas_IA")[i]));
		}

		N = tree1->GetEntries();
		for(int j=0; j<N; j++ )
		{
			tree1->GetEntry(j);
			Map.at("Raps")			.push_back(double(fRaps));
			Map.at("Xsec_AnAn")		.push_back(double(fXsec));
		}

		Maps.push_back(Map);
	}
	return Maps;
}

void drawbBK(TString Case, TLegend * leg)
{
	std::vector<TString> Name	 = {"bBK_GG",	"bBK_A"};

	std::vector<int> colors = {2,8,4};
	auto bBKs = readbBKtheory();

	for (int i = 0; i < bBKs.size(); ++i)
	{
		auto Map = bBKs[i];
		auto gr_AnAn = new TGraph(Map.at("Raps").size(),	Map.at("Raps").data(),	Map.at("Xsec_AnAn").data());
		auto gr_Simgas = new TGraph(Map.at("Ws").size(),	Map.at("Ws").data(),	Map.at("Sigmas").data());
		auto gr_R = new TGraph(Map.at("Xs").size(),	Map.at("Xs").data(),	Map.at("R").data());
	
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

std::vector< std::map<TString, std::vector<double>> > readCDtheory() 
{
	const double JpsiMass   = 3.096916;
	std::vector<TString> fileName_Sigma	 = {"CD_BGK_Jpsi_Sigma",	"CD_GBW_Jpsi_Sigma",	"CD_IIM_Jpsi_Sigma"};
	std::vector<TString> fileName_Xsec	 = {"CD_BGK_Jpsi_Xsec",		"CD_GBW_Jpsi_Xsec",		"CD_IIM_Jpsi_Xsec"};

	std::vector< std::map<TString, std::vector<double>> > Maps;
	
	for (int i = 0; i < fileName_Sigma.size(); ++i)
	{
		std::map<TString, std::vector<double>> Map ={
			{"Raps",			{}},
			{"Xsec_AnAn",	{}},
			{"Xs",			{}},
			{"Ws",			{}},
			{"Sigmas",		{}},
			{"R",			{}}
		};

		TTree *tree = new TTree(); 
		tree->ReadFile(Form("./inputfiles/%s.dat",fileName_Sigma[i].Data()), "W:Sigma");

		TTree *tree1 = new TTree(); 
		tree1->ReadFile(Form("./inputfiles/%s.dat",fileName_Xsec[i].Data()), "Raps:Xsec");

		float fWs,fSigmas,fRaps,fXsec;

		tree->SetBranchAddress("W",			&fWs);
		tree->SetBranchAddress("Sigma",		&fSigmas);

		tree1->SetBranchAddress("Raps", 	&fRaps);
		tree1->SetBranchAddress("Xsec",		&fXsec);

		int N = tree->GetEntries();
		for(int j=0; j<N; j++ )
		{
			tree->GetEntry(j);
			Map.at("Ws")			.push_back(double(fWs));
			Map.at("Xs")			.push_back(pow(JpsiMass,2)/pow(Map.at("Ws")[j],2));
			Map.at("Sigmas")		.push_back(double(fSigmas)*0.001);
		}

		getImpulseApprox(Map);
		for (int i = 0; i < N; ++i)
		{
			Map.at("R")		.push_back(sqrt(Map.at("Sigmas")[i]/Map.at("Sigmas_IA")[i]));
		}

		N = tree1->GetEntries();
		for(int j=0; j<N; j++ )
		{
			tree1->GetEntry(j);
			Map.at("Raps")			.push_back(double(fRaps));
			Map.at("Xsec_AnAn")		.push_back(double(fXsec)*0.001);
		}

		Maps.push_back(Map);
	}
	return Maps;
}

void drawCD(TString Case, TLegend * leg)
{
	std::vector<TString> Name	 = {"CD_BGK",	"CD_GBW",	"CD_IIM"};

	std::vector<int> colors = {2,8,4};
	auto CDs = readCDtheory();

	for (int i = 0; i < CDs.size(); ++i)
	{
		auto Map = CDs[i];
		auto gr_AnAn = new TGraph(Map.at("Raps").size(),	Map.at("Raps").data(),	Map.at("Xsec_AnAn").data());
		auto gr_Simgas = new TGraph(Map.at("Ws").size(),	Map.at("Ws").data(),	Map.at("Sigmas").data());
		auto gr_R = new TGraph(Map.at("Xs").size(),	Map.at("Xs").data(),	Map.at("R").data());
	
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