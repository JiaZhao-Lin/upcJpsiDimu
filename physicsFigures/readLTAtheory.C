std::vector< std::map<TString, std::vector<double>> > readLTAtheory() 
{
	std::vector<TString> fileName = {"LTA_Jpsi_weak_shadowing", "LTA_Jpsi_strong_shadowing", "EPS09_central_Jpsi"};
	std::vector< std::map<TString, std::vector<double>> > Maps;
	
	for (int i = 0; i < fileName.size(); ++i)
	{
		std::map<TString, std::vector<double>> Map ={
			{"Rap",			{}},
			{"Xsec_AnAn",	{}},
			{"Xsec_0n0n",	{}},
			{"Xsec_0nXnSum",{}},
			{"Xsec_XnXn",	{}},
		};

		TTree *tree = new TTree(); 
		tree->ReadFile(Form("./inputfiles/%s.dat",fileName[i].Data()), "y:AnAn:0n0n:0nXnSum:XnXn");

		float fRap,fXsec_AnAn,fXsec_0n0n,fXsec_0nXnSum,fXsec_XnXn;

		tree->SetBranchAddress("y", 		&fRap);
		tree->SetBranchAddress("AnAn",		&fXsec_AnAn);
		tree->SetBranchAddress("0n0n",		&fXsec_0n0n);
		tree->SetBranchAddress("0nXnSum",	&fXsec_0nXnSum);
		tree->SetBranchAddress("XnXn",		&fXsec_XnXn);

		int N = tree->GetEntries();

		for(int j=0; j<N; j++ )
		{
			tree->GetEntry(j);
			Map.at("Rap")			.push_back(double(fRap));
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
	std::vector<TString> Name = {"LTA weak shadowing", "LTA strong shadowing", "EPS09 central"};
	std::vector<int> colors = {2,8,4};
	auto LTAs = readLTAtheory();

	for (int i = 0; i < LTAs.size(); ++i)
	{
		auto Map = LTAs[i];
		auto gr_AnAn = new TGraph(Map.at("Rap").size(),	Map.at("Rap").data(),	Map.at("Xsec_AnAn").data());
		auto gr_0n0n = new TGraph(Map.at("Rap").size(),	Map.at("Rap").data(),	Map.at("Xsec_0n0n").data());
		auto gr_0nXnSum = new TGraph(Map.at("Rap").size(),	Map.at("Rap").data(),	Map.at("Xsec_0nXnSum").data());
		auto gr_XnXn = new TGraph(Map.at("Rap").size(),	Map.at("Rap").data(),	Map.at("Xsec_XnXn").data());
	
		if (Case == "AnAn"){
			gr_AnAn->Draw("same");
			gr_AnAn->SetLineColor(colors[i]);
			leg->AddEntry(gr_AnAn,	Name[i].Data(),		"l" );
		}
		else if (Case == "0n0n"){
			gr_0n0n->Draw("same");
			gr_0n0n->SetLineColor(colors[i]);
			leg->AddEntry(gr_0n0n,	Name[i].Data(),		"l" );
		}
		else if (Case == "0nXnSum"){
			gr_0nXnSum->Draw("same");
			gr_0nXnSum->SetLineColor(colors[i]);
			leg->AddEntry(gr_0nXnSum,	Name[i].Data(),	"l" );
		}
		else if (Case == "XnXn"){
			gr_XnXn->Draw("same");
			gr_XnXn->SetLineColor(colors[i]);
			leg->AddEntry(gr_XnXn,	Name[i].Data(),		"l" );
		}
	}
}

std::vector< std::map<TString, std::vector<double>> > readLTAtheory_Sigmas_R() 
{
	std::vector<TString> fileName = {"LTA_Jpsi_Sigma_R"};
	std::vector< std::map<TString, std::vector<double>> > Maps;
	
	for (int i = 0; i < fileName.size(); ++i)
	{
		std::map<TString, std::vector<double>> Map ={
			{"Xs",			{}},
			{"R",			{}},
			{"R_High",		{}},
			{"R_Low",		{}},
			{"Ws",			{}},
			{"Sigmas",		{}},
			{"Sigmas_High",	{}},
			{"Sigmas_Low",	{}}
		};

		TTree *tree = new TTree(); 
		tree->ReadFile(Form("./inputfiles/%s.dat",fileName[i].Data()), "Xs:R_High:R_Low:Ws:Sigmas_High:Sigmas_Low");

		float fXs,fR_Low,fR_High,fWs, fSigmas_Low,fSigmas_High;

		tree->SetBranchAddress("Xs", 			&fXs);
		tree->SetBranchAddress("R_High",		&fR_High);
		tree->SetBranchAddress("R_Low",			&fR_Low);
		tree->SetBranchAddress("Ws",			&fWs);
		tree->SetBranchAddress("Sigmas_High",	&fSigmas_High);
		tree->SetBranchAddress("Sigmas_Low",	&fSigmas_Low);

		int N = tree->GetEntries();

		for(int j=0; j<N; j++ )
		{
			tree->GetEntry(j);
			Map.at("Xs")			.push_back(double(fXs));
			Map.at("R")				.push_back(double((fR_High+fR_Low)/2));
			Map.at("R_High")		.push_back(double(fR_High));
			Map.at("R_Low")			.push_back(double(fR_Low));
			Map.at("Ws")			.push_back(double(fWs));
			Map.at("Sigmas")		.push_back(double((fSigmas_High+fSigmas_Low)/2));
			Map.at("Sigmas_High")	.push_back(double(fSigmas_High));
			Map.at("Sigmas_Low")	.push_back(double(fSigmas_Low));
		}

		Maps.push_back(Map);
	}
	return Maps;
}

void drawLTA_Sigmas_R(TString Case, TLegend * leg)
{
	std::vector<TString> Name = {"LTA"};
	std::vector<int> colors = {2,3,4};
	auto LTAs = readLTAtheory_Sigmas_R();

	for (int i = 0; i < LTAs.size(); ++i)
	{
		auto Map = LTAs[i];
		auto gr_SigmasVsW = new TGraphAsymmErrors(Map.at("Ws").size(),	
													Map.at("Ws").data(),			Map.at("Sigmas").data(),
													0,								0,
													0,								0);
		auto gr_RvsX	 = new TGraphAsymmErrors(Map.at("Xs").size(),	
													Map.at("Xs").data(),			Map.at("R").data(),
													0,								0,
													0,								0);

		if (Case == "Sigmas")
		{
			// gr_SigmasVsW->SetFillColor(colors[i]);
	  //  		gr_SigmasVsW->SetFillStyle(3005);
	   		gr_SigmasVsW->Draw("same");
			leg->AddEntry(gr_SigmasVsW,	Name[i].Data(),		"l" );
		}
		else if (Case == "R")
		{
			// gr_RvsX->SetFillColor(colors[i]);
	  //  		gr_RvsX->SetFillStyle(3005);
	   		gr_RvsX->Draw("same");
			leg->AddEntry(gr_RvsX,	Name[i].Data(),		"l" );
		}
	}
}