#include "./ImpulseApprox.C"
#include "../common/DataReader.C"
#include "./ParamConverter.C"

enum class TheoryList
{
	Sigma_R_CGC,
	Sigma_IA,
	DSigmaDy_LTA,
	Sigma_R_LTA,
	Sigma_R_GG,
	DSigmaDy_Sigma_R_bBK,
	DSigmaDy_Sigma_R_CD
};

struct Theory
{
	TString Case;
	Theory(TString Case_) : Case{Case_} {}
	virtual void Draw(TLegend * leg) = 0;
	virtual ~Theory() = default;
};

struct Sigma_IA : Theory
{
	TGraph* gr_Sigma_IA	{nullptr};

	Sigma_IA(TString Case_ ) : Theory{Case_} {}
	virtual ~Sigma_IA() override
	{
		delete gr_Sigma_IA;
	}

	virtual void Draw(TLegend * leg) override
	{
		std::vector<double> temp_W;	// For plotting IA
		for (int i = 5; i < 500; ++i){	temp_W.push_back(i);	}
		auto IA_pair = ImpulseApprox::GetImpulseApprox(temp_W);

		gr_Sigma_IA	= new TGraphErrors(temp_W.size(),	temp_W.data(),	IA_pair.first.data());
		gr_Sigma_IA	->SetMarkerColor(2);
		gr_Sigma_IA	->SetLineStyle(5);
		gr_Sigma_IA	->SetLineColor(1);
		gr_Sigma_IA	->SetLineWidth(2);
		gr_Sigma_IA	->Draw("lsame");

		leg->AddEntry(gr_Sigma_IA,	"Impulse Approximation",	"l");
	}
};

struct Sigma_R_CGC : Theory
{
	TGraph * gr_SigmaVsW	{nullptr};
	TGraph * gr_RvsX		{nullptr};
	
	Sigma_R_CGC(TString Case_ ) : Theory(Case_) {}
	virtual ~Sigma_R_CGC()	override
	{
		delete gr_SigmaVsW;
		delete gr_RvsX;
	}

	virtual void Draw(TLegend * leg) override
	{
		std::vector<TString> Name = {"CGC IPsat"};
		std::vector<int> colors = {11};
		std::vector<int> styles = {5};

		std::map<TString, std::vector<double>> Map ={
			{"X",			{}},
			{"W",			{}},
			{"Sigma",		{}},
			{"R",			{}},
		};

		SingleDataReader SDR1("./inFiles/CGC_Jpsi_Sigma_NoFluct.txt", 
							{"W", "Sigma", "Sigma_InCoh"});

		Map.at("W")	 = SDR1.GetVec("W");	
		Map.at("Sigma") = SDR1.GetVec("Sigma");
		
		for(int j=0; j<Map.at("W").size(); j++ )
		{
			Map.at("X")			.push_back(ParamConverter::W2x(Map.at("W")[j]));
		}

		auto IA_pair = ImpulseApprox::GetImpulseApprox(Map.at("W"));
		for (int i = 0; i < Map.at("Sigma").size(); ++i)
		{
			Map.at("R")		.push_back(sqrt(Map.at("Sigma")[i]/IA_pair.first[i]));
		}

		auto gr_SigmaVsW 	= new TGraph(Map.at("W").size(),	
									Map.at("W").data(),			Map.at("Sigma").data());
		auto gr_RvsW 		= new TGraph(Map.at("X").size(),	
									Map.at("X").data(),			Map.at("R").data());

		if (Case == "Sigma")
		{
			// gr_SigmaVsW->SetFillColor(colors[0]);
			// gr_SigmaVsW->SetFillStyle(3005);
			
			gr_SigmaVsW->SetLineColor(colors[0]);
			gr_SigmaVsW->SetLineStyle(styles[0]);
			gr_SigmaVsW->SetLineWidth(2);
			gr_SigmaVsW->Draw("lsame");
			leg->AddEntry(gr_SigmaVsW,	Name[0].Data(),		"l" );

		}
		else if (Case == "R")
		{
			// gr_RvsW->SetFillColor(colors[0]);
			// gr_RvsW->SetFillStyle(3005);
			gr_RvsW->SetLineColor(colors[0]);
			gr_RvsW->SetLineStyle(styles[0]);
			gr_RvsW->SetLineWidth(2);
			gr_RvsW->Draw("lsame");
			leg->AddEntry(gr_RvsW,	Name[0].Data(),		"l" );
		}
	}
};


struct DSigmaDy_LTA : Theory
{
	TGraph* gr_AnAn		{nullptr};
	TGraph* gr_0n0n		{nullptr};
	TGraph* gr_0nXnSum	{nullptr};
	TGraph* gr_XnXn		{nullptr};

	DSigmaDy_LTA(TString Case_ ) : Theory{Case_} {}
	virtual ~DSigmaDy_LTA() override
	{
		delete gr_AnAn;
		delete gr_0n0n;
		delete gr_0nXnSum;
		delete gr_XnXn;
	}

	virtual void Draw(TLegend * leg) override
	{
		std::vector<TString> Name = {"LTA_SS", "LTA_WS", "EPS09 central"};
		std::vector<int> colors = {1,2,1};
		std::vector<int> styles = {2,5,3};

		MultiDataReader MDR({"./inFiles/LTA_Jpsi_strong_shadowing.dat", "./inFiles/LTA_Jpsi_weak_shadowing.dat", "./inFiles/EPS09_central_Jpsi.dat"},
							{"Dy", "AnAn", "0n0n", "0nXnSum", "XnXn"});

		for (int i = 0; i < MDR.SDRs.size()-1; ++i)
		{
			gr_AnAn 	= new TGraph(MDR.GetVec(i, "Dy").size(),	MDR.GetVec(i, "Dy").data(),	MDR.GetVec(i, "AnAn").data());
			gr_0n0n 	= new TGraph(MDR.GetVec(i, "Dy").size(),	MDR.GetVec(i, "Dy").data(),	MDR.GetVec(i, "0n0n").data());
			gr_0nXnSum 	= new TGraph(MDR.GetVec(i, "Dy").size(),	MDR.GetVec(i, "Dy").data(),	MDR.GetVec(i, "0nXnSum").data());
			gr_XnXn 	= new TGraph(MDR.GetVec(i, "Dy").size(),	MDR.GetVec(i, "Dy").data(),	MDR.GetVec(i, "XnXn").data());
		
			if (Case == "AnAn")
			{
				gr_AnAn->SetLineColor(colors[i]);
				gr_AnAn->SetLineStyle(styles[i]);
				gr_AnAn->SetLineWidth(2);
				gr_AnAn->Draw("lsame");
				leg->AddEntry(gr_AnAn,	Name[i].Data(),		"l" );
			}
			else if (Case == "0n0n")
			{
				gr_0n0n->SetLineColor(1);
				gr_0n0n->SetLineStyle(styles[i]);
				gr_0n0n->SetLineWidth(2);
				gr_0n0n->Draw("lsame");
				leg->AddEntry(gr_0n0n,	" ",		"l" );
				//leg->AddEntry(gr_0n0n,	Name[i].Data(),		"l" );
			}
			else if (Case == "0nXnSum")
			{
				gr_0nXnSum->SetLineColor(2);
				gr_0nXnSum->SetLineStyle(styles[i]);
				gr_0nXnSum->SetLineWidth(2);
				gr_0nXnSum->Draw("lsame");
				leg->AddEntry(gr_0nXnSum,	" ",	"l" );
				//leg->AddEntry(gr_0nXnSum,	Name[i].Data(),	"l" );
			}
			else if (Case == "XnXn")
			{
				gr_XnXn->SetLineColor(4);
				gr_XnXn->SetLineStyle(styles[i]);
				gr_XnXn->SetLineWidth(2);
				gr_XnXn->Draw("lsame");
				leg->AddEntry(gr_XnXn,	" ",		"l" );
				//leg->AddEntry(gr_XnXn,	Name[i].Data(),		"l" );
			}
		}
	}
};

struct Sigma_R_LTA : Theory
{
	TGraph* gr_SigmaVsW_StrSuppr	{nullptr};
	TGraph* gr_SigmaVsW_WeakSuppr	{nullptr};
	TGraph* gr_RvsX_StrSuppr		{nullptr};
	TGraph* gr_RvsX_WeakSuppr		{nullptr};

	Sigma_R_LTA(TString Case_ ) : Theory{Case_} {}
	virtual ~Sigma_R_LTA() override
	{
		delete gr_SigmaVsW_StrSuppr;
		delete gr_SigmaVsW_WeakSuppr;
		delete gr_RvsX_StrSuppr;
		delete gr_RvsX_WeakSuppr;
	}

	virtual void Draw(TLegend * leg)	override
	{
		std::vector<TString> Name = {"LTA_SS","LTA_WS"};
		std::vector<int> colors = {2,2,2};
		std::vector<int> styles = {1,7,3};

		SingleDataReader SDR1("./inFiles/LTA_Jpsi_Sigma_R_1.dat", 
							{"Dy", "X_Left", "R_StrSuppr_Left", "R_WeakSuppr_Left", "W_Left", "Sigma_StrSuppr_Left", "Sigma_WeakSuppr_Left", "X_Right", "R_StrSuppr_Right", "R_WeakSuppr_Right", "W_Right", "Sigma_StrSuppr_Right", "Sigma_WeakSuppr_Right"});
		SingleDataReader SDR2("./inFiles/LTA_Jpsi_Sigma_R_2.dat",
							{"Dy", "W_Left", "Sigma_StrSuppr_Left", "Sigma_WeakSuppr_Left", "W_Right", "Sigma_StrSuppr_Right", "Sigma_WeakSuppr_Right"});

		std::map<TString, std::vector<double>> Map ={
			{"X",				{}},
			{"R_WeakSuppr",		{}},
			{"R_StrSuppr",		{}},
			{"W",				{}},
			{"Sigma_WeakSuppr",{}},
			{"Sigma_StrSuppr",	{}}
		};

		/*Using R with IA calculated from our code for now. There is 5-10% difference in calculated form factor value compares to Vadim*/
		Map.at("X")				=	SDR1.GetVec("X_Right");
		// Map.at("R_WeakSuppr")		=	SDR1.GetVec("R_WeakSuppr_Right");
		// Map.at("R_StrSuppr")		=	SDR1.GetVec("R_StrSuppr_Right");
		Map.at("W")				=	SDR2.GetVec("W_Right");
		Map.at("Sigma_WeakSuppr")	=	SDR2.GetVec("Sigma_WeakSuppr_Right");
		Map.at("Sigma_StrSuppr")	=	SDR2.GetVec("Sigma_StrSuppr_Right");
		auto X_temp 				= SDR1.GetVec("X_Left");			std::reverse(X_temp.begin(),X_temp.end());
		// auto R_WeakSuppr_temp 	= SDR1.GetVec("R_WeakSuppr_Left");	std::reverse(R_WeakSuppr_temp.begin(),R_WeakSuppr_temp.end());
		// auto R_StrSuppr_temp 	= SDR1.GetVec("R_StrSuppr_Left");	std::reverse(R_StrSuppr_temp.begin(),R_StrSuppr_temp.end());
		auto Ws_temp 				= SDR2.GetVec("W_Left");			std::reverse(Ws_temp.begin(),Ws_temp.end());
		auto Sigma_WeakSuppr_temp 	= SDR2.GetVec("Sigma_WeakSuppr_Left");	std::reverse(Sigma_WeakSuppr_temp.begin(),Sigma_WeakSuppr_temp.end());
		auto Sigma_StrSuppr_temp 	= SDR2.GetVec("Sigma_StrSuppr_Left");	std::reverse(Sigma_StrSuppr_temp.begin(),Sigma_StrSuppr_temp.end());

		Map.at("X")			.insert(Map.at("X").begin(),			X_temp.begin(),	X_temp.end());
		// Map.at("R_WeakSuppr")	.insert(Map.at("R_WeakSuppr").begin(),	R_WeakSuppr_temp.begin(),	R_WeakSuppr_temp.end());
		// Map.at("R_StrSuppr")	.insert(Map.at("R_StrSuppr").begin(),	R_StrSuppr_temp.begin(),	R_StrSuppr_temp.end());
		Map.at("W")			.insert(Map.at("W").begin(),			Ws_temp.begin(),	Ws_temp.end());
		Map.at("Sigma_WeakSuppr")	.insert(Map.at("Sigma_WeakSuppr").begin(),	Sigma_WeakSuppr_temp.begin(),	Sigma_WeakSuppr_temp.end());
		Map.at("Sigma_StrSuppr")	.insert(Map.at("Sigma_StrSuppr").begin(),	Sigma_StrSuppr_temp.begin(),	Sigma_StrSuppr_temp.end());
		
		auto IA_pair = ImpulseApprox::GetImpulseApprox(Map.at("W"));
		for (int i = 0; i < Map.at("W").size(); ++i)
		{
			Map.at("R_WeakSuppr")		.push_back(sqrt(Map.at("Sigma_WeakSuppr")[i]/IA_pair.first[i]));
			Map.at("R_StrSuppr")		.push_back(sqrt(Map.at("Sigma_StrSuppr")[i]/IA_pair.first[i]));
		}

		gr_SigmaVsW_StrSuppr 	= new TGraph(Map.at("W").size(),	
													Map.at("W").data(),			Map.at("Sigma_StrSuppr").data());
		gr_SigmaVsW_WeakSuppr 	= new TGraph(Map.at("W").size(),	
													Map.at("W").data(),			Map.at("Sigma_WeakSuppr").data());
		gr_RvsX_StrSuppr	 	= new TGraph(Map.at("X").size(),	
													Map.at("X").data(),			Map.at("R_StrSuppr").data());
		gr_RvsX_WeakSuppr	 	= new TGraph(Map.at("X").size(),	
													Map.at("X").data(),			Map.at("R_WeakSuppr").data());

		if (Case == "Sigma")
		{
			// gr_SigmaVsW_StrSuppr->SetFillColor(colors[0]);
			// gr_SigmaVsW_StrSuppr->SetFillStyle(3005);
			gr_SigmaVsW_StrSuppr->SetLineColor(colors[0]);
			gr_SigmaVsW_StrSuppr->SetLineStyle(styles[0]);
			gr_SigmaVsW_StrSuppr->SetLineWidth(2);
			gr_SigmaVsW_StrSuppr->Draw("lsame");
			leg->AddEntry(gr_SigmaVsW_StrSuppr,	Name[0].Data(),		"l" );

			// gr_SigmaVsW_WeakSuppr->SetFillColor(colors[1]);
			// gr_SigmaVsW_WeakSuppr->SetFillStyle(3005);
			gr_SigmaVsW_WeakSuppr->SetLineColor(colors[1]);
			gr_SigmaVsW_WeakSuppr->SetLineStyle(styles[1]);
			gr_SigmaVsW_WeakSuppr->SetLineWidth(2);
			gr_SigmaVsW_WeakSuppr->Draw("lsame");
			leg->AddEntry(gr_SigmaVsW_WeakSuppr,	Name[1].Data(),		"l" );
		}
		else if (Case == "R")
		{
			// gr_RvsX_StrSuppr->SetFillColor(colors[0]);
			// gr_RvsX_StrSuppr->SetFillStyle(3005);
			gr_RvsX_StrSuppr->SetLineColor(colors[0]);
			gr_RvsX_StrSuppr->SetLineStyle(styles[0]);
			gr_RvsX_StrSuppr->SetLineWidth(2);
			gr_RvsX_StrSuppr->Draw("lsame");
			leg->AddEntry(gr_RvsX_StrSuppr,	Name[0].Data(),		"l" );

			// gr_RvsX_WeakSuppr->SetFillColor(colors[1]);
			// gr_RvsX_WeakSuppr->SetFillStyle(3005);
			gr_RvsX_WeakSuppr->SetLineColor(colors[1]);
			gr_RvsX_WeakSuppr->SetLineStyle(styles[1]);
			gr_RvsX_WeakSuppr->SetLineWidth(2);
			gr_RvsX_WeakSuppr->Draw("lsame");
			leg->AddEntry(gr_RvsX_WeakSuppr,	Name[1].Data(),		"l" );
		}
	}
};

struct Sigma_R_GG : Theory
{
	TGraph * gr_SigmaVsW	{nullptr};
	TGraph * gr_RvsX		{nullptr};
	
	Sigma_R_GG(TString Case_ ) : Theory(Case_) {}
	virtual ~Sigma_R_GG()	override
	{
		delete gr_SigmaVsW;
		delete gr_RvsX;
	}

	virtual void Draw(TLegend * leg) override
	{
		std::vector<TString> Name = {"GG-hs"};
		std::vector<int> colors = {6};
		std::vector<int> styles = {5};

		std::map<TString, std::vector<double>> Map ={
			{"X",			{}},
			{"W",			{}},
			{"Sigma",		{}},
			{"R",			{}},
		};

		SingleDataReader SDR1("./inFiles/GG-hs_Jpsi_Sigma.txt", 
							{"X", "Sigma"});

		Map.at("X")	 = SDR1.GetVec("X");	
		Map.at("Sigma") = SDR1.GetVec("Sigma");
		int N = Map.at("X").size()-4;	
		Map.at("X").resize(N);	Map.at("Sigma").resize(N);
		
		for(int j=0; j<Map.at("X").size(); j++ )
		{
			Map.at("W")			.push_back(ParamConverter::x2W(Map.at("X")[j]));
			Map.at("Sigma")[j] *=	0.001;
		}

		auto IA_pair = ImpulseApprox::GetImpulseApprox(Map.at("W"));
		for (int i = 0; i < N; ++i)
		{
			Map.at("R")		.push_back(sqrt(Map.at("Sigma")[i]/IA_pair.first[i]));
		}

		auto gr_SigmaVsW 	= new TGraph(Map.at("W").size(),	
									Map.at("W").data(),			Map.at("Sigma").data());
		auto gr_RvsW 		= new TGraph(Map.at("X").size(),	
									Map.at("X").data(),			Map.at("R").data());

		if (Case == "Sigma")
		{
			// gr_SigmaVsW->SetFillColor(colors[0]);
			// gr_SigmaVsW->SetFillStyle(3005);
			gr_SigmaVsW->SetLineColor(colors[0]);
			gr_SigmaVsW->SetLineStyle(styles[0]);
			gr_SigmaVsW->SetLineWidth(2);
			gr_SigmaVsW->Draw("lsame");
			leg->AddEntry(gr_SigmaVsW,	Name[0].Data(),		"l" );

		}
		else if (Case == "R")
		{
			// gr_RvsW->SetFillColor(colors[0]);
			// gr_RvsW->SetFillStyle(3005);
			gr_RvsW->SetLineColor(colors[0]);
			gr_RvsW->SetLineStyle(styles[0]);
			gr_RvsW->SetLineWidth(2);
			gr_RvsW->Draw("lsame");
			leg->AddEntry(gr_RvsW,	Name[0].Data(),		"l" );
		}
	}
};

struct DSigmaDy_Sigma_R_bBK : Theory
{
	TGraph*	gr_AnAn		{nullptr};
	TGraph*	gr_Simgas	{nullptr};
	TGraph*	gr_R		{nullptr};

	DSigmaDy_Sigma_R_bBK(TString Case_ ) : Theory(Case_) {}
	virtual ~DSigmaDy_Sigma_R_bBK()	override
	{
		delete gr_AnAn;
		delete gr_Simgas;
		delete gr_R;
	}

	virtual void Draw(TLegend * leg) override
	{
		std::vector<TString> Name	 = {"bBK_GG",	"bBK_A"};

		std::vector<int> colors = {4, 4};
		std::vector<int> styles = {1, 7};

		MultiDataReader MDR1({"./inFiles/bBK_GG_Jpsi_Sigma.txt", "./inFiles/bBK_A_Jpsi_Sigma.txt"},
							{"X", "W", "Sigma"});
		MultiDataReader MDR2({"./inFiles/bBK_GG_Jpsi_Xsec.txt", "./inFiles/bBK_A_Jpsi_Xsec.txt"},
							{"Dy", "DSigmaDy_AnAn"});

		for (int i = 0; i < MDR1.SDRs.size(); ++i)
		{
			std::map<TString, std::vector<double>> Map ={
				{"W",			{MDR1.GetVec(i, "W")}},
				{"Sigma",		{MDR1.GetVec(i, "Sigma")}},
				{"Dy",			{MDR2.GetVec(i, "Dy")}},
				{"DSigmaDy_AnAn",{MDR2.GetVec(i, "DSigmaDy_AnAn")}},
				{"X",			{MDR1.GetVec(i, "X")}},
				{"R",			{}}
			};

			auto IA_pair = ImpulseApprox::GetImpulseApprox(Map.at("W"));
			for (int j = 0; j < Map.at("W").size(); ++j)
			{
				Map.at("Sigma")[j] *= 0.001;
				Map.at("R")		.push_back(sqrt(Map.at("Sigma")[j]/IA_pair.first[j]));
			}

			gr_AnAn 	= new TGraph(Map.at("Dy").size(),	Map.at("Dy").data(),	Map.at("DSigmaDy_AnAn").data());
			gr_Simgas 	= new TGraph(Map.at("W").size(),	Map.at("W").data(),	Map.at("Sigma").data());
			gr_R 		= new TGraph(Map.at("X").size(),	Map.at("X").data(),	Map.at("R").data());
		
			if (Case == "DSigmaDy_AnAn"){
				gr_AnAn->SetLineColor(colors[i]);
				gr_AnAn->SetLineStyle(styles[i]);
				gr_AnAn->SetLineWidth(2);
				gr_AnAn->Draw("lsame");
				leg->AddEntry(gr_AnAn,	Name[i].Data(),		"l" );
			}
			else if (Case == "Sigma"){
				gr_Simgas->SetLineColor(colors[i]);
				gr_Simgas->SetLineStyle(styles[i]);
				gr_Simgas->SetLineWidth(2);
				gr_Simgas->Draw("lsame");
				leg->AddEntry(gr_Simgas,	Name[i].Data(),		"l" );
			}
			else if (Case == "R"){
				gr_R->SetLineColor(colors[i]);
				gr_R->SetLineStyle(styles[i]);
				gr_R->SetLineWidth(2);
				gr_R->Draw("lsame");
				leg->AddEntry(gr_R,	Name[i].Data(),		"l" );
			}
		}
	}
};

struct DSigmaDy_Sigma_R_CD : Theory
{
	TGraph* gr_AnAn		{nullptr};
	TGraph* gr_Simgas	{nullptr};
	TGraph* gr_R		{nullptr};

	DSigmaDy_Sigma_R_CD(TString Case_ ) : Theory(Case_) {}
	virtual ~DSigmaDy_Sigma_R_CD()	override
	{
		delete gr_AnAn;
		delete gr_Simgas;
		delete gr_R;
	}

	virtual void Draw(TLegend * leg)	override
	{
		std::vector<TString> Name	 = {"CD_BGK",	"CD_GBW",	"CD_IIM"};
		std::vector<int> colors = {8, 8, 8};
		std::vector<int> styles = {1, 7, 3};

		MultiDataReader MDR1({"./inFiles/CD_BGK_Jpsi_Sigma.dat", "./inFiles/CD_GBW_Jpsi_Sigma.dat", "./inFiles/CD_IIM_Jpsi_Sigma.dat"},
							{"W", "Sigma"});
		MultiDataReader MDR2({"./inFiles/CD_BGK_Jpsi_Xsec.dat", "./inFiles/CD_GBW_Jpsi_Xsec.dat", "./inFiles/CD_IIM_Jpsi_Xsec.dat"},
							{"Dy", "DSigmaDy"});

		for (int i = 0; i < MDR1.SDRs.size(); ++i)
		{
			std::map<TString, std::vector<double>> Map ={
				{"W",			{MDR1.GetVec(i, "W")}},
				{"Sigma",		{MDR1.GetVec(i, "Sigma")}},
				{"Dy",		{MDR2.GetVec(i, "Dy")}},
				{"DSigmaDy",		{MDR2.GetVec(i, "DSigmaDy")}},
				{"X",			{}},
				{"R",			{}}
			};
			auto IA_pair = ImpulseApprox::GetImpulseApprox(Map.at("W"));
			for (int j = 0; j < Map.at("W").size(); ++j)
			{
				Map.at("Sigma")[j] *= 0.001;
				Map.at("X")	.push_back(ParamConverter::W2x( Map.at("W")[j] ) );
				Map.at("R")		.push_back(sqrt(Map.at("Sigma")[j]/IA_pair.first[j]));
			}
			for (int j = 0; j < Map.at("Dy").size(); ++j)
			{
				Map.at("DSigmaDy")[j] *= 0.001;
			}

			auto gr_AnAn	= new TGraph(Map.at("Dy").size(),	Map.at("Dy").data(),	Map.at("DSigmaDy").data());
			auto gr_Simgas 	= new TGraph(Map.at("W").size(),	Map.at("W").data(),	Map.at("Sigma").data());
			auto gr_R 		= new TGraph(Map.at("X").size(),	Map.at("X").data(),	Map.at("R").data());
		
			if (Case == "AnAn"){
				gr_AnAn->SetLineColor(colors[i]);
				gr_AnAn->SetLineStyle(styles[i]);
				gr_AnAn->SetLineWidth(2);
				gr_AnAn->Draw("lsame");
				leg->AddEntry(gr_AnAn,	Name[i].Data(),		"l" );
			}
			else if (Case == "Sigma"){
				gr_Simgas->SetLineColor(colors[i]);
				gr_Simgas->SetLineStyle(styles[i]);
				gr_Simgas->SetLineWidth(2);
				gr_Simgas->Draw("lsame");
				leg->AddEntry(gr_Simgas,	Name[i].Data(),		"l" );
			}
			else if (Case == "R"){
				gr_R->SetLineColor(colors[i]);
				gr_R->SetLineStyle(styles[i]);
				gr_R->SetLineWidth(2);
				gr_R->Draw("lsame");
				leg->AddEntry(gr_R,	Name[i].Data(),		"l" );
			}
		}
	}
};
