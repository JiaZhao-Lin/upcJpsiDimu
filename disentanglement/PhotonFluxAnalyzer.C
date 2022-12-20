#ifndef PhotonFluxAnalyzer_H
#define PhotonFluxAnalyzer_H

#include "Analyzer.h"
#include "ParamConverter.C"

struct PhotonFluxAnalyzer : Analyzer
{
protected:
    std::map<TString, std::vector<double>> PhotonFluxMap	=   //map template for loading flux
    {
        {"Energy_Table_AnAn"    ,{}}    ,{"Rap_Table_AnAn"      ,{}}    ,{"dNdk_Table_AnAn"     ,{}}    ,{"dNdy_Table_AnAn"     ,{}},
        {"Energy_Table_0n0n"    ,{}}    ,{"Rap_Table_0n0n"      ,{}}    ,{"dNdk_Table_0n0n"     ,{}}    ,{"dNdy_Table_0n0n"     ,{}},
        {"Energy_Table_0nXnSum" ,{}}    ,{"Rap_Table_0nXnSum"   ,{}}    ,{"dNdk_Table_0nXnSum"  ,{}}    ,{"dNdy_Table_0nXnSum"  ,{}},
        {"Energy_Table_XnXn"    ,{}}    ,{"Rap_Table_XnXn"      ,{}}    ,{"dNdk_Table_XnXn"     ,{}}    ,{"dNdy_Table_XnXn"     ,{}},
        {"biter_0n0n"           ,{}}    ,{"PofPhotonB_0n0n"     ,{}}    ,{"PofHadronB_0n0n"     ,{}}    ,{"PofB_0n0n"           ,{}},
        {"biter_0nXnSum"        ,{}}    ,{"PofPhotonB_0nXnSum"  ,{}}    ,{"PofHadronB_0nXnSum"  ,{}}    ,{"PofB_0nXnSum"        ,{}},
        {"biter_XnXn"           ,{}}    ,{"PofPhotonB_XnXn"     ,{}}    ,{"PofHadronB_XnXn"     ,{}}    ,{"PofB_XnXn"           ,{}}
    };
    std::map<TString, std::vector<double>> PhotonFluxMapClean = PhotonFluxMap;   //remember this template and use it to reset the map
    std::map<TString, std::vector<double>> InterPhotonFluxMap;   //use this map to remember the interpolated flux based on the rap
    bool getFluxErr = false;

public:
    TString inFileDir;
    TString subCase;
    std::vector<TString> CasesName = {"AnAn",	"0n0n", "0nXnSum", "XnXn"};

    PhotonFluxAnalyzer(AnalysisData& data_, TString inFileDir_, TString subCase_) : Analyzer{data_}, inFileDir(inFileDir_), subCase(subCase_) {};

    void LoadPhotonFlux(TString inFileDir_, TString subCase_)
    {
        for (int i = 0; i < CasesName.size(); ++i)
        {
            TString FluxFileName  = Form("Flux_%s%s.txt", CasesName[i].Data(), subCase_.Data());
            TString PofBFileName  = Form("PofB_%s%s.txt", CasesName[i].Data(), subCase_.Data());
            cout<<"LoadPhotonFlux-------->Loading Photon Flux From " + FluxFileName + " in the Dir:"<<inFileDir_<<endl;

            ifstream myfile(Form("%s%s", inFileDir_.Data(), FluxFileName.Data()));

            if (myfile.is_open())
            {
                std::string line;
                // Read one line at a time into the variable line:
                while(std::getline(myfile, line))
                {
                    std::vector<double>   	lineData;
                    std::stringstream  		lineStream(line);

                    double value;
                    // Read an integer at a time from the line
                    while(lineStream >> value)
                    {
                        // Add the integers from a line to a 1D array (vector)
                        lineData.push_back(value);
                    }
                    // cout<<line<<endl;
                    // When all the integers have been read, add the 1D array
                    PhotonFluxMap.at("Energy_Table_" + CasesName[i]).push_back(lineData[0]);
                    PhotonFluxMap.at("dNdk_Table_" + CasesName[i]) 	.push_back(lineData[1]);
                    PhotonFluxMap.at("dNdy_Table_" + CasesName[i])	.push_back(lineData[2]);

                    PhotonFluxMap.at("Rap_Table_" + CasesName[i])	.push_back( ParamConverter::w2y( lineData[0] )	);
                }
            }
            else throw std::runtime_error( "ERROR!!! Unable to open Flux file!!!");

            if (i == 0) continue;
            ifstream myfile1(Form("%s%s", inFileDir_.Data(), PofBFileName.Data()));
            if (myfile1.is_open())
            {
                std::string line;
                // Read one line at a time into the variable line:
                while(std::getline(myfile1, line))
                {
                    std::vector<double>   	lineData;
                    std::stringstream  		lineStream(line);

                    double value;
                    // Read an integer at a time from the line
                    while(lineStream >> value)
                    {
                        // Add the integers from a line to a 1D array (vector)
                        lineData.push_back(value);
                    }
                    // cout<<line<<endl;
                    // When all the integers have been read, add the 1D array
                    PhotonFluxMap.at("biter_" + CasesName[i])		.push_back(lineData[0]);
                    PhotonFluxMap.at("PofPhotonB_" + CasesName[i]) 	.push_back(lineData[1]);
                    PhotonFluxMap.at("PofHadronB_" + CasesName[i]) 	.push_back(lineData[2]);
                    PhotonFluxMap.at("PofB_" + CasesName[i])  		.push_back(lineData[3]);
                }
            }
            else cout << "ERROR!!! Unable to open PofB file!!!";
        }
        cout<<"LoadPhotonFlux-------->DONE "<<endl<<endl;
    }

    double Interpolate(const double Egamma, TString Case)
    {
        const int nstep = 100;
        const double Emin = 1.000000e-05, Emax = 4.777493e+02; //Emax = 2.627452e+02 for 2.76TeV

        double lnEmin=log(Emin);
        double lnEmax=log(Emax);
        double lEgamma=log(Egamma);
        double dlnE=(lnEmax-lnEmin)/nstep; 

        //       >> Egamma between Ilt and Ilt+1
        int Ilt = int((lEgamma-lnEmin)/dlnE);
        //       >> ln(Egamma) for first point 
        double lnElt = lnEmin + Ilt*dlnE; 
        //       >> Interpolate
        double flux_r = PhotonFluxMap.at("dNdy_Table_" + Case)[Ilt-1] + ((lEgamma-lnElt)/dlnE)*(PhotonFluxMap.at("dNdy_Table_" + Case)[Ilt]- PhotonFluxMap.at("dNdy_Table_" + Case)[Ilt-1]);
        flux_r = flux_r/Egamma;

        return flux_r;
    }


    // calculate the photon flux from the rap given in the data.
    // Will add the flux result to the data if it is not there already, otherwise it will do NOTHING
    // This way the function can be called multiple times without any problem and be used to calculate the flux for uncertainty
    void InterpolateFlux(TString inFileDir_, TString subCase_)
    {
        // if (PhotonFluxMap.at("Energy_Table_AnAn").size() == 0){	loadPhotonFlux(inFileDir, subCase);	}

        //Reset the photon flux map to the clean one
        PhotonFluxMap.clear();
        InterPhotonFluxMap.clear();
        PhotonFluxMap = PhotonFluxMapClean;
        LoadPhotonFlux( inFileDir_, subCase_);

        //Calculate photon energy from the rap
        int n_data = data.GetSize("Rap");
        std::vector<double> w_temp;
        for (int j = 0; j < n_data; ++j)
        {
            double w = ParamConverter::y2w( data.Get("Rap", j) );
            w_temp      .push_back( w );
        }

        //Add the photon energy if it is not in the map already
        if ( !data.IsMapKeyExist("w") ) data.Add("w", w_temp);

        //Calculate photon flux from the rap
        for (int i = 0; i < CasesName.size(); ++i)
        {
            TString Case = CasesName[i];
            std::vector<double> dNdy_temp;

            for (int j = 0; j < n_data; ++j)
            {
                double w = data.Get("w", j);
                dNdy_temp   .push_back(	w * Interpolate(w, Case)	);
            }

            //Remember this interpolated flux
            InterPhotonFluxMap["dNdy_" + Case] = dNdy_temp;

            //Add the photon flux if it is not in the map already
            if ( !data.IsMapKeyExist("dNdy_" + Case) )  data.Add("dNdy_" + Case, dNdy_temp);
        }
        
        cout<<"InterpolateFlux-------->DONE"<<endl<<endl;
    }

    //Calculate the error of photon flux based on the file given by the simulation
    void CalculateFluxErr()
    {
        std::vector<TString> inFileDirList  = {"../simulation/flux/", "../simulation/flux/", "../simulation/flux/", "../simulation/flux/",};
        std::vector<TString> subCaseList    = {"_SigNN68p3R6p70a0p59", "_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p53", "_SigNN68p3R6p64a0p59"};
        std::vector<std::vector<double>> dNdy_Errs  (CasesName.size(),    std::vector<double> (data.GetSize("Rap"),   0)  );
        std::vector<std::vector<double>> dNdy_Uncers(CasesName.size(),    std::vector<double> (data.GetSize("Rap"),   0)  );

        assert(inFileDirList.size() == subCaseList.size());
        for (int i = 0; i < inFileDirList.size(); ++i)
        {
            TString inFileDir_  = inFileDirList[i];
            TString subCase_    = subCaseList[i];

            InterpolateFlux(inFileDir_,  subCase_);

            for (int j = 0; j < CasesName.size(); ++j)
            {
                TString Case = CasesName[j];

                for (int k = 0; k < data.GetSize("Rap"); ++k)
                {
                    auto temp_flux = InterPhotonFluxMap.at("dNdy_" + Case)[k];
                    auto temp_flux_default = data.Get("dNdy_" + Case, k);
                    auto temp_err = abs(temp_flux - temp_flux_default);
                    dNdy_Errs[j][k] = max(dNdy_Errs[j][k],  temp_err);
                    dNdy_Uncers[j][k] = dNdy_Errs[j][k] / temp_flux_default * 100.0;
                }
            }
        }

        for (int i = 0; i < CasesName.size(); ++i)
        {
            data.Add("dNdy_" + CasesName[i] + "_Err", dNdy_Errs[i]);
            data.Add("dNdy_" + CasesName[i] + "_Uncer", dNdy_Uncers[i]);
        }
    }

    void SetFluxErr(const bool getFluxErr_ = false)
    {
        getFluxErr = getFluxErr_;
    }

    void Handle() override
    {
        cout << endl << "+++PhotonFluxAnalyzer::Handling..." <<endl;
        InterpolateFlux(inFileDir, subCase);

        if (getFluxErr) CalculateFluxErr();

        // Analyzer::Handle();
    }
};


// void PhotonFluxAnalyzer()
// {
// 	//------------------Standard Test-------------------------------------------
//     // AnalysisData data("Test");
//     // AnalysisDataObserver obs;
//     // data.Subscribe(&obs);
//     // data.Add("Rap",         {1.75, -1.75, 2, -2, 2.25, -2.25});
// 	// struct PhotonFluxAnalyzer ana(data, "../simulation/flux/", "_SigNN68p3R6p67a0p56");
//     // ana.Handle();

// 	// plotFlux(TestMap,	"0n0n");
// 	// plotFlux(TestMap,	"0nXnSum");
// 	// plotFlux(TestMap,	"XnXn");
// 	// plotFlux(TestMap,	"AnAn");
// 	// plotPofB();
// 	//--------------------------------------------------------------------------

// 	//--------------------------------------------------------------------------
// 	//------------------Flux For Drawing ALICE and LHCb values------------------
// 	// std::map<TString, std::vector<double>>	TestMap = 
// 	// {
// 	// 	{"Raps",	{ 0, -3.875, -3.625, -4.25, -3.75 }},
// 	// 	{"dNdy",	{}}
// 	// };

// 	// auto TestMap1 = TestMap;
// 	// auto TestMap2 = TestMap;
// 	// InterpolateFlux(TestMap, "flux/", "_SigNN68p3R6p67a0p56");
// 	// plotFlux(TestMap,	"0n0n");
// 	// plotFlux(TestMap,	"0nXnSum");
// 	// plotFlux(TestMap,	"XnXn");
// 	// plotFlux(TestMap,	"AnAn");

// 	// InterpolateFlux(TestMap1, "flux/", "_SigNN68p3R6p64a0p53");
// 	// InterpolateFlux(TestMap2, "flux/", "_SigNN68p3R6p70a0p59");
// 	// auto v1 = getFluxUncer(TestMap.at("dNdy_AnAn"),	TestMap1.at("dNdy_AnAn"));
// 	// auto v2 = getFluxUncer(TestMap.at("dNdy_AnAn"),	TestMap2.at("dNdy_AnAn"));
// 	// for (int i = 0; i < v1.size(); ++i)
// 	// {
// 	// 	auto uncer = (v1[i] > v2[i]) ? v1[i]: v2[i];
// 	// 	cout << "Raps:" << TestMap.at("Raps")[i] << "	Flux:" << TestMap.at("dNdy_AnAn")[i]; cout<< "	Uncer:" << uncer <<endl;
// 	// }

// 	// CompareFlux("AnAn",		{"_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59", "_SigNN68p3R6p64a0p59", "_SigNN68p3R6p70a0p53"},	{"#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59", "#sigma_{NN} = 68.3, R = 6.64, a = 0.59", "#sigma_{NN} = 68.3, R = 6.70, a = 0.53"});
// 	// CompareFlux("0n0n",		{"_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59", "_SigNN68p3R6p64a0p59", "_SigNN68p3R6p70a0p53"},	{"#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59", "#sigma_{NN} = 68.3, R = 6.64, a = 0.59", "#sigma_{NN} = 68.3, R = 6.70, a = 0.53"});
// 	// CompareFlux("0nXnSum",	{"_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59", "_SigNN68p3R6p64a0p59", "_SigNN68p3R6p70a0p53"},	{"#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59", "#sigma_{NN} = 68.3, R = 6.64, a = 0.59", "#sigma_{NN} = 68.3, R = 6.70, a = 0.53"});
// 	// CompareFlux("XnXn",		{"_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59", "_SigNN68p3R6p64a0p59", "_SigNN68p3R6p70a0p53"},	{"#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59", "#sigma_{NN} = 68.3, R = 6.64, a = 0.59", "#sigma_{NN} = 68.3, R = 6.70, a = 0.53"});

// 	// CompareFlux("AnAn",		{"_SigNN68p3R6p67a0p56", "_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59"},	{"#sigma_{NN} = 68.3, R = 6.67, a = 0.56", "#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59"});
// 	// CompareFlux("0n0n",		{"_SigNN68p3R6p67a0p56", "_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59"},	{"#sigma_{NN} = 68.3, R = 6.67, a = 0.56", "#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59"});
// 	// CompareFlux("0nXnSum",	{"_SigNN68p3R6p67a0p56", "_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59"},	{"#sigma_{NN} = 68.3, R = 6.67, a = 0.56", "#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59"});
// 	// CompareFlux("XnXn",		{"_SigNN68p3R6p67a0p56", "_SigNN68p3R6p64a0p53", "_SigNN68p3R6p70a0p59"},	{"#sigma_{NN} = 68.3, R = 6.67, a = 0.56", "#sigma_{NN} = 68.3, R = 6.64, a = 0.53", "#sigma_{NN} = 68.3, R = 6.70, a = 0.59"});
// 	//--------------------------------------------------------------------------
// }

#endif