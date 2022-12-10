#ifndef SHADOWINGRATIOANALYZER_H
#define SHADOWINGRATIOANALYZER_H

#include "Analyzer.h"
#include "ImpulseApprox.C"

struct ShadowingRatioAnalyzer : Analyzer
{
    bool Simple_R_SysErr = false;

    ShadowingRatioAnalyzer(AnalysisData& data_) : Analyzer{data_} {};

    void Set_Simple_R_SysErr(bool val) { Simple_R_SysErr = val; }

    void PrintHandling() const
    {
        cout << endl << std::setfill('+') << std::setw(100) << Form("%s::Handle Is Handling The Business... ", typeid(this).name()) << endl;
        cout << setfill(' ');
    }

    void Calculate_R()
    {
        int n_data = data.GetSize("Sigma");
        std::vector<double> temp;
        for (int i = 0; i < n_data; ++i)
        {
            temp.push_back( sqrt( data.Get("Sigma", i)/data.Get("Sigma_IA", i) ) );
        }
        data.Add("R", temp);
    }

    void Calculate_R_Err()
    {
        int n_data = data.GetSize("Sigma");
        std::vector<double> temp;
        for (int i = 0; i < n_data; ++i)
        {
            temp.push_back( 0.5 * data.Get("R", i) * data.Get("Sigma_Err", i)/data.Get("Sigma", i) );
        }
        data.Add("R_Err", temp);
    }

    void Calculate_Simple_R_SysErr()
    {
        int n_data = data.GetSize("Sigma");
        std::vector<double> temp;
        for (int i = 0; i < n_data; ++i)
        {
            temp.push_back( 0.5 * data.Get("R", i) * TMath::Hypot( data.Get("Sigma_SysErr", i)/data.Get("Sigma", i),
																data.Get("Sigma_IA_Err", i)/data.Get("Sigma_IA", i) ) );
        }
        data.Add("R_SysErr", temp);
    }
    
    void Handle() override
    {
        cout << endl << "+++ShadowingRatioAnalyzer::Handling..." <<endl;

        //Calculate Sigma_IA
        auto IA_pair = ImpulseApprox::GetImpulseApprox(data.Get("W"));
        data.Add("Sigma_IA",    IA_pair.first);
        data.Add("Sigma_IA_Err",IA_pair.second);

        //Calculate R must be first!
        Calculate_R();
        Calculate_R_Err();
        if (Simple_R_SysErr) Calculate_Simple_R_SysErr();

        //Moving to the next Analyzer
        Analyzer::Handle();
    }
};


// void ShadowingRatioAnalyzer()
// {
//     AnalysisData data("Test");
//     AnalysisDataObserver obs;
//     data.Subscribe(&obs);
//     data.Add("Sigma",       {1,2,3});
//     data.Add("Rap",         {1,2,3});
//     data.Add("W",           {10,22,33});
//     data.Add("Sigma_Err",   {0.1,0.2,0.3});
//     struct ShadowingRatioAnalyzer R_Ana(data);
//     R_Ana.Handle();
// }

#endif