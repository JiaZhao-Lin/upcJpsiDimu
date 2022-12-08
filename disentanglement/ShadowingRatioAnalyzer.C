#ifndef SHADOWINGRATIOANALYZER_H
#define SHADOWINGRATIOANALYZER_H

#include "ImpulseApproxAnalyzer.C"

struct ShadowingRatioAnalyzer : ImpulseApproxAnalyzer
{
    ShadowingRatioAnalyzer(AnalysisData& data_) : ImpulseApproxAnalyzer(data_) {};


    void Calculate_R()
    {
        int n_data = data.GetSize("Rap");
        std::vector<double> temp;
        for (int i = 0; i < n_data; ++i)
        {
            temp.push_back( sqrt( data.Get("Sigma", i)/data.Get("Sigma_IA", i) ) );
        }
        data.Add("R", temp);
    }

    void Calculate_R_Err()
    {
        int n_data = data.GetSize("Rap");
        std::vector<double> temp;
        for (int i = 0; i < n_data; ++i)
        {
            temp.push_back( 0.5 * data.Get("R", i) * data.Get("Sigma_Err", i)/data.Get("Sigma", i) );
        }
        data.Add("R_Err", temp);
    }
    
    
    void Handle() override
    {
        cout << endl << "+++ShadowingRatioAnalyzer::Handling..." <<endl;
        //Calculate R must be first!
        GetImpulseApprox();
        Calculate_R();
        Calculate_R_Err();

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