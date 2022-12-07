#ifndef ConversionAnalyzer_H
#define ConversionAnalyzer_H

#include "Analyzer.h"
#include "ParamConverter.C"

struct ConversionAnalyzer : Analyzer
{
    ConversionAnalyzer(AnalysisData& data_) : Analyzer(data_) {};

    void AddRap()
    {
        std::vector<double> temp;
        int n_Dy = data.Get("Dy").size();
        for (int i = 0; i < n_Dy; i++)
        {
            temp.push_back(data.Get("Dy", i)); temp.push_back( -data.Get("Dy", i) );
        }
        data.Add("Rap", temp);
    }

    void Conversions()
    {
        std::vector<double> W_temp;
        std::vector<double> X_temp;
        std::vector<double> rap_temp = data.Get("Rap");
        for (int i = 0; i < rap_temp.size(); i++)
        {
            W_temp.push_back( ParamConverter::y2W( rap_temp[i]) );
            X_temp.push_back( ParamConverter::y2x( rap_temp[i]) );
        }

        data.Add("W", W_temp);
        data.Add("X", X_temp);
    }

    void Handle() override
    {
        cout << endl << "+++ConversionAnalyzer::Handling..." <<endl;
        
        AddRap();
        Conversions();

        Analyzer::Handle();
    }
};



#endif