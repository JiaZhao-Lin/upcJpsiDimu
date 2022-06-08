#include "PhysParameters.h"
#include "../simulation/getPhotonFlux.C"


std::vector<std::vector<double>> fit2D(std::map<TString, std::vector<double>> Map_Xsec, const double flux_uncer = 1.0)
{
   std::vector<std::vector<double>> X_RapBin, Y_RapBin, Z_RapBin, Z_Err_RapBin;
   std::vector<double> Raps,   Raps_Err;
   std::vector<double> Sigmas, Sigmas_Err;

   //getting xsection from the map
   std::vector<double> Rap = Map_Xsec["Rap"], RapErr = Map_Xsec["RapErr"];
   std::vector<double> Xsec_0n0n = Map_Xsec["Xsec_0n0n"], XsecErr_0n0n = Map_Xsec["XsecErr_0n0n"];
   std::vector<double> Xsec_0nXnSum = Map_Xsec["Xsec_0nXnSum"], XsecErr_0nXnSum = Map_Xsec["XsecErr_0nXnSum"];
   std::vector<double> Xsec_XnXn = Map_Xsec["Xsec_XnXn"], XsecErr_XnXn = Map_Xsec["XsecErr_XnXn"];

   for (int i = 0; i < Rap.size(); ++i)
   {
      Raps        .push_back(Rap[i]);              Raps        .push_back( - Rap[i] );
      Raps_Err    .push_back(RapErr[i]);           Raps_Err    .push_back( RapErr[i] );
   }

   //Getting flux from getPhotonFlux.C
   std::map<TString, std::vector<double>> Temp_Map = {{"Raps",Raps}};
   InterpolateFlux(Temp_Map,"../simulation/flux/");
   std::vector<TVector2> flux_0n0n;
   std::vector<TVector2> flux_0nXnSum;
   std::vector<TVector2> flux_XnXn;
   for (int i = 0; i < Raps.size(); i+=2)
   {
      TVector2 temp_0n0n = {Temp_Map.at("dNdy_0n0n")[i], Temp_Map.at("dNdy_0n0n")[i+1]};
      TVector2 temp_0nXnSum = {Temp_Map.at("dNdy_0nXnSum")[i], Temp_Map.at("dNdy_0nXnSum")[i+1]};
      TVector2 temp_XnXn = {Temp_Map.at("dNdy_XnXn")[i], Temp_Map.at("dNdy_XnXn")[i+1]};

      flux_0n0n.push_back(temp_0n0n);
      flux_0nXnSum.push_back(temp_0nXnSum);
      flux_XnXn.push_back(temp_XnXn);
   }

   //Start filling points
   for (int i = 0; i < flux_0n0n.size(); ++i)
   {
      std::vector<double> X, Y, Z, Z_Err;
      X.push_back(flux_0n0n[i].X()*flux_uncer);      X.push_back(flux_0nXnSum[i].X()*flux_uncer);      X.push_back(flux_XnXn[i].X()*flux_uncer);
      Y.push_back(flux_0n0n[i].Y()*flux_uncer);      Y.push_back(flux_0nXnSum[i].Y()*flux_uncer);      Y.push_back(flux_XnXn[i].Y()*flux_uncer);
      Z.push_back(Xsec_0n0n[i]);          Z.push_back(Xsec_0nXnSum[i]);          Z.push_back(Xsec_XnXn[i]);
      Z_Err.push_back(XsecErr_0n0n[i]);   Z_Err.push_back(XsecErr_0nXnSum[i]);   Z_Err.push_back(XsecErr_XnXn[i]);
      X_RapBin.push_back(X);
      Y_RapBin.push_back(Y);
      Z_RapBin.push_back(Z);
      Z_Err_RapBin.push_back(Z_Err);
   }

   for (int i = 0; i < X_RapBin.size(); ++i)
   {
      cout<<"fit2D: Proceessing ------------------>Rap: "<<Rap[i]<<" <------------------------------- "<<endl;
      TGraph2DErrors* gr = new TGraph2DErrors(X_RapBin[i].size(), &X_RapBin[i][0], &Y_RapBin[i][0], &Z_RapBin[i][0],0,0, &Z_Err_RapBin[i][0]);
      TF2 * f = new TF2("func","[0] * x + [1] * y");
      f->SetParameters(0.001,0.5);

      auto c = new TCanvas();
      gr->Fit(f);
      gr->SetTitle("TGraph2D TF2 Fit; dN_{1}/dy; dN_{2}/dy; d#sigma/dy");
      gr->SetMarkerColor(kBlue);
      gr->SetMarkerSize(1.5);
      gr->SetLineWidth(5);
      gr->Draw("pez");
      // f->SetMarkerColor(kBlue);
      f->Draw("same surf");

      drawLatex(0.05, 0.95, Form("(#sigma(y = %.2f), #sigma(y = %.2f)) = (%.4f #pm %.4f, %.4f #pm %.4f)", Rap[i], -Rap[i],
                              f->GetParameter(0), f->GetParError(0), f->GetParameter(1), f->GetParError(1)),      42,       0.05,      1);
      
      c->SaveAs(Form("outplots/fit2D_%d.pdf",i));
      cout<<endl;
      
      Sigmas      .push_back(f->GetParameter(0));  Sigmas      .push_back(f->GetParameter(1));
      Sigmas_Err  .push_back(f->GetParError(0));   Sigmas_Err  .push_back(f->GetParError(1));

      delete gr; delete c; delete f;
   }
   return {Raps, Raps_Err, Sigmas, Sigmas_Err};
}