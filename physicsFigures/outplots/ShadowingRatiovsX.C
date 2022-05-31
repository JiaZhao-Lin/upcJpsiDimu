#ifdef __CLING__
#pragma cling optimize(0)
#endif
void ShadowingRatiovsX()
{
//=========Macro generated from canvas: c1/c1
//=========  (Thu May 26 11:17:56 2022) by ROOT version 6.26/02
   Double_t inf = 10e6;
   TCanvas *c1 = new TCanvas("c1", "c1",10,53,800,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-4.750615,0.06987952,-1.81166,1.154217);
   c1->SetFillColor(10);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogx();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.12);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.12);
   c1->SetFrameFillStyle(4000);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderSize(2);
   c1->SetFrameFillStyle(4000);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderSize(2);
   
   TH2D *htem2d__1 = new TH2D("htem2d__1","",10,4e-05,0.011,10,0.2,1.1);
   htem2d__1->SetLineWidth(2);
   htem2d__1->SetMarkerStyle(20);
   htem2d__1->GetXaxis()->SetTitle("x");
   htem2d__1->GetXaxis()->CenterTitle(true);
   htem2d__1->GetXaxis()->SetLabelFont(42);
   htem2d__1->GetXaxis()->SetLabelOffset(0.008);
   htem2d__1->GetXaxis()->SetLabelSize(0.04);
   htem2d__1->GetXaxis()->SetTitleSize(0.07);
   htem2d__1->GetXaxis()->SetTickLength(0.04);
   htem2d__1->GetXaxis()->SetTitleOffset(0.69);
   htem2d__1->GetXaxis()->SetTitleFont(42);
   htem2d__1->GetYaxis()->SetTitle("R^{Pb}_{g}(x, #mu^{2}=2.4 GeV^{2})");
   htem2d__1->GetYaxis()->SetLabelFont(42);
   htem2d__1->GetYaxis()->SetLabelOffset(0.008);
   htem2d__1->GetYaxis()->SetLabelSize(0.04);
   htem2d__1->GetYaxis()->SetTitleSize(0.06);
   htem2d__1->GetYaxis()->SetTickLength(0.02);
   htem2d__1->GetYaxis()->SetTitleOffset(0.85);
   htem2d__1->GetYaxis()->SetTitleFont(42);
   htem2d__1->GetZaxis()->SetLabelFont(42);
   htem2d__1->GetZaxis()->SetLabelOffset(0.008);
   htem2d__1->GetZaxis()->SetLabelSize(0.05);
   htem2d__1->GetZaxis()->SetTitleSize(0.06);
   htem2d__1->GetZaxis()->SetTickLength(0.02);
   htem2d__1->GetZaxis()->SetTitleOffset(0.95);
   htem2d__1->GetZaxis()->SetTitleFont(42);
   htem2d__1->Draw("");
   
   Double_t Graph0_fx3001[6] = {
   0.0001072038,
   0.003550104,
   8.349044e-05,
   0.004558424,
   6.502242e-05,
   0.005853132};
   Double_t Graph0_fy3001[6] = {
   inf,
   inf,
   inf,
   inf,
   inf,
   inf};
   Double_t Graph0_felx3001[6] = {
   5e-06,
   0.00015,
   4e-06,
   0.0002,
   3e-06,
   0.00025};
   Double_t Graph0_fely3001[6] = {
   0.03043005,
   0.04730371,
   0.02776808,
   0.04352832,
   0.02522762,
   0.04818824};
   Double_t Graph0_fehx3001[6] = {
   5e-06,
   0.00015,
   4e-06,
   0.0002,
   3e-06,
   0.00025};
   Double_t Graph0_fehy3001[6] = {
   0.03043005,
   0.04730371,
   0.02776808,
   0.04352832,
   0.02522762,
   0.04818824};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6,Graph0_fx3001,Graph0_fy3001,Graph0_felx3001,Graph0_fehx3001,Graph0_fely3001,Graph0_fehy3001);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = 1182;
   color = new TColor(ci, 0.7, 0.7, 0.7, " ", 0.7);
   grae->SetFillColor(ci);
   grae->SetMarkerStyle(24);
   
   TH1F *Graph_Graph03001 = new TH1F("Graph_Graph03001","Graph",100,5.582018e-05,0.006707243);
   Graph_Graph03001->SetMinimum(100);
   Graph_Graph03001->SetMaximum(100);
   Graph_Graph03001->SetDirectory(0);
   Graph_Graph03001->SetStats(0);
   Graph_Graph03001->SetLineWidth(2);
   Graph_Graph03001->SetMarkerStyle(20);
   Graph_Graph03001->GetXaxis()->SetLabelFont(42);
   Graph_Graph03001->GetXaxis()->SetLabelOffset(0.008);
   Graph_Graph03001->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph03001->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph03001->GetXaxis()->SetTickLength(0.02);
   Graph_Graph03001->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph03001->GetXaxis()->SetTitleFont(42);
   Graph_Graph03001->GetYaxis()->SetLabelFont(42);
   Graph_Graph03001->GetYaxis()->SetLabelOffset(0.008);
   Graph_Graph03001->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph03001->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph03001->GetYaxis()->SetTickLength(0.02);
   Graph_Graph03001->GetYaxis()->SetTitleOffset(0.95);
   Graph_Graph03001->GetYaxis()->SetTitleFont(42);
   Graph_Graph03001->GetZaxis()->SetLabelFont(42);
   Graph_Graph03001->GetZaxis()->SetLabelOffset(0.008);
   Graph_Graph03001->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph03001->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph03001->GetZaxis()->SetTickLength(0.02);
   Graph_Graph03001->GetZaxis()->SetTitleOffset(0.95);
   Graph_Graph03001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph03001);
   
   grae->Draw("2");
   
   Double_t Graph1_fx1001[6] = {
   0.0001072038,
   0.003550104,
   8.349044e-05,
   0.004558424,
   6.502242e-05,
   0.005853132};
   Double_t Graph1_fy1001[6] = {
   inf,
   inf,
   inf,
   inf,
   inf,
   inf};
   Double_t Graph1_fex1001[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph1_fey1001[6] = {
   100,
   100,
   100,
   100,
   100,
   100};
   TGraphErrors *gre = new TGraphErrors(6,Graph1_fx1001,Graph1_fy1001,Graph1_fex1001,Graph1_fey1001);
   gre->SetName("Graph1");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph11001 = new TH1F("Graph_Graph11001","Graph",100,5.852018e-05,0.006431943);
   Graph_Graph11001->SetMinimum(100);
   Graph_Graph11001->SetMaximum(1000);
   Graph_Graph11001->SetDirectory(0);
   Graph_Graph11001->SetStats(0);
   Graph_Graph11001->SetLineWidth(2);
   Graph_Graph11001->SetMarkerStyle(20);
   Graph_Graph11001->GetXaxis()->SetLabelFont(42);
   Graph_Graph11001->GetXaxis()->SetLabelOffset(0.008);
   Graph_Graph11001->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph11001->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph11001->GetXaxis()->SetTickLength(0.02);
   Graph_Graph11001->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph11001->GetXaxis()->SetTitleFont(42);
   Graph_Graph11001->GetYaxis()->SetLabelFont(42);
   Graph_Graph11001->GetYaxis()->SetLabelOffset(0.008);
   Graph_Graph11001->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph11001->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph11001->GetYaxis()->SetTickLength(0.02);
   Graph_Graph11001->GetYaxis()->SetTitleOffset(0.95);
   Graph_Graph11001->GetYaxis()->SetTitleFont(42);
   Graph_Graph11001->GetZaxis()->SetLabelFont(42);
   Graph_Graph11001->GetZaxis()->SetLabelOffset(0.008);
   Graph_Graph11001->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph11001->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph11001->GetZaxis()->SetTickLength(0.02);
   Graph_Graph11001->GetZaxis()->SetTitleOffset(0.95);
   Graph_Graph11001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph11001);
   
   gre->Draw("pez");
   
   Double_t Graph2_fx1002[1] = {
   0.001};
   Double_t Graph2_fy1002[1] = {
   0.65};
   Double_t Graph2_fex1002[1] = {
   0};
   Double_t Graph2_fey1002[1] = {
   0.03};
   gre = new TGraphErrors(1,Graph2_fx1002,Graph2_fy1002,Graph2_fex1002,Graph2_fey1002);
   gre->SetName("Graph2");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);
   gre->SetLineColor(4);
   gre->SetLineWidth(2);
   gre->SetMarkerColor(4);
   gre->SetMarkerStyle(24);
   
   TH1F *Graph_Graph21002 = new TH1F("Graph_Graph21002","Graph",100,0.0009,1.101);
   Graph_Graph21002->SetMinimum(0.614);
   Graph_Graph21002->SetMaximum(0.686);
   Graph_Graph21002->SetDirectory(0);
   Graph_Graph21002->SetStats(0);
   Graph_Graph21002->SetLineWidth(2);
   Graph_Graph21002->SetMarkerStyle(20);
   Graph_Graph21002->GetXaxis()->SetLabelFont(42);
   Graph_Graph21002->GetXaxis()->SetLabelOffset(0.008);
   Graph_Graph21002->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph21002->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph21002->GetXaxis()->SetTickLength(0.02);
   Graph_Graph21002->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph21002->GetXaxis()->SetTitleFont(42);
   Graph_Graph21002->GetYaxis()->SetLabelFont(42);
   Graph_Graph21002->GetYaxis()->SetLabelOffset(0.008);
   Graph_Graph21002->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph21002->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph21002->GetYaxis()->SetTickLength(0.02);
   Graph_Graph21002->GetYaxis()->SetTitleOffset(0.95);
   Graph_Graph21002->GetYaxis()->SetTitleFont(42);
   Graph_Graph21002->GetZaxis()->SetLabelFont(42);
   Graph_Graph21002->GetZaxis()->SetLabelOffset(0.008);
   Graph_Graph21002->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph21002->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph21002->GetZaxis()->SetTickLength(0.02);
   Graph_Graph21002->GetZaxis()->SetTitleOffset(0.95);
   Graph_Graph21002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph21002);
   
   gre->Draw("pez");
   
   TLegend *leg = new TLegend(0.17,0.6,0.55,0.8,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph1","CMS","lp");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","ALICE","lp");
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   TLatex *   tex = new TLatex(0.15,0.86,"Pb+Pb UPC #sqrt{s_{NN}} = 5.02 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.42,0.22,"Nuclear suppression factor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
