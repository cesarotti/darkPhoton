#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <sstream> 

void noisePlots()

{
  int n = 14;
  double xtalNum[14] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 50};
  double reso1g[14] = {31.95, 23.65, 18.93, 16.55, 14.98, 13.05, 11.06, 9.93, 9.12, 7.67, 4.593, 2.33, 2.04, .494};
  double reso1gn1[14] = {29.26, 23.25, 19.84, 16.54, 14.55, 13.21, 11.80, 10.51, 9.39, 8.96, 7.055, 7.01, 6.72, 6.50};
  double reso1gn3[14] = {29.83, 23.53, 19.57, 17.80, 16.18, 14.72, 13.99, 13.97, 14.10, 13.84, 13.93, 13.64, 13.61, 15.58};

  TGraph *gr1g = new TGraph (n, xtalNum, reso1g);
  TGraph *gr1gn1 = new TGraph (n, xtalNum, reso1gn1);
  TGraph *gr1gn3 = new TGraph (n, xtalNum, reso1gn3);

  TCanvas *canvas = new TCanvas ("canvas", "Energy Resolution v. Crystal Num", 1000, 500);
  canvas->Divide(2, 2);
  canvas->cd(1);
  gr1gn3->SetMarkerStyle(21); gr1gn3->SetMarkerColor(55);
  gr1gn3->Draw("APL");
  gr1gn1->SetMarkerStyle(21); gr1gn1->SetMarkerColor(13); gr1gn1->Draw("PL");
  gr1g->SetMarkerStyle(21); gr1g->SetMarkerColor(80); gr1g->Draw("PL");

  canvas->cd(2); gr1g->SetTitle("No Noise"); gr1g->Draw("APL");
  canvas->cd(3); gr1gn1->SetTitle("1 MeV Noise"); gr1gn1->Draw("APL");
  canvas->cd(4); gr1gn3->SetTitle("3 MeV Noise"); gr1gn3->Draw("APL");
  
  

}
