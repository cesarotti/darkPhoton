#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <sstream> 

void dfsPlots1000()

{
  int n = 15;
  double xtalNum[15] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 50};
  double reso1g[15] = {31.98, 23.65, 18.93, 16.55, 14.98, 13.04, 11.05, 9.385, 9.070, 7.653, 4.527, 3.476, 2.517, 2.723, 4.993};
  double reso1gn1[15] = {36.67, 26.45, 19.31, 17.25, 16.78, 13.55, 12.18, 11.00, 10.44, 8.671, 4.912, 3.324, 2.068, 1.410, 1.146};
  double reso1gn3[15] = {34.09, 24.80, 20.46, 18.67, 16.08, 13.83, 12.59, 11.16, 9.525, 8.640, 5.26, 3.372, 2.829, 2.374, 2.251};

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
