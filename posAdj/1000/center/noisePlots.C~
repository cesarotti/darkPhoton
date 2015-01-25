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
  double reso500[14] = {19.82, 15.14, 11.86, 9.789, 8.455, 7.747, 6.837, 5.582, 4.925, 4.279, 2.395, 1.357, .703, .474};
  double reso500n1[14] = {19.75, 15.03, 11.81, 10.02, 8.457, 7.448, 6.923, 6.563, 6.098, 6.043, 5.802, 5.592, 5.695, 5.937};
  double reso500n3[14] = {19.26, 15.19, 12.78, 11.57, 11.29, 10.97, 10.59, 10.84, 10.79, 10.68, 11.28, 11.62, 11.88, 14.13};

  TGraph *gr500 = new TGraph (n, xtalNum, reso500);
  TGraph *gr500n1 = new TGraph (n, xtalNum, reso500n1);
  TGraph *gr500n3 = new TGraph (n, xtalNum, reso500n3);

  TCanvas *canvas = new TCanvas ("canvas", "Energy Resolution v. Crystal Num", 1000, 500);
  canvas->Divide(2, 2);
  canvas->cd(1);
  gr500n3->SetMarkerStyle(21); gr500n3->SetMarkerColor(55);
  gr500n3->Draw("APL");
  gr500n1->SetMarkerStyle(21); gr500n1->SetMarkerColor(13); gr500n1->Draw("PL");
  gr500->SetMarkerStyle(21); gr500->SetMarkerColor(80); gr500->Draw("PL");

  canvas->cd(2); gr500->SetTitle("No Noise"); gr500->Draw("APL");
  canvas->cd(3); gr500n1->SetTitle("1 MeV Noise"); gr500n1->Draw("APL");
  canvas->cd(4); gr500n3->SetTitle("3 MeV Noise"); gr500n3->Draw("APL");
  
  

}
