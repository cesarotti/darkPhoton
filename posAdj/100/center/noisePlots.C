#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <sstream> 

void noisePlots()
/*
  This is a plot for the clustering heuristic of ordering the energies, then choosing how many crystals to count. 
*/
{
  int n = 13;
  double xtalNum[13] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 25, 50};
  double reso100[13] = {-18.81, -10.4, -6.81, -4.713, -3.39, -2.448, -1.81, -1.353, -1.046, -0.8456, -.583, -.583, -.578};
  double reso100n3[13] = {-18.81, -10.38, -6.8, -4.716, -3.409, -2.49, -1.898, -1.402, -1.135, -.9051, -.634, -.642, -.6382};
  double reso100n1[13] = {-18.79, -10.38, -6.81, -4.728, -3.395, -2.458, -1.827, -1.397, -1.086, -.863, -.6285, -.601, -.573};

  TGraph *gr100 = new TGraph (n, xtalNum, reso100);
  TGraph *gr100n3 = new TGraph (n, xtalNum, reso100n3);
  TGraph *gr100n1 = new TGraph (n, xtalNum, reso100n1);

  TCanvas *canvas = new TCanvas ("canvas", "Energy Resolution v. Crystal Num", 1000, 500);
  canvas->Divide(2,2);
  canvas->cd(1);
  gr100n3->SetMarkerStyle(21); gr100n3->SetMarkerColor(100);
  gr100n3->Draw("APL");
  canvas->cd(1);
  gr100n1->SetMarkerStyle(21); gr100n1->SetMarkerColor(55);  gr100n1->Draw("PL");canvas->cd(1);
  gr100->SetMarkerStyle(21); gr100->SetTitle("All Curves"); gr100->Draw("PL");


  canvas->cd(2); gr100->SetMarkerStyle(21); gr100->SetTitle("No Noise");
  gr100->GetXaxis()->SetTitle("# of Crystals in Cluster"); 
  gr100->GetYaxis()->SetTitle("Difference in resolved v. true"); gr100->Draw("APL");
  canvas->cd(3); gr100n1->SetMarkerStyle(21); gr100n1->SetMarkerColor(55);
  gr100n1->SetTitle("1 MeV Noise"); gr100n1->Draw("APL");
  canvas->cd(4); gr100n3->SetMarkerStyle(21); gr100n3->SetMarkerColor(100);
  gr100n3->SetTitle("3 MeV Noise"); gr100n3->Draw("APL");


}
