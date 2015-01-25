#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <sstream> 

void dfsPlots100()
/*
  This is a plot for the clustering heuristic of ordering the energies, then choosing how many crystals to count. 
*/
{
  int n = 10;
  double xtalNum[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  double reso100[10] = {8.516, 4.64, 3.051, 2.012, 1.769, 1.571, 1.386, 1.314, 1.174, 1.180};
  double reso100n3[10] = {9.337, 4.522, 3.549, 2.644, 2.495, 2.317, 2.121, 2.015, 1.842, 1.862};
  double reso100n1[10] = {8.309, 4.236,  2.606, 2.094, 1.725, 1.370, 1.228, 1.106, 1.02, .969};

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
