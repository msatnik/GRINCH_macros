#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include <fstream>
#include <string>
#include "TProfile.h"
#include "TGraphErrors.h"

// Global params

Double_t xmin_mirror1 = -0.25;
Double_t xmax_mirror1 = -0.12;
Double_t xmin_mirror2 = -0.2;
Double_t xmax_mirror2 = 0.01;
Double_t xmin_mirror3 = -0.05;
Double_t xmax_mirror3 = 0.15;
Double_t xmin_mirror4 = 0.09;
Double_t xmax_mirror4 = 0.22;

Double_t intercept_mirror1 = 0.170;
Double_t slope_mirror1=1.361;
Double_t intercept_mirror2 = 0.005;
Double_t slope_mirror2 = 1.424;
Double_t intercept_mirror3 = 0.0095;
Double_t slope_mirror3 = 1.339;
Double_t intercept_mirror4 = -0.270;
Double_t slope_mirror4 = 1.4658;

//// first try at fits
// Double_t intercept_mirror1 = 0.143;
// Double_t slope_mirror1=1.252;
// Double_t intercept_mirror2 = 0.006;
// Double_t slope_mirror2 = 1.432;
// Double_t intercept_mirror3 = 0.01;
// Double_t slope_mirror3 = 1.334;
// Double_t intercept_mirror4 = -0.244;
// Double_t slope_mirror4 = 1.323;



//// sbs9
Double_t SIGMA_mirror1 = 0.013;
Double_t SIGMA_mirror2 = 0.012;
Double_t SIGMA_mirror3 = 0.013;
Double_t SIGMA_mirror4 = 0.011;

//// sbs8
// Double_t SIGMA_mirror1 = 0.02;
// Double_t SIGMA_mirror2 = 0.0114; 
// Double_t SIGMA_mirror3 = 0.0130;
// Double_t SIGMA_mirror4 = 0.0107;

Double_t nSIGMA = 3;

// functions
Double_t linear(Double_t *x, Double_t *par);

// MAIN
void Plot_grinchdx_fits(){ 

  TFile *f1 = TFile::Open("../output/sbs14.root"); // Load rootfile
  TFile *f2 = TFile::Open("../output/sbs14.root"); // Load rootfile

  // Load Histograms
  TH2D* hist1 = (TH2D*)f1->Get("h_BBgr_clusxdiff_trth_allclus_mirror2");
  TH2D* hist2 = (TH2D*)f2->Get("h_BBgr_clusxdiff_trth_allclus_mirror2");  
  
  //combine histograms
  TH2D* combinedhist = (TH2D*)(hist1->Clone("combinedhist"));
  combinedhist ->Add(hist2);

  // set up linear fits
  TF1 *fit_mirror1 = new TF1("fit_mirror1", linear, xmin_mirror1,xmax_mirror1,2);
  fit_mirror1 ->SetParameters(intercept_mirror1,slope_mirror1);

  TF1 *fit_mirror2 = new TF1("fit_mirror2", linear, xmin_mirror2,xmax_mirror2,2);
  fit_mirror2 ->SetParameters(intercept_mirror2,slope_mirror2);

  TF1 *fit_mirror3 = new TF1("fit_mirror3", linear, xmin_mirror3,xmax_mirror3,2);
  fit_mirror3 ->SetParameters(intercept_mirror3,slope_mirror3);

  TF1 *fit_mirror4 = new TF1("fit_mirror4", linear, xmin_mirror4,xmax_mirror4,2);
  fit_mirror4 ->SetParameters(intercept_mirror4,slope_mirror4);



  // set up TGraphErrors 
  TH2D* hist = (TH2D*)(combinedhist->Clone("hist"));
  const int nBinsX = hist ->GetNbinsX();
  const int nBinsY = hist ->GetNbinsY();
  const int nBins =  nBinsX;
  const double xmin = hist ->GetXaxis()->GetXmin();
  const double xmax = hist ->GetXaxis()->GetXmax();
  const double ymin = hist ->GetYaxis()->GetXmin();
  const double ymax = hist ->GetYaxis()->GetXmax();
  // arrays to store fit results
  std::vector<double> binCenters(nBins);
  std::vector<double> means(nBins);
  std::vector<double> errorbars(nBins);
  std::vector<double> zeros(nBins);
  // loop throught the bins
  for (int i = 0; i<nBins; i++){
    // get bin center
    double binCenter = hist->GetXaxis()->GetBinCenter(i+1);
    double mean=0; 
    double errorbar = 0;
    if (binCenter >= xmin_mirror1 && binCenter <= xmax_mirror1)
      {
	mean = slope_mirror1*binCenter + intercept_mirror1;
	errorbar = nSIGMA*SIGMA_mirror1;
	//cout<<"mirror 1"<<endl;
      }
    else if (binCenter > xmin_mirror2 && binCenter <= xmax_mirror2)
      {
	mean = slope_mirror2*binCenter + intercept_mirror2;
	errorbar = nSIGMA*SIGMA_mirror2;
	//cout<<"mirror2"<<endl;
      }
    else if (binCenter > xmin_mirror3 && binCenter <= xmax_mirror3)
      {
	mean = slope_mirror3*binCenter + intercept_mirror3;
	errorbar = nSIGMA*SIGMA_mirror3;
	//cout<<"mirror3"<<endl;
      }
    else if (binCenter > xmin_mirror4 && binCenter <= xmax_mirror4)
      {
	mean = slope_mirror4*binCenter + intercept_mirror4;
	errorbar = nSIGMA*SIGMA_mirror4;
	//cout<<"mirror4"<<endl;
      }
    else if (binCenter < xmin_mirror1 || binCenter > xmax_mirror4){ 
      mean = 0; 
      errorbar=0;
      // cout<< "ouside bounds" <<endl;
    }
    else { 
      cout<<"didn't code logic right in the loop: binCenter = "<<binCenter<<endl;
    }
    //cout<<"mean "<<mean<<endl;
      means[i]=mean;
      errorbars[i] = errorbar;
      zeros[i]=0;
      binCenters[i] = binCenter;
  
  }
 // convert vectors to arrays suitable for TGraphErrors 
  int nPoints = nBins;
  double* xValues =  &binCenters[0];
  double* yValues = &means[0];
  double* yErrors = &errorbars[0];
  double* xErrors = &zeros[0];
  // make TGraphErrors
  TGraphErrors* graph =  new TGraphErrors(nPoints, xValues, yValues, xErrors, yErrors);
  graph ->SetMarkerStyle(6);
  // done setting up TGraphErrors





  TCanvas* canvas = new TCanvas("canvas","dx vs th");
  combinedhist ->SetContour(256);
  combinedhist->Draw("colz");
  fit_mirror1 ->Draw("same");
  fit_mirror2 ->Draw("same");
  fit_mirror3->Draw("same");
  fit_mirror4 ->Draw("same");
  graph->Draw("same");

}//end MAIN

Double_t linear(Double_t *x, Double_t *par)
{
  Double_t fit =  par[0] + par[1]*x[0];
  return fit;
}
