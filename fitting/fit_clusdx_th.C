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

// Global Params
Double_t yMin = -0.4;
Double_t yMax= 0.4;
Double_t xMin = 0.02;
Double_t xMax= 0.15;

// Names of the files we want to look at. I'm going to look at both sbs8 and sbs9 at the same time. 
TString FileName1="../output/sbs9_oct2_nofit_bigbins.root";
TString FileName2="../output/sbs8_oct2_nofit_bigbins.root";

// Name of the histogram we want to look at. 
TString HistName = "h_BBgr_clusxdiff_trth_allclus";

// To be overwritten
// Guesses for the linear fits on each to give the function somewhere to start
double MeanPredictionOffset = 0;
double MeanPredictionSlope = 1.35;

Double_t xMinMirror1 = -0.24;
Double_t xMaxMirror1= -0.17;
double MeanPredictionOffsetMirror1 = 0.15;
double MeanPredictionSlopeMirror1 = 1.25;

Double_t xMinMirror2 = -0.20;
Double_t xMaxMirror2= -0.05;
double MeanPredictionOffsetMirror2 = 0.005;
double MeanPredictionSlopeMirror2 = 1.42;

Double_t xMinMirror3 = 0.02;
Double_t xMaxMirror3= 0.15;
double MeanPredictionOffsetMirror3 = 0.009;
double MeanPredictionSlopeMirror3 = 1.35;

Double_t xMinMirror4 = 0.14;
Double_t xMaxMirror4= 0.21;
double MeanPredictionOffsetMirror4 = -0.24;
double MeanPredictionSlopeMirror4 = 1.28;

// Function Declarations 
Double_t linear(Double_t *x, Double_t *par);
Double_t trad(Double_t *x, Double_t *par);
Double_t gausFit(Double_t *x, Double_t *par);
TGraphErrors* FitBinsGausTH2D(TH2D* hist, TH1D* sigmahist, double xMin, double xMax, double MeanPredictionSlope, double MeanPredictionOffset);

// MAIN
void fit_clusdx_th(int mirrorNumber=3){


  if(mirrorNumber == 1){
    xMin = xMinMirror1;
    xMax = xMaxMirror1;
    MeanPredictionOffset = MeanPredictionOffsetMirror1;
    MeanPredictionSlope = MeanPredictionSlopeMirror1;
  }else if (mirrorNumber == 2){
    xMin = xMinMirror2;
    xMax = xMaxMirror2;
    MeanPredictionOffset = MeanPredictionOffsetMirror2;
    MeanPredictionSlope = MeanPredictionSlopeMirror2;
  }else if (mirrorNumber ==3){
    xMin = xMinMirror3;
    xMax = xMaxMirror3;
    MeanPredictionOffset = MeanPredictionOffsetMirror3;
    MeanPredictionSlope = MeanPredictionSlopeMirror3;
  }else if (mirrorNumber ==4){
    xMin = xMinMirror4;
    xMax = xMaxMirror4;
    MeanPredictionOffset = MeanPredictionOffsetMirror4;
    MeanPredictionSlope = MeanPredictionSlopeMirror4;
  }else{
    cout<<"error with mirror number input"<<endl;
    return;
  }

  

  TFile *f1 = TFile::Open(FileName1); // Load rootfile
  TFile *f2 = TFile::Open(FileName2); // Load rootfile

  // Set Up Canvas

  // Load Histograms
  TH2D* hist1 = (TH2D*)f1->Get(HistName);
  TH2D* hist2 = (TH2D*)f2->Get(HistName);


  TH2D* combinedhist = (TH2D*)(hist1->Clone("combinedhist"));
  combinedhist ->Add(hist2);
  

  TH1D* sigmahist =  new TH1D("sigmahist",";sigma from gaussian fit;",100,0,0.05);
  
  TH2D* origHist = (TH2D*)(combinedhist->Clone("origHist"));

  
  TH2D* clonedHist = (TH2D*)(origHist->Clone("clonedHist"));
  
   

    TCanvas* canvas =  new TCanvas("canvas","2D histo with gaussian fits");
    TGraphErrors* graph = new TGraphErrors();
    //// This function goes and fits each vertical slice of the TH2D over a given range. A prediction of what the fit is going to be is also needed. It returns a TGraphErrors that is the mean of each gaussian with the error bars being 1 sigma of those gausians 
    graph = FitBinsGausTH2D(clonedHist,sigmahist,xMin,xMax,MeanPredictionSlope, MeanPredictionOffset);
    // customize graph
    graph ->SetTitle("Vertical Bin Means with Gaussian fit");
    graph ->GetXaxis()->SetTitle("x axis");
    graph ->GetYaxis()->SetTitle("mean");
    graph ->SetMarkerStyle(6);
    combinedhist->Draw("colz");    graph ->Draw("P SAME");



    TF1 *fit_linear = new TF1("fit_linear", "pol1", xMin, xMax);
    fit_linear ->SetParameters(0,1.35);
    TFitResultPtr r = graph ->Fit("fit_linear", "MQR","",xMin,xMax);
    Int_t fitStatus = r;
    //if (fitStatus !=0) {cout<<"fit error" <<endl;}
    Double_t linear_offset = fit_linear->GetParameter(0);
    Double_t linear_slope = fit_linear->GetParameter(1);
    canvas ->Update();
    
    cout<<endl;
    cout<<"linear fit: "<<endl;
    cout<< "offset "<<linear_offset<<endl; 
    cout<<"slope "<< linear_slope<<endl;
    cout<<endl;    
  
    TCanvas* canvas2 =  new TCanvas("canvas2","sigmas");
    sigmahist ->Draw();


       
}//end main

//FUNCTIONS

Double_t gausFit(Double_t *x, Double_t *par)
{
  Double_t height = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];

  Double_t fit = height * exp(-0.5 * pow((x[0] - mean)/sigma,2));
  return fit;
}


TGraphErrors* FitBinsGausTH2D(TH2D* hist, TH1D* sigmahist, double xMin, double xMax, double MeanPredictionSlope = 1.35, double MeanPredictionOffset  = 0)
{
  // Define parameters
  const int nBinsX = hist->GetNbinsX();
  const double ymin = hist->GetYaxis()->GetXmin();
  const double ymax = hist->GetYaxis()->GetXmax();

  // Vectors to store fit results
  std::vector<double> binCenters;
  std::vector<double> means;
  std::vector<double> sigmas;
  std::vector<double> half_sigmas;
  std::vector<double> zeros;

  // Loop through the x bins and only consider the bins within the specified xMin to xMax range
  for (int i = 1; i <= nBinsX; i++) 
  {
    double binCenter = hist->GetXaxis()->GetBinCenter(i);
    
    // Skip bins outside the specified range
    if (binCenter < xMin || binCenter > xMax)
      continue;

    // Get 1D projection of the current bin in the Y direction
    TH1D* binHist = hist->ProjectionY("_py", i, i);
    
    // Define the Gaussian fit function
    TF1* gaussian = new TF1("gaussian", gausFit, ymin, ymax, 3);
    
    // Initial guesses for the fit parameters
    double mean_prediction = MeanPredictionSlope * binCenter + MeanPredictionOffset;
    gaussian->SetParameters(1000, mean_prediction, 0.015);
    gaussian->SetParLimits(1, mean_prediction - 0.1, mean_prediction + 0.1);
    gaussian->SetParLimits(2, 0.005, 0.1);

    // Perform the fit
    binHist->Fit(gaussian, "Q R");

    // Store fit results
    double mean = gaussian->GetParameter(1);
    double sigma = gaussian->GetParameter(2);
    
    binCenters.push_back(binCenter);
    means.push_back(mean);
    sigmas.push_back(sigma);
    sigmahist->Fill(sigma);
    half_sigmas.push_back(sigma * 0.5);
    zeros.push_back(0);

    // Clean up
    delete gaussian;
  }

  // Convert vectors to arrays suitable for TGraphErrors 
  int nPoints = binCenters.size();
  double* xValues = &binCenters[0];
  double* yValues = &means[0];
  double* yErrors = &sigmas[0];
  double* xErrors = &zeros[0];
  
  // Create the TGraphErrors
  TGraphErrors* graph = new TGraphErrors(nPoints, xValues, yValues, xErrors, yErrors);

  return graph;
}

