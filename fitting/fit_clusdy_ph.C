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
Double_t xMin = -0.06;
Double_t xMax= 0.06;

TString FileName1="../output/sbs9_oct9.root";
TString FileName2="../output/sbs8_oct9.root";

TString HistName = "h_BBgr_clusydiff_trph_allclus_mirror2";
int PolyOrder = 3;

TString HistNameMirror1 =  "h_BBgr_clusydiff_trph_allclus_mirror1";
int PolyOrderMirror1 = 3;

TString HistNameMirror2 =  "h_BBgr_clusydiff_trph_allclus_mirror2";
int PolyOrderMirror2 = 3;

TString HistNameMirror3 =  "h_BBgr_clusydiff_trph_allclus_mirror3";
int PolyOrderMirror3 = 3;

TString HistNameMirror4 =  "h_BBgr_clusydiff_trph_allclus_mirror4";
int PolyOrderMirror4 = 3;


double MeanPredictionSlope = -2.44;
double MeanPredictionOffset = 0;

// Function Declarations 
Double_t linear(Double_t *x, Double_t *par);
Double_t trad(Double_t *x, Double_t *par);
Double_t gausFit(Double_t *x, Double_t *par);
Double_t poly2(Double_t *x, Double_t *par);
Double_t poly3(Double_t *x, Double_t *par);
Double_t poly4(Double_t *x, Double_t *par);
Double_t poly5(Double_t *x, Double_t *par);
TGraphErrors* FitBinsGausTH2D(TH2D* hist, TH1D* sigmahist, double xMin, double xMax, double MeanPredictionSlope, double MeanPredictionOffset);

// MAIN
void fit_clusdy_ph(int MirrorNumber =3){

  if(MirrorNumber == 1){
    HistName = HistNameMirror1;
    PolyOrder = PolyOrderMirror1;
  }else if(MirrorNumber == 2){
    HistName = HistNameMirror2;
    PolyOrder = PolyOrderMirror2;
  }else if(MirrorNumber == 3){
    HistName = HistNameMirror3;
    PolyOrder = PolyOrderMirror3;
  }else if(MirrorNumber == 4){
    HistName = HistNameMirror4;
    PolyOrder = PolyOrderMirror4;
  }else{
    cout<<"error with mirror number input"<<endl;
    return;
  }


  TFile *f1 = TFile::Open(FileName1); // Load rootfile
  TFile *f2 = TFile::Open(FileName2); // Load rootfile

 
  // Load Histograms
  TH2D* hist1 = (TH2D*)f1->Get(HistName);
  TH2D* hist2 = (TH2D*)f2->Get(HistName);

  // Add the histos from the two kinematics together
  TH2D* combinedhist = (TH2D*)(hist1->Clone("combinedhist"));
  combinedhist ->Add(hist2);

  

  TH1D* sigmahist =  new TH1D("sigmahist",";sigma from gaussian fit;",25,0,0.05);
  
  

    TCanvas* canvas =  new TCanvas("canvas","2D histo with gaussian fits");


    TGraphErrors* graph = new TGraphErrors();
    //// This function goes and fits each vertical slice of the TH2D over a given range. A prediction of what the fit is going to be is also needed. It returns a TGraphErrors that is the mean of each gaussian with the error bars being 1 sigma of those gausians 
    graph = FitBinsGausTH2D(combinedhist,sigmahist, xMin, xMax, MeanPredictionSlope, MeanPredictionOffset);
    // customize graph
    graph ->SetTitle("Vertical Bin Means with Gaussian fit");
    graph ->GetXaxis()->SetTitle("x axis");
    graph ->GetYaxis()->SetTitle("mean");
    graph ->SetMarkerStyle(6);
    //graph ->SetMarkerSize(1.2);
    //draw graph
    combinedhist->Draw("colz");
    graph ->Draw("P SAME");


    TF1 *fit_poly3 = new TF1("fit_poly3", "pol3", xMin, xMax);
    fit_poly3 ->SetParameters(0.04,-1.7,0.9,0);
    //fit_poly3 ->SetParLimits(2,0,0);
    // fit_poly3 ->SetParLimits(2,-10,10);
    // fit_poly3 ->SetParLimits(1,-10,10);
    TFitResultPtr r = graph ->Fit("fit_poly3", "MQR","",xMin,xMax);
    Int_t fitStatus = r;
    //if (fitStatus !=0) {cout<<"fit error" <<endl;}
    Double_t poly3_par0 = fit_poly3->GetParameter(0);
    Double_t poly3_par1 = fit_poly3->GetParameter(1);
    Double_t poly3_par2 = fit_poly3->GetParameter(2);
    Double_t poly3_par3 = fit_poly3->GetParameter(3);
    Double_t poly3_par0_error = fit_poly3->GetParError(0);
    Double_t poly3_par1_error = fit_poly3->GetParError(1);
    Double_t poly3_par2_error = fit_poly3->GetParError(2);
    Double_t poly3_par3_error = fit_poly3->GetParError(3);

    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
   // Add entries to the legend for the fit parameters
    legend->AddEntry(fit_poly3, "3rd order poly", "l");
    legend->AddEntry("", Form("par0 = %.4f +- %.4f ", poly3_par0, poly3_par0_error), "");
    legend->AddEntry("", Form("par1 = %.4f +- %.4f ", poly3_par1, poly3_par1_error), "");
    legend->AddEntry("", Form("par2 = %.4f +- %.4f ", poly3_par2, poly3_par2_error), "");
    legend->AddEntry("", Form("par3 = %.4f +- %.4f ", poly3_par3, poly3_par3_error), "");
    legend->Draw();
    canvas ->Update();
    cout<<endl;
    cout<<"3rd order poly fit: "<<endl;
    cout<< "par0 "<<poly3_par0<<endl; 
    cout<<"par1 "<< poly3_par1<<endl;
    cout<< "par2 "<<poly3_par2<<endl; 
    cout<<"par3 "<< poly3_par3<<endl;
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
    gaussian->SetParameters(1000, mean_prediction, 0.02);
    gaussian->SetParLimits(1, mean_prediction - 0.1, mean_prediction + 0.1);
    gaussian->SetParLimits(2, 0.005, 0.1);

    // Perform the fit
    binHist->Fit(gaussian, "Q R");

    // Store fit results
    double mean = gaussian->GetParameter(1);
    double sigma = gaussian->GetParameter(2);

    ////  if the stdev got filled to the smallest value possible, it probably didn't fit right. Going to make it bigger so that it doesn't affect the fit so much. 
    if (sigma == 0.005){
      sigma = 0.1;
    }
    
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

