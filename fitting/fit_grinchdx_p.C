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
Double_t yMin = -0.2;
Double_t yMax= 0.2;
Double_t xMin = 2.4;
Double_t xMax= 3.5;

// Function Declarations 
Double_t linear(Double_t *x, Double_t *par);
Double_t quadratic(Double_t *x, Double_t *par);
Double_t trad(Double_t *x, Double_t *par);
Double_t gausFit(Double_t *x, Double_t *par);
TGraphErrors* FitBinsGausTH2D(TH2D* hist);

// MAIN
void fit_grinchdx_p(){ 

  TFile *f1 = TFile::Open("../output/sbs8_sept28_bigmirrorcuts.root"); // Load rootfile

  // Set Up Canvas

  // Load Histograms
  TH2D *origHist = (TH2D*)f1->Get("h_BBgr_clusxdiff_trp_allclus_mirror3");
  
   
    TCanvas* canvas =  new TCanvas("canvas","2D histo with gaussian fits");
    TGraphErrors* graph = new TGraphErrors();
    // TH2D *testHist = (TH2D*)(origHist->Clone("clonedHist"));
    graph = FitBinsGausTH2D(origHist);
    // customize graph
    graph ->SetTitle("Vertical Bin Means with Gaussian fit");
    graph ->GetXaxis()->SetTitle("x axis");
    graph ->GetYaxis()->SetTitle("mean");
    graph ->SetMarkerStyle(6);
    graph ->SetMarkerSize(1);
    //draw graph
    origHist->Draw("colz");
    //fit_poly2 ->Draw("SAME");
    graph ->Draw("P SAME");

    

    // TF1 *fit_linear = new TF1("fit_linear", linear, xMin, xMax, 2);
    // fit_linear ->SetParameters(0,-0.5);
    // fit_linear ->SetParLimits(1,-1,0.5);
    // fit_linear ->SetParLimits(0,-0.1,0.1);
    // TFitResultPtr r = graph ->Fit("fit_linear", "MQR","",xMin,xMax);
    // Int_t fitStatus = r;
    // //if (fitStatus !=0) {cout<<"fit error" <<endl;}
    // Double_t linear_offset = fit_linear->GetParameter(0);
    // Double_t linear_slope = fit_linear->GetParameter(1);
    // canvas ->Update();
    // cout<<"linear fit: "<<endl;
    // cout<< "offset "<<linear_offset<<endl; 
    // cout<<"slope "<< linear_slope<<endl;
    // cout<<endl; 


 
    
  

    TF1 *fit_trad = new TF1("fit_trad", trad, xMin, xMax, 2);
    fit_trad ->SetParameters(-0.12,0.2);
    // fit_trad ->SetLineColor(kBlue);
    TFitResultPtr r1 = graph ->Fit("fit_trad", "MQR","",xMin,xMax);
    Double_t trad_offset =  fit_trad->GetParameter(0);
    Double_t trad_p1 = fit_trad->GetParameter(1);
    cout<<"traditional fit :"<<endl;
    cout<<"offset " << trad_offset<<endl; 
    cout<< "p1 "<<trad_p1<<endl; 
    

       
}//end main

//FUNCTIONS

Double_t linear(Double_t *x, Double_t *par)
{
  Double_t fit =  par[0] + par[1]*x[0];
  return fit;
}

Double_t quadratic(Double_t *x, Double_t *par)
{
  Double_t fit =  par[0] + par[1]*x[0] + par[2]*x[0];
  return fit;
}

Double_t exp_decay(Double_t *x, Double_t *par)
{
  Double_t fit =  par[0] + par[1]*x[0] + par[2]*x[0];
  return fit;
}

Double_t trad(Double_t *x, Double_t *par)
{
  Double_t fit =  par[0] +  par[1]* pow(x[0],-0.5);
  return fit;
}


Double_t gausFit(Double_t *x, Double_t *par)
{
  Double_t height = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];

  Double_t fit = height * exp(-0.5 * pow((x[0] - mean)/sigma,2));
  return fit;
}


TGraphErrors* FitBinsGausTH2D(TH2D* hist)
//// Fits each vertical bin in a TH2D to a gaussian and returns a TGraphErros with the means and the sigma as the vertical error bars. 
{
  // Fit each vertical bin to a gaussian and plot the means on top of the 2D histo. 
  // Define parameters
  const int nBinsX = hist ->GetNbinsX();
  const int nBinsY = hist ->GetNbinsY();
  const int nBins =  nBinsX;
  const double xmin = hist ->GetXaxis()->GetXmin();
  const double xmax = hist ->GetXaxis()->GetXmax();
  const double ymin = hist ->GetYaxis()->GetXmin();
  const double ymax = hist ->GetYaxis()->GetXmax();

  // cout<<endl;
  // cout<<"nBins: "<<nBins<<endl;
  // cout<<endl;

  // arrays to store fit results
  std::vector<double> binCenters(nBins);
  std::vector<double> means(nBins);
  std::vector<double> sigmas(nBins);
  std::vector<double> half_sigmas(nBins);
  std::vector<double> zeros(nBins);

  //loop through the bins
  for (int i = 0; i<nBins; i++) // 
    {
      // get 1d histo of y projection
      TH1D* binHist =  hist->ProjectionY("_py",i +1, i+1);
      // get bin center
      double binCenter = hist->GetXaxis()->GetBinCenter(i+1);
      // define fit and gaussian fcn
      TF1* gaussian =  new TF1("gaussian",gausFit,ymin,ymax,3); // this uses the default "gaus" from root. Could potentially be worth it to make your own gaus function as sometimes the default can be weird
      //cout<<"ymin "<<ymin<<" ymax: "<<ymax <<endl;
       double mean_prediction = 0*binCenter;
      gaussian ->SetParameters(1000,mean_prediction,0.01);
      gaussian->SetParLimits(1,mean_prediction - 0.05, mean_prediction +0.05);
      gaussian ->SetParLimits(2,0.001,0.05);

      binHist->Fit(gaussian, "Q");// "quiet" option
      //cout<<"bincenter: " <<binCenter<< " mean pred: "<<mean_prediction<<endl;
      //get fit results
      double mean = gaussian -> GetParameter(1);
      double sigma = gaussian ->GetParameter(2);
      //store results 
      binCenters[i] = binCenter;
      means[i] = mean;
      sigmas[i] = sigma;
      half_sigmas[i] = sigma*0.5;
      zeros[i] = 0;
      // cout<<"bin center: "<<binCenter<<endl;
      // cout<<"mean: "<<mean<<endl;
      // cout<<"sigma: "<<sigma<<endl;
      // cout<<endl;
      //// clean up
      delete gaussian;
    }// end loop over bins

  // convert vectors to arrays suitable for TGraphErrors 
  int nPoints = nBins;
  double* xValues =  &binCenters[0];
  double* yValues = &means[0];
  double* yErrors = &half_sigmas[0];
  double* xErrors = &zeros[0];
    
  // make TGraphErrors
  TGraphErrors* graph =  new TGraphErrors(nPoints, xValues, yValues, xErrors, yErrors);
    
  return graph;
}


