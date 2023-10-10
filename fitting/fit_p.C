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
Double_t yMin = -7;
Double_t yMax= 10;
Double_t xMin = 8 ;
Double_t xMax= 16;

// Function Declarations 
Double_t poly2(Double_t *x, Double_t *par);
Double_t trad(Double_t *x, Double_t *par);
Double_t gausFit(Double_t *x, Double_t *par);
TGraphErrors* FitBinsGausTH2D(TH2D* hist);

// MAIN
void fit_clusdx_th(){ 

  TFile *f1 = TFile::Open("../output/sbs9_oct2_nofit.root"); // Load rootfile
  TFile *f2 = TFile::Open("../output/sbs8_oct2_nofit.root"); // Load rootfile

  // Set Up Canvas

  // Load Histograms
  TH2D* hist1 = (TH2D*)f1->Get("h_BBgr_clusxdiff_trth_allclus_mirror2");
  TH2D* hist2 = (TH2D*)f2->Get("h_BBgr_clusxdiff_trth_allclus_mirror2");

  TH2D* hist3 = (TH2D*)f1->Get("h_BBgr_clusxdiff_trth_allclus_mirror3");
  TH2D* hist4 = (TH2D*)f2->Get("h_BBgr_clusxdiff_trth_allclus_mirror3");


  // TH2D* hist3 = (TH2D*)f1->Get("h_BBtr_p_th_cut_mirror3");
  // TH2D* hist4 = (TH2D*)f2->Get("h_BBtr_p_th_cut_mirror3");


  TH2D* combinedhist = (TH2D*)(hist1->Clone("combinedhist"));
  combinedhist ->Add(hist2);
  TH2D* combinedhist2 = (TH2D*)(hist3->Clone("combinedhist2"));
  combinedhist2 ->Add(hist4);
  TH2D* combinedhist3 = (TH2D*)(hist1->Clone("combinedhist3"));
  combinedhist3 ->Add(hist2);
  combinedhist3 ->Add(hist3);
  combinedhist3 ->Add(hist4);
  
  
 //  // Clone Histogram and limit Y range. (Maybe should also limit x range?) 
//   // Could maybe use a clone command instead of doing it manually like this 
//   TH2D* clonedHist = new TH2D("clonedHist", origHist->GetTitle(),
//                                 origHist->GetXaxis()->GetNbins(),
//                                 origHist->GetXaxis()->GetXmin(),
//                                 origHist->GetXaxis()->GetXmax(),
//                                 origHist->GetYaxis()->FindBin(yMax) - origHist->GetYaxis()->FindBin(yMin) + 1,
//                                 yMin, yMax);

//  // Loop over the X and Y bins of the original TH2D
//     for (int i = 1; i <= origHist->GetXaxis()->GetNbins(); i++) {
//         for (int j = origHist->GetYaxis()->FindBin(yMin); j <= origHist->GetYaxis()->FindBin(yMax); j++) {
//             double binContent = origHist->GetBinContent(i, j);
//             clonedHist->SetBinContent(i, j - origHist->GetYaxis()->FindBin(yMin) + 1, binContent);
//         }
//     }

//  TH2D* clonedHist2 = new TH2D("clonedHist2", origHist->GetTitle(),
//                                 origHist->GetXaxis()->GetNbins(),
//                                 origHist->GetXaxis()->GetXmin(),
//                                 origHist->GetXaxis()->GetXmax(),
//                                 origHist->GetYaxis()->FindBin(yMax) - origHist->GetYaxis()->FindBin(yMin) + 1,
//                                 yMin, yMax);

//  // Loop over the X and Y bins of the original TH2D
//     for (int i = 1; i <= origHist->GetXaxis()->GetNbins(); i++) {
//         for (int j = origHist->GetYaxis()->FindBin(yMin); j <= origHist->GetYaxis()->FindBin(yMax); j++) {
//             double binContent = origHist->GetBinContent(i, j);
//             clonedHist2->SetBinContent(i, j - origHist->GetYaxis()->FindBin(yMin) + 1, binContent);
//         }
//     }

// TH2D* clonedHist3 = new TH2D("clonedHist3", origHist->GetTitle(),
//                                 origHist->GetXaxis()->GetNbins(),
//                                 origHist->GetXaxis()->GetXmin(),
//                                 origHist->GetXaxis()->GetXmax(),
//                                 origHist->GetYaxis()->FindBin(yMax) - origHist->GetYaxis()->FindBin(yMin) + 1,
//                                 yMin, yMax);

//  // Loop over the X and Y bins of the original TH2D
//     for (int i = 1; i <= origHist->GetXaxis()->GetNbins(); i++) {
//         for (int j = origHist->GetYaxis()->FindBin(yMin); j <= origHist->GetYaxis()->FindBin(yMax); j++) {
//             double binContent = origHist->GetBinContent(i, j);
//             clonedHist3->SetBinContent(i, j - origHist->GetYaxis()->FindBin(yMin) + 1, binContent);
//         }
//     }

//   //  TH2D *clonedHist = (TH2D*)(origHist->Clone("clonedHist"));
//   //  TH2D *clonedHist2 = (TH2D*)(origHist->Clone("clonedHist2"));
   

    TCanvas* canvas =  new TCanvas("canvas","2D histo with gaussian fits");
    // combinedhist ->Draw("colz");
    //combinedhist2 ->Draw("colz same");
    combinedhist3 ->Draw("colz");


    // TGraphErrors* graph = new TGraphErrors();
    // // TH2D *testHist = (TH2D*)(origHist->Clone("clonedHist"));
    // graph = FitBinsGausTH2D(clonedHist3);
    // // customize graph
    // graph ->SetTitle("Vertical Bin Means with Gaussian fit");
    // graph ->GetXaxis()->SetTitle("x axis");
    // graph ->GetYaxis()->SetTitle("mean");
    // graph ->SetMarkerStyle(20);
    // graph ->SetMarkerSize(1.2);
    // //draw graph
    // clonedHist3->Draw("colz");
    // //fit_poly2 ->Draw("SAME");
    // graph ->Draw("P SAME");



    // TF1 *fit_poly2 = new TF1("fit_poly2", poly2, xMin, xMax, 2);
    // fit_poly2 ->SetParameters(4,-0.2);
    // TFitResultPtr r = graph ->Fit("fit_poly2", "MQR","",xMin,xMax);
    // Int_t fitStatus = r;
    // //if (fitStatus !=0) {cout<<"fit error" <<endl;}
    // Double_t poly_offset = fit_poly2->GetParameter(0);
    // Double_t poly_slope = fit_poly2->GetParameter(1);
    // canvas ->Update();
    
    // cout<<"poly fit: "<<endl;
    // cout<< "offset "<<poly_offset<<endl; 
    // cout<<"slope "<< poly_slope<<endl;
    // cout<<endl; 
    
  


       
}//end main

//FUNCTIONS

Double_t poly2(Double_t *x, Double_t *par)
{
  Double_t fit =  par[0] + par[1]*x[0];
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
      TF1* gaussian =  new TF1("gaussian","gaus",xmin,xmax); // this uses the default "gaus" from root. Could potentially be worth it to make your own gaus function as sometimes the default can be weird.
      binHist->Fit(gaussian, "Q");// "quiet" option
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
}// end TGraphErrors* FitBinsGausTH2D(TH2D* hist)


