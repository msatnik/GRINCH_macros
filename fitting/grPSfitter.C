#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "TMath.h"

#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <math.h>

Double_t poly(Double_t *x, Double_t *par);
Double_t Gaus(Double_t *x, Double_t *par);
Double_t poly_Gaus(Double_t *x, Double_t *par);
Double_t skewed_gaus(Double_t *x, Double_t *par);
Double_t poly_skewed_gaus(Double_t *x, Double_t *par);
Double_t poly5_skewed_gaus(Double_t *x, Double_t *par);
Double_t poly5(Double_t *x, Double_t *par);
Double_t poly4(Double_t *x, Double_t *par);
Double_t poly4_skewed_gaus(Double_t *x, Double_t *par);
Double_t poly_2Gaus(Double_t *x, Double_t *par);
Double_t LandauTest(Double_t *x, Double_t *par);
Double_t Landau(Double_t *x, Double_t *par);
Double_t poly4_2Landau(Double_t *x, Double_t *par);
Double_t Landau2(Double_t *x, Double_t *par);
Double_t poly4_Landau_Gaus(Double_t *x, Double_t *par);
Double_t skewed_2gaus(Double_t *x, Double_t *par);


// class CustomFit : public TF1{
// public: CustomFit() : TF1("CustomFit", gaus, 0,0.2,3);
// };


void grPSfitter(){// MAIN
  TFile *f1 = TFile::Open("../output/sbs9_gr30_ps20.root"); // Load rootfile

  // Set up Canvas
  TCanvas *c1 = new TCanvas();
  //c1->Divide(2,2);

  // Load Histograms
  TH1D *h_BBps_e = (TH1D*)f1->Get("h_BBps_e");
  TH1D *h_BBps_grcut = (TH1D*)f1->Get("h_BBps_grcut");
  TH1D *h_BBps_granticut = (TH1D*)f1->Get("h_BBps_granticut");

  
  // initial values for fits
  Double_t par0_init = 500;
  Double_t par1_init = -10;
  Double_t par2_init = 10;
  Double_t par3_init = -10;
  Double_t par4_init = 10;
  Double_t par5_init = -10;
  Double_t par6_init = 10;
  Double_t height_init = 5695;
  Double_t mean_init = 0.09;
  Double_t sigma_init = 0.02;
  Double_t alpha_init = 0.5;
  Double_t height2_init = 1763;
  Double_t mean2_init = 0.52;
  Double_t sigma2_init = 0.23;
  Double_t alpha2_init = 0.5;
						
  Double_t constant_init = 30900;
  Double_t mpv_init =  0.08;
  Double_t sigma_landau_init = 0.012;
  Double_t constant2_init = 10040;
  Double_t mpv2_init =  0.46;
  Double_t sigma2_landau_init = 0.17;


  TF1 *fit_landau_test =  new TF1("fit_landau_test", LandauTest, 0.02, 1.4 , 3);
  fit_landau_test ->SetParameters(30900, 0.08, 0.012);

  TF1 *fit_landau_gaus_test =  new TF1("fit_landau_gaus_test", poly4_Landau_Gaus, 0.02, 1.4, 11);
  fit_landau_gaus_test ->SetParameters(0, 0, 0, 0,0,constant_init, mpv_init,sigma_landau_init, height2_init, mean2_init, sigma2_init);

  TF1* fit_landau2_test = new TF1("fit_landau2_test", Landau2, 0.02,1.4,6);
  fit_landau2_test ->SetParameters(constant_init, mpv_init,sigma_landau_init, constant2_init, mpv2_init, sigma2_landau_init);

    TF1 *fit_poly4_landau_test =  new TF1("fit_poly4_landau_test", poly4_2Landau, 0.02, 1.4, 11);
    fit_poly4_landau_test ->SetParameters(0, 0, 0, 0,0,constant_init, mpv_init,sigma_landau_init, constant2_init, mpv2_init, sigma2_landau_init);

    TF1 *fit_gaus_test =  new TF1("fit_gaus_test", poly_2Gaus, 0.02, 1.4, 11);
    fit_gaus_test ->SetParameters(0,0,0,0,0,height_init, mean_init, sigma_init, height2_init, mean2_init, sigma2_init);

    TF1 *fit_skewed_2gaus_test =  new TF1("fit_skewed_2gaus_test",skewed_2gaus, 0.02, 1.4,8);
    fit_skewed_2gaus_test ->SetParameters(height_init, mean_init, sigma_init, alpha_init, height2_init, mean2_init, sigma2_init, alpha2_init);


  TF1 *fit_BBps_granticut = new TF1("fit_BBps_granticut",poly_2Gaus, 0.02,1.4,11);
  fit_BBps_granticut ->SetParameters(0, 0, 0, 0,0,height_init, mean_init,sigma_init, height2_init, mean2_init, sigma2_init);
  TFitResultPtr r = h_BBps_granticut ->Fit(fit_BBps_granticut,"RQ0");
  Int_t fitStatus = r;
  if (fitStatus !=0) {cout<<"fit error poly_2Gaus" <<endl;}
  Double_t par0 = fit_BBps_granticut ->GetParameter(0);
  Double_t par1 = fit_BBps_granticut ->GetParameter(1);
  Double_t par2 = fit_BBps_granticut ->GetParameter(2);
  Double_t par3 = fit_BBps_granticut ->GetParameter(3);
  Double_t par4 = fit_BBps_granticut ->GetParameter(4);
  Double_t par0_err = fit_BBps_granticut ->GetParError(0);
  Double_t par1_err = fit_BBps_granticut ->GetParError(1);
  Double_t par2_err = fit_BBps_granticut ->GetParError(2);
  Double_t par3_err = fit_BBps_granticut ->GetParError(3);
  Double_t par4_err = fit_BBps_granticut ->GetParError(4);
  Double_t height= fit_BBps_granticut ->GetParameter(5);
  Double_t mean = fit_BBps_granticut ->GetParameter(6);
  Double_t sigma = fit_BBps_granticut ->GetParameter(7);
  Double_t height_err= fit_BBps_granticut ->GetParError(5);
  Double_t mean_err = fit_BBps_granticut ->GetParError(6);
  Double_t sigma_err = fit_BBps_granticut ->GetParError(7);
  Double_t height2= fit_BBps_granticut ->GetParameter(8);
  Double_t mean2 = fit_BBps_granticut ->GetParameter(9);
  Double_t sigma2= fit_BBps_granticut ->GetParameter(10);
  Double_t height2_err= fit_BBps_granticut ->GetParError(8);
  Double_t mean2_err = fit_BBps_granticut ->GetParError(9);
  Double_t sigma2_err = fit_BBps_granticut ->GetParError(10);

  
  TF1 *fit_result_poly_BBps_granticut = new TF1("fit_result_poly_BBps_granticut", poly4, 0.02,1.4,5);
  fit_result_poly_BBps_granticut ->SetParameters(par0,par1,par2,par3,par4);
  TF1 *fit_result_Gaus1_BBps_granticut = new TF1("fit_result_Gaus1_BBps_granticut", Gaus, 0.02,0.3,3);
  fit_result_Gaus1_BBps_granticut ->SetParameters(height, mean, sigma);
  TF1 *fit_result_Gaus2_BBps_granticut = new TF1("fit_result_Gaus2_BBps_granticut", Gaus, 0.02,1.4,3);
  fit_result_Gaus2_BBps_granticut ->SetParameters(height2, mean2, sigma2);

  fit_result_poly_BBps_granticut ->SetLineColor(kMagenta - 6);
  fit_result_Gaus1_BBps_granticut ->SetLineColor(kViolet - 6);
  fit_result_Gaus2_BBps_granticut ->SetLineColor(kGreen -6);
  fit_BBps_granticut ->SetLineColor(kRed);

  h_BBps_granticut ->Draw("hist");
  fit_result_poly_BBps_granticut ->Draw("same");
  fit_result_Gaus1_BBps_granticut ->Draw("same");
  fit_result_Gaus2_BBps_granticut ->Draw("same");
  fit_BBps_granticut ->Draw("same");



    TF1 *fit_BBps_granticut_landau = new TF1("fit_BBps_granticut_landau",poly4_2Landau, 0.02,1.4,11);
    fit_BBps_granticut_landau ->SetParameters(par0_init, par1_init, par2_init, par3_init,par4_init, constant_init, mpv_init, sigma_landau_init, constant2_init, mpv2_init, sigma2_landau_init);
    TH1D *h_BBps_granticut_copy =  (TH1D*)(h_BBps_granticut->Clone("h_BBps_granticut_copy"));
  TFitResultPtr r2 = h_BBps_granticut_copy ->Fit(fit_BBps_granticut_landau,"RQ0");
  Int_t fitStatus2 = r2;
  if (fitStatus2 !=0) {cout<<"fit error poly4_2Landau" <<endl;}
  Double_t par0_landau = fit_BBps_granticut_landau ->GetParameter(0);
  Double_t par1_landau = fit_BBps_granticut_landau ->GetParameter(1);
  Double_t par2_landau = fit_BBps_granticut_landau->GetParameter(2);
  Double_t par3_landau = fit_BBps_granticut_landau->GetParameter(3);
  Double_t par4_landau = fit_BBps_granticut_landau->GetParameter(4);
  Double_t par0_landau_err = fit_BBps_granticut_landau ->GetParError(0);
  Double_t par1_landau_err = fit_BBps_granticut_landau->GetParError(1);
  Double_t par2_landau_err = fit_BBps_granticut_landau ->GetParError(2);
  Double_t par3_landau_err = fit_BBps_granticut_landau ->GetParError(3);
  Double_t par4_landau_err = fit_BBps_granticut_landau ->GetParError(4);
  Double_t constant= fit_BBps_granticut_landau ->GetParameter(5);
  Double_t mpv = fit_BBps_granticut_landau ->GetParameter(6);
  Double_t sigma_landau = fit_BBps_granticut_landau ->GetParameter(7);
  Double_t constant_err= fit_BBps_granticut_landau ->GetParError(5);
  Double_t mpv_err = fit_BBps_granticut_landau->GetParError(6);
  Double_t sigma_landau_err = fit_BBps_granticut_landau ->GetParError(7);
  Double_t constant2= fit_BBps_granticut_landau->GetParameter(8);
  Double_t mpv2 = fit_BBps_granticut_landau ->GetParameter(9);
  Double_t sigma2_landau= fit_BBps_granticut_landau ->GetParameter(10);
  Double_t constnat2_err= fit_BBps_granticut_landau ->GetParError(8);
  Double_t mpv2_err = fit_BBps_granticut_landau ->GetParError(9);
  Double_t sigma2_landau_err = fit_BBps_granticut_landau ->GetParError(10);
 
  TF1 *fit_result_poly_BBps_granticut_landau = new TF1("fit_result_poly_BBps_granticut_landau", poly4, 0.02,1.4,5);
  fit_result_poly_BBps_granticut ->SetParameters(par0_landau,par1_landau,par2_landau,par3_landau,par4_landau);
  TF1 *fit_result_Landau1_BBps_granticut = new TF1("fit_result_Landau1_BBps_granticut", Landau, 0.02,0.3,3);
  fit_result_Landau1_BBps_granticut ->SetParameters(constant,mpv, sigma_landau);
  TF1 *fit_result_Landau2_BBps_granticut = new TF1("fit_result_Landau2_BBps_granticut", Landau, 0.02,1.4,3);
  fit_result_Landau2_BBps_granticut ->SetParameters(constant2, mpv2, sigma2_landau);

  // fit_result_poly_BBps_granticut_landau ->SetLineColor(kMagenta - 6);
  // fit_result_Landau1_BBps_granticut ->SetLineColor(kViolet - 6);
  // fit_result_Landau2_BBps_granticut ->SetLineColor(kGreen -6);
  // fit_BBps_granticut_landau ->SetLineColor(kRed);

  // h_BBps_granticut_copy ->Draw("hist");
  // fit_result_poly_BBps_granticut_landau ->Draw("same");
  // fit_result_Landau1_BBps_granticut ->Draw("same");
  // fit_result_Landau2_BBps_granticut ->Draw("same");
  // fit_BBps_granticut_landau ->Draw("same");


  //  TF1 *fit_BBps_granticut_landau2 = new TF1("fit_BBps_granticut_landau2",Landau2, 0.02,1.4,6);
  //   fit_BBps_granticut_landau2 ->SetParameters( constant_init, mpv_init, sigma_landau_init, constant2_init, mpv2_init, sigma2_landau_init);
  //   TH1D *h_BBps_granticut_copy2 =  (TH1D*)(h_BBps_granticut->Clone("h_BBps_granticut_copy2"));
  // TFitResultPtr r3 = h_BBps_granticut_copy2 ->Fit(fit_BBps_granticut_landau2,"RQ0");
  // Int_t fitStatus3 = r3;
  // if (fitStatus3 !=0) {cout<<"fit error" <<endl;}

    TF1 *fit_BBps_granticut_landau_gaus = new TF1("fit_BBps_granticut_landau_gaus",poly4_Landau_Gaus, 0.02,1.4,11);
    fit_BBps_granticut_landau_gaus ->SetParameters( 0, 0, 0, 0,0, constant_init, mpv_init, sigma_landau_init, height2_init, mean2_init, sigma2_init);
    TH1D *h_BBps_granticut_copy3 =  (TH1D*)(h_BBps_granticut->Clone("h_BBps_granticut_copy3"));
  TFitResultPtr r3_1 = h_BBps_granticut->Fit(fit_BBps_granticut_landau_gaus,"RQ0");
  Int_t fitStatus3_1 = r3_1;
  if (fitStatus3_1 !=0) {cout<<"fit error poly4_Landau_Gaus" <<endl;}


  TF1 *fit_BBps_granticut_skewed_2gaus = new TF1("fit_BBps_granticut_skewed_2gaus" ,skewed_2gaus, 0.02, 1.4,8);
  fit_BBps_granticut_skewed_2gaus ->SetParameters(height_init, mean_init, sigma_init, alpha_init, height2_init, mean2_init, sigma2_init, alpha2_init);
  TFitResultPtr r3_2= h_BBps_granticut->Fit(fit_BBps_granticut_skewed_2gaus,"RQ0");
  Int_t fitStatus3_2 = r3_2;
  if (fitStatus3_2 !=0) {cout<<"fit error skewed 2 gaus" <<endl;}


  return;
}// END MAIN



// FIT FUNCTIONS
Double_t poly(Double_t *x, Double_t *par)
{
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) + par[5] * pow(x[0],5) +  par[6] * pow(x[0],6);
  return fit;
}

Double_t Gaus(Double_t *x, Double_t *par)
{
  Double_t height = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];

  Double_t fit = height * exp(-0.5 * pow((x[0] - mean)/sigma,2));
  return fit;
}

// Double_t GausTest(Double_t *x, Double_t *par){
//   Double_t mean = par[0];
//   Double_t sigma = par[1];
//   Double_t fit = TMath::Gaus(x[0], mean, sigma);
//   return fit;
// }

Double_t LandauTest(Double_t *x, Double_t *par)
{
  Double_t amp = par[0];
  Double_t mpv = par[1];
  Double_t sigma = par[2];
  // Double_t fit = amp * ::ROOT::Math::landau_pdf ( (x[0] - mpv)/sigma) ;
  Double_t fit = amp * ::TMath::Landau(x[0], mpv, sigma);
  return fit; 
}

Double_t Landau(Double_t *x, Double_t *par)
{
  Double_t amp = par[0];
  Double_t mpv = par[1];
  Double_t sigma = par[2];
  // Double_t fit = amp * ::ROOT::Math::landau_pdf ( (x[0] - mpv)/sigma) ;
  Double_t fit = amp * ::TMath::Landau(x[0], mpv, sigma);
  return fit; 
}

Double_t Landau2(Double_t *x, Double_t *par)
{
  Double_t constant = par[0];
  Double_t mpv = par[1];
  Double_t sigma = par[2];
   Double_t constant2 = par[3];
  Double_t mpv2 = par[4];
  Double_t sigma2= par[5];
  Double_t fit = constant* ::TMath::Landau(x[0],mpv,sigma) + constant2* ::TMath::Landau(x[0], mpv2, sigma2) ;
  return fit; 
}

Double_t poly4_2Landau(Double_t *x, Double_t *par)
{
  Double_t constant = par[5];
  Double_t mpv = par[6];
  Double_t sigma = par[7];
  Double_t constant2 = par[8];
  Double_t mpv2 = par[9];
  Double_t sigma2 = par[10];
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) + constant *::TMath::Landau(x[0],mpv,sigma) + constant2*::TMath::Landau(x[0],mpv2,sigma2) ;
  return fit;
}

Double_t poly4_Landau_Gaus(Double_t *x, Double_t *par)
{
  Double_t constant = par[5];
  Double_t mpv = par[6];
  Double_t sigma = par[7];
  Double_t height = par[8];
  Double_t mean = par[9];
  Double_t sigma_gaus= par[10];
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) +constant* ::TMath::Landau(x[0],mpv,sigma) + height * exp(-0.5 * pow((x[0] - mean)/sigma_gaus,2));
  return fit;
}

Double_t poly_Gaus(Double_t *x, Double_t *par)
{
   Double_t height = par[7];
   Double_t mean = par[8];
   Double_t sigma = par[9];

   Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) + par[5] * pow(x[0],5) +  par[6] * pow(x[0],6)+  height * exp(-0.5 * pow((x[0] - mean)/sigma,2)) ;
   return fit; 
}

Double_t poly_2Gaus(Double_t *x, Double_t *par) // need to keep to 11 params or under
{
   Double_t height = par[5];
   Double_t mean = par[6];
   Double_t sigma = par[7];

   Double_t height2 = par[8];
   Double_t mean2 = par[9];
   Double_t sigma2 = par[10];

   Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) +  height * exp(-0.5 * pow((x[0] - mean)/sigma,2)) + height2 * exp(-0.5 * pow((x[0] - mean2)/sigma2,2) );
   return fit; 
}

Double_t skewed_gaus(Double_t *x, Double_t *par)
{
  Double_t amp = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t alpha = par[3];
  Double_t fit = amp * exp(-0.5 * pow((x[0] - mean)/sigma,2)) * (1 + TMath::Erf( 0.7071 * alpha * (x[0] - mean)/sigma) );
  return fit;
}

Double_t skewed_2gaus(Double_t *x, Double_t *par)
{
  Double_t amp = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t alpha = par[3];
  Double_t amp2 = par[4];
  Double_t mean2 = par[5];
  Double_t sigma2 = par[6];
  Double_t alpha2 = par[7];
  Double_t fit = amp * exp(-0.5 * pow((x[0] - mean)/sigma,2)) * (1 + TMath::Erf( 0.7071 * alpha * (x[0] - mean)/sigma) ) + amp2 * exp(-0.5 * pow((x[0] - mean2)/sigma2,2)) * (1 + TMath::Erf( 0.7071 * alpha2 * (x[0] - mean2)/sigma2) );
  return fit;
}


Double_t poly_skewed_gaus(Double_t *x, Double_t *par)
{
  Double_t amp = par[7];
  Double_t mean = par[8];
  Double_t sigma = par[9];
  Double_t alpha = par[10];
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) + par[5] * pow(x[0],5) +  par[6] * pow(x[0],6) +  amp * exp(-0.5 * pow((x[0] - mean)/sigma,2)) * (1 + TMath::Erf( 0.7071 * alpha * (x[0] - mean)/sigma) );
  return fit;
}


Double_t poly5(Double_t *x, Double_t *par)
{
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) + par[5] * pow(x[0],5);
  return fit;
}

Double_t poly4(Double_t *x, Double_t *par)
{
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) ;
  return fit;
}

Double_t poly5_skewed_gaus(Double_t *x, Double_t *par)
{
  Double_t amp = par[6];
  Double_t mean = par[7];
  Double_t sigma = par[8];
  Double_t alpha = par[9];
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) + par[5] * pow(x[0],5) +  amp * exp(-0.5 * pow((x[0] - mean)/sigma,2)) * (1 + TMath::Erf( 0.7071 * alpha * (x[0] - mean)/sigma) );
  return fit;
}

Double_t poly4_skewed_gaus(Double_t *x, Double_t *par)
{
  Double_t amp = par[5];
  Double_t mean = par[6];
  Double_t sigma = par[7];
  Double_t alpha = par[8];
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) +  amp * exp(-0.5 * pow((x[0] - mean)/sigma,2)) * (1 + TMath::Erf( 0.7071 * alpha * (x[0] - mean)/sigma) );
  return fit;
}
