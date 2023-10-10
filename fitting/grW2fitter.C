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

void grW2fitter(){ //MAIN

  TFile *f1 = TFile::Open("../output/sbs8_mirror2.root"); // Load rootfile

  // Set up Canvas
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,1);

  TCanvas *c2 = new TCanvas();
  c2->Divide(2,1);

  // Load histograms
  TH1D *h_W2_gr_anticut= (TH1D*)f1->Get("h_W2_gr_anticut"); //put back to gr
  TH1D *h_W2_ps_anitcut = (TH1D*)f1->Get("h_W2_ps_anticut");
  TH1D *h_W2_gr_ps_anticut = (TH1D*)f1->Get("h_W2_gr_ps_anticut");
  TH1D *h_W2 = (TH1D*)f1->Get("h_W2");
  TH1D *h_W2_gr_cut = (TH1D*)f1->Get("h_W2_gr_cut");//put back to gr

  // Initial values for the fits
  Double_t par0_init = 50;
  Double_t par1_init = -200;
  Double_t par2_init = 200;
  Double_t par3_init = -200;
  Double_t par4_init = 200;
  Double_t par5_init = -200;
  Double_t par6_init = 200;
  Double_t height_init = 600;
  Double_t mean_init = 0.9;
  Double_t sigma_init = 0.1;
  Double_t amp_sk_init = 600;
  Double_t mean_sk_init = 0.78;
  Double_t sigma_sk_init = 0.1;
  Double_t alpha_init = 1;

  Double_t par0_cut_init = 50;
  Double_t par1_cut_init = -200;
  Double_t par2_cut_init = 200;
  Double_t par3_cut_init = -200;
  Double_t par4_cut_init = 200;
  Double_t par5_cut_init = -200;
  Double_t par6_cut_init = 200;
  Double_t height_cut_init = 3000;
  Double_t mean_cut_init = 0.9;
  Double_t sigma_cut_init = 0.1;
  Double_t amp_sk_cut_init = 3000;
  Double_t mean_sk_cut_init = 0.78;
  Double_t sigma_sk_cut_init = 0.1;
  Double_t alpha_cut_init = 1;

  Double_t par0, par1, par2, par3, par4,par5,par6, height, mean, sigma;

  Double_t width = h_W2_gr_anticut ->GetBinWidth(1);
  Int_t bin = 0.4/width;
  Double_t bincontent_W2_gr = h_W2_gr_anticut ->GetBinContent(bin);
  Double_t bincontent_W2_gr_ps = h_W2_gr_ps_anticut ->GetBinContent(bin);
  Double_t bincontent_W2_gr_cut = h_W2_gr_cut ->GetBinContent(bin);

  cout<<"bin content W2 gr: "<<bincontent_W2_gr<<endl;
  cout<<"bin content W2 gr ps: "<<bincontent_W2_gr_ps<<endl;
  cout<< "bin content W2 gr cut: "<<bincontent_W2_gr_cut <<endl;

  Double_t W2_scale_gr_ps = bincontent_W2_gr / bincontent_W2_gr_ps;
  Double_t W2_scale_gr_cut = bincontent_W2_gr_cut/ bincontent_W2_gr_ps;

  cout<<"scale for backgorund to gr anticut: "<< W2_scale_gr_ps <<endl;
  cout<<"scale for backgorund to gr cut: "<< W2_scale_gr_cut <<endl;

  TH1D *h_W2_gr_ps_anticut_scaled =  (TH1D*)(h_W2_gr_ps_anticut->Clone("h_W2_gr_ps_anticut_scaled"));
  h_W2_gr_ps_anticut_scaled ->Scale(W2_scale_gr_ps);
  
  TH1D *h_W2_gr_ps_cut_scaled =  (TH1D*)(h_W2_gr_ps_anticut->Clone("h_W2_gr_ps_cut_scaled"));
  h_W2_gr_ps_cut_scaled ->Scale(W2_scale_gr_cut);

  TH1D *h_W2_gr_anticut_subtracted = (TH1D*)(h_W2_gr_anticut->Clone("h_W2_gr_anticut_subtracted"));
  h_W2_gr_anticut_subtracted ->Add(h_W2_gr_ps_anticut_scaled, -1);
  

  TF1 *fit_background_anticut = new TF1("fit_background_anticut", poly, 0,1.4,7);
  fit_background_anticut ->SetParameters(par0_init,par1_init,par2_init,par3_init,par4_init, par5_init, par6_init);
  TFitResultPtr r = h_W2_gr_ps_anticut_scaled ->Fit(fit_background_anticut,"RQ0");
  Int_t fitStatus = r;
  if (fitStatus !=0) {cout<<"fit error" <<endl;}
  else{
    Double_t par0 = fit_background_anticut ->GetParameter(0);
    Double_t par1 = fit_background_anticut ->GetParameter(1);
    Double_t par2 = fit_background_anticut ->GetParameter(2);
    Double_t par3 = fit_background_anticut ->GetParameter(3);
    Double_t par4 = fit_background_anticut ->GetParameter(4);
    Double_t par5 = fit_background_anticut ->GetParameter(5);
    Double_t par6 = fit_background_anticut ->GetParameter(6);
    Double_t par0_err = fit_background_anticut ->GetParError(0);
    Double_t par1_err = fit_background_anticut ->GetParError(1);
    Double_t par2_err = fit_background_anticut ->GetParError(2);
    Double_t par3_err = fit_background_anticut ->GetParError(3);
    Double_t par4_err = fit_background_anticut ->GetParError(4);
    Double_t par5_err = fit_background_anticut ->GetParError(5);
    Double_t par6_err = fit_background_anticut ->GetParError(6);
    cout<< "par0: "<< par0 <<" +/- "<<par0_err<<endl;
    cout<< "par1: "<< par1<<" +/- "<<par1_err<<endl;
    cout<< "par2: "<< par2<<" +/- "<<par2_err<<endl;
    cout<< "par3: "<< par3<<" +/- "<<par3_err<<endl;
    cout<< "par4: "<< par4<<" +/- "<<par4_err<<endl;
    cout<< "par5: "<< par5<<" +/- "<<par5_err<<endl;
    cout<< "par6: "<< par6<<" +/- "<<par6_err<<endl;
  }

  TF1 *fit_W2_gr_subtraction =  new TF1("fit_W2_gr_subtraction", Gaus, 0.6,1.1,3);
  fit_W2_gr_subtraction->SetParameters(600,0.9,1);
  TFitResultPtr r2 = h_W2_gr_anticut_subtracted ->Fit(fit_W2_gr_subtraction,"RQ0");
  Int_t fitStatus2 = r2;
  if (fitStatus2 !=0) { cout<<"fit error" <<endl; }
  else{
    Double_t height= fit_W2_gr_subtraction ->GetParameter(0);
    Double_t mean = fit_W2_gr_subtraction ->GetParameter(1);
    Double_t sigma = fit_W2_gr_subtraction ->GetParameter(2);
    Double_t height_err= fit_W2_gr_subtraction ->GetParError(0);
    Double_t mean_err = fit_W2_gr_subtraction ->GetParError(1);
    Double_t sigma_err = fit_W2_gr_subtraction ->GetParError(2);
    cout<<"height: "<<height<<" +/- "<<height_err<<endl;
    cout<<"mean: "<<mean<<" +/- "<<mean_err<<endl;
    cout<<"sigma: "<<sigma<<" +/- "<<sigma_err<<endl;
  }

  TF1 *fit_W2_gr_subtraction_skew =  new TF1("fit_W2_gr_subtraction_skew", skewed_gaus, 0.6, 1.2, 4);
  fit_W2_gr_subtraction_skew->SetParameters(600,0.9,1,0);
  TFitResultPtr r2_1 = h_W2_gr_anticut_subtracted ->Fit(fit_W2_gr_subtraction_skew,"RQ0");
  Int_t fitStatus2_1 = r2_1;
  if (fitStatus2_1 !=0) { cout<<"fit error" <<endl; }
  else{
  }

  TF1 *fit_W2_gr_anticut_overall =  new TF1("fit_W2_gr_anticut_overall", poly_Gaus, 0, 1.4,10);
  fit_W2_gr_anticut_overall ->SetParameters(par0_init, par1_init, par2_init, par3_init, par4_init, par5_init, par6_init, height_init, mean_init, sigma_init);
  fit_W2_gr_anticut_overall ->SetParNames("p0", "p1", "p2", "p3", "p4", "p5", "p6","height","mean","sigma");
  //fit_W2_gr_anticut_overall ->SetParameters(par0, par1, par2, par3, par4, height, mean, sigma);
  TFitResultPtr r3 = h_W2_gr_anticut ->Fit(fit_W2_gr_anticut_overall,"RQ0");
  Int_t fitStatus3 = r3;
   if (fitStatus3 !=0) { cout<<"fit error" <<endl; }
     Double_t par0_overall = fit_W2_gr_anticut_overall ->GetParameter(0);
     Double_t par1_overall  = fit_W2_gr_anticut_overall ->GetParameter(1);
     Double_t par2_overall  = fit_W2_gr_anticut_overall ->GetParameter(2);
     Double_t par3_overall  = fit_W2_gr_anticut_overall ->GetParameter(3);
     Double_t par4_overall  = fit_W2_gr_anticut_overall ->GetParameter(4);
     Double_t par5_overall  = fit_W2_gr_anticut_overall ->GetParameter(5);
     Double_t par6_overall  = fit_W2_gr_anticut_overall ->GetParameter(6);
     Double_t par0_overall_err = fit_W2_gr_anticut_overall ->GetParError(0);
     Double_t par1_overall_err  = fit_W2_gr_anticut_overall ->GetParError(1);
     Double_t par2_overall_err  = fit_W2_gr_anticut_overall ->GetParError(2);
     Double_t par3_overall_err  = fit_W2_gr_anticut_overall ->GetParError(3);
     Double_t par4_overall_err  = fit_W2_gr_anticut_overall ->GetParError(4);
     Double_t par5_overall_err  = fit_W2_gr_anticut_overall ->GetParError(5);
     Double_t par6_overall_err  = fit_W2_gr_anticut_overall ->GetParError(6);
     cout<<endl;
     cout<<"overall fit"<<endl;
     cout<< "par0: "<< par0_overall <<" +/- "<<par0_overall_err<<endl;
     cout<< "par1: "<< par1_overall <<" +/- "<<par1_overall_err<<endl;
     cout<< "par2: "<< par2_overall <<" +/- "<<par2_overall_err <<endl;
     cout<< "par3: "<< par3_overall <<" +/- "<<par3_overall_err <<endl;
     cout<< "par4: "<< par4_overall <<" +/- "<<par4_overall_err<<endl;
     cout<< "par5: "<< par5_overall <<" +/- "<<par5_overall_err<<endl;
     cout<< "par6: "<< par6_overall <<" +/- "<<par6_overall_err<<endl;
     Double_t height_overall = fit_W2_gr_anticut_overall ->GetParameter(7);
     Double_t mean_overall  = fit_W2_gr_anticut_overall ->GetParameter(8);
     Double_t sigma_overall  =fit_W2_gr_anticut_overall ->GetParameter(9);
      Double_t height_overall_err = fit_W2_gr_anticut_overall ->GetParError(7);
     Double_t mean_overall_err  = fit_W2_gr_anticut_overall ->GetParError(8);
     Double_t sigma_overall_err  =fit_W2_gr_anticut_overall ->GetParError(9);
     cout<<"height: "<<height_overall <<" +/- "<<height_overall_err<<endl;
     cout<<"mean: "<<mean_overall <<" +/- "<<mean_overall_err<<endl;
     cout<<"sigma: "<<sigma_overall <<" +/- "<<sigma_overall_err<<endl;   

     TF1 *gaus_result =  new TF1("gaus_result", Gaus, 0,1.4,3);
     gaus_result ->SetParameters(height_overall, mean_overall, sigma_overall);
     gaus_result ->SetParNames("height", "mean","sigma");

     TF1 *poly_result = new TF1("poly_result",poly,0,1.4,7);
     poly_result ->SetParameters(par0_overall, par1_overall, par2_overall, par3_overall, par4_overall, par5_overall,par6_overall);
     
     c1->cd(2);
     
     h_W2_gr_anticut ->Draw("hist");
     gaus_result ->SetLineColor(kGreen -6);
     gaus_result ->SetFillColor(kGreen - 6);
     gaus_result ->SetFillStyle(3003);
     poly_result ->SetLineColor(kMagenta-6);
     poly_result ->SetFillColor(kMagenta-6);
     poly_result ->SetFillStyle(3003);
     fit_W2_gr_anticut_overall ->SetLineColor(kRed);
     fit_background_anticut ->SetLineColor(kViolet -6);
     fit_background_anticut ->SetLineStyle(2);

     poly_result ->Draw("same");
     gaus_result ->Draw("same");
     fit_W2_gr_anticut_overall ->Draw("same");
     // fit_background_anticut ->Draw("same");

     Double_t integral_gaus = gaus_result ->Integral(0.55,1.25); // need to account for bin width????
     //cout<<"integral of gaus: "<<integral_gaus <<endl; 
     Double_t counts_gaus =  integral_gaus/width;
     cout<< "counts gaus: "<< counts_gaus <<endl;

     Double_t integral_poly =  poly_result ->Integral(0.55,1.25);
     Double_t counts_poly = integral_poly/width;
     cout<< "counts poly: " <<counts_poly <<endl;
     //cout<< "sum: "<< integral_gaus/width + integral_poly/width <<endl;
 
     Double_t integral_overall = fit_W2_gr_anticut_overall ->Integral(0.55,1.25);
     //cout<< "integral of overall fit: "<< integral_overall<<endl;
     cout<<"counts overall: "<<integral_overall/width<<endl;

   

   TF1 *fit_W2_gr_anticut_overall_skew =  new TF1("fit_W2_gr_anticut_overall_skew", poly_skewed_gaus, 0, 1.4,11);
   fit_W2_gr_anticut_overall_skew ->SetParameters(par0_init, par1_init, par2_init, par3_init, par4_init, par5_init, par6_init, amp_sk_init, mean_sk_init, sigma_sk_init,alpha_init);
    TH1D *h_W2_gr_anticut_copy =  (TH1D*)(h_W2_gr_anticut->Clone("h_W2_gr_anticut_copy"));
    TFitResultPtr r3_1 = h_W2_gr_anticut_copy ->Fit(fit_W2_gr_anticut_overall_skew,"RQ0");
  Int_t fitStatus3_1 = r3_1;
   if (fitStatus3_1 !=0) { cout<<"fit error" <<endl; }
     Double_t par0_sk_overall = fit_W2_gr_anticut_overall_skew ->GetParameter(0);
     Double_t par1_sk_overall  = fit_W2_gr_anticut_overall_skew ->GetParameter(1);
     Double_t par2_sk_overall  = fit_W2_gr_anticut_overall_skew ->GetParameter(2);
     Double_t par3_sk_overall  = fit_W2_gr_anticut_overall_skew ->GetParameter(3);
     Double_t par4_sk_overall  = fit_W2_gr_anticut_overall_skew ->GetParameter(4);
     Double_t par5_sk_overall  = fit_W2_gr_anticut_overall_skew ->GetParameter(5);
     Double_t par6_sk_overall  = fit_W2_gr_anticut_overall_skew ->GetParameter(6);
     Double_t par0_sk_overall_err = fit_W2_gr_anticut_overall_skew ->GetParError(0);
     Double_t par1_sk_overall_err  = fit_W2_gr_anticut_overall_skew ->GetParError(1);
     Double_t par2_sk_overall_err  = fit_W2_gr_anticut_overall_skew->GetParError(2);
     Double_t par3_sk_overall_err  = fit_W2_gr_anticut_overall_skew ->GetParError(3);
     Double_t par4_sk_overall_err  = fit_W2_gr_anticut_overall_skew ->GetParError(4);
     Double_t par5_sk_overall_err  = fit_W2_gr_anticut_overall_skew ->GetParError(5);
     Double_t par6_sk_overall_err  = fit_W2_gr_anticut_overall_skew->GetParError(6);
     cout<<endl;
     cout<<"overall fit with skewed gaus"<<endl;
     cout<< "par0: "<< par0_sk_overall <<" +/- "<<par0_sk_overall_err<<endl;
     cout<< "par1: "<< par1_sk_overall <<" +/- "<<par1_sk_overall_err<<endl;
     cout<< "par2: "<< par2_sk_overall <<" +/- "<<par2_sk_overall_err <<endl;
     cout<< "par3: "<< par3_sk_overall <<" +/- "<<par3_sk_overall_err <<endl;
     cout<< "par4: "<< par4_sk_overall <<" +/- "<<par4_sk_overall_err<<endl;
     cout<< "par5: "<< par5_sk_overall <<" +/- "<<par5_sk_overall_err<<endl;
     cout<< "par6: "<< par6_sk_overall <<" +/- "<<par6_sk_overall_err<<endl;
     Double_t amp_sk_overall = fit_W2_gr_anticut_overall_skew ->GetParameter(7);
     Double_t mean_sk_overall  = fit_W2_gr_anticut_overall_skew  ->GetParameter(8);
     Double_t sigma_sk_overall  =fit_W2_gr_anticut_overall_skew  ->GetParameter(9);
     Double_t alpha_overall  =fit_W2_gr_anticut_overall_skew  ->GetParameter(10);
     Double_t amp_sk_overall_err = fit_W2_gr_anticut_overall_skew  ->GetParError(7);
     Double_t mean_sk_overall_err  = fit_W2_gr_anticut_overall_skew  ->GetParError(8);
     Double_t sigma_sk_overall_err  =fit_W2_gr_anticut_overall_skew  ->GetParError(9);
     Double_t alpha_overall_err  =fit_W2_gr_anticut_overall_skew  ->GetParError(10);
     cout<<"amp skew: "<<amp_sk_overall <<" +/- "<<amp_sk_overall_err<<endl;
     cout<<"mean skew: "<<mean_sk_overall <<" +/- "<<mean_sk_overall_err<<endl;
     cout<<"sigma skew: "<<sigma_sk_overall <<" +/- "<<sigma_sk_overall_err<<endl;  
     cout<<"alpha: "<<alpha_overall <<" +/- "<<alpha_overall_err<<endl;  

     TF1 *skewed_gaus_result =  new TF1("skewed_gaus_result", skewed_gaus, 0,1.4,4);
     skewed_gaus_result ->SetParameters(amp_sk_overall, mean_sk_overall, sigma_sk_overall,alpha_overall);
     //gaus_result ->SetParNames("height", "mean","sigma");

     TF1 *skewed_poly_result = new TF1("skewed_poly_result",poly,0,1.4,7);
     skewed_poly_result ->SetParameters(par0_sk_overall, par1_sk_overall, par2_sk_overall, par3_sk_overall, par4_sk_overall, par5_sk_overall,par6_sk_overall);

     c2->cd(2);
   
     skewed_gaus_result ->SetLineColor(kGreen -6);
     skewed_gaus_result ->SetFillColor(kGreen - 6);
     skewed_gaus_result ->SetFillStyle(3003);
     skewed_poly_result ->SetLineColor(kMagenta-6);
     skewed_poly_result ->SetFillColor(kMagenta-6);
     skewed_poly_result ->SetFillStyle(3003);
     fit_W2_gr_anticut_overall_skew ->SetLineColor(kRed);
     fit_background_anticut ->SetLineColor(kViolet -6);
     fit_background_anticut ->SetLineStyle(2);

     h_W2_gr_anticut_copy ->Draw("hist");
     skewed_poly_result ->Draw("same");
     skewed_gaus_result ->Draw("same");
     fit_W2_gr_anticut_overall_skew ->Draw("same");
     //fit_background_anticut ->Draw("same");
     
      Double_t integral_gaus_skew = skewed_gaus_result ->Integral(0.55,1.25); 
      Double_t counts_gaus_skew = integral_gaus_skew/width;
     cout<< "counts skewed gaus: "<< counts_gaus_skew <<endl;

     Double_t integral_poly_skew =  skewed_poly_result ->Integral(0.55,1.25);
     Double_t counts_poly_skew = integral_poly_skew/width;
     cout<< "counts poly skewed: " << counts_poly_skew<<endl;
    
     Double_t integral_overall_skew = fit_W2_gr_anticut_overall_skew ->Integral(0.55,1.25);
     cout<<"counts overall skewed: "<<integral_overall_skew/width<<endl;
  

   TF1 *fit_W2_gr_anticut_overall_skew5 =  new TF1("fit_W2_gr_anticut_overall_skew5", poly5_skewed_gaus, 0, 1.4,10);
   fit_W2_gr_anticut_overall_skew5 ->SetParameters(par0_init, par1_init, par2_init, par3_init, par4_init, par5_init, amp_sk_init, mean_sk_init, sigma_sk_init,alpha_init);
    TH1D *h_W2_gr_anticut_copy5 =  (TH1D*)(h_W2_gr_anticut->Clone("h_W2_gr_anticut_copy5"));
    TFitResultPtr r3_5 = h_W2_gr_anticut_copy5 ->Fit(fit_W2_gr_anticut_overall_skew5,"RQ0");
  Int_t fitStatus3_5 = r3_5;
   if (fitStatus3_5 !=0) { cout<<"fit error overall skewed poly5" <<endl; }
   else{
   }

   TF1 *fit_W2_gr_anticut_overall_skew4 =  new TF1("fit_W2_gr_anticut_overall_skew4", poly4_skewed_gaus, 0, 1.4,9);
   fit_W2_gr_anticut_overall_skew4 ->SetParameters(par0_init, par1_init, par2_init, par3_init, par4_init, amp_sk_init, mean_sk_init, sigma_sk_init,alpha_init);
    TH1D *h_W2_gr_anticut_copy4 =  (TH1D*)(h_W2_gr_anticut->Clone("h_W2_gr_anticut_copy4"));
    TFitResultPtr r3_4 = h_W2_gr_anticut_copy4 ->Fit(fit_W2_gr_anticut_overall_skew4,"RQ0");
  Int_t fitStatus3_4 = r3_4;
   if (fitStatus3_4 !=0) { cout<<"fit error overall skewed poly4" <<endl; }
   else{
   }


   TF1 *fit_background_anticut_big = new TF1("fit_background_anticut_big", poly, 0,1.4,7);
   fit_background_anticut_big ->SetParameters(par0_init,par1_init,par2_init,par3_init,par4_init, par5_init, par6_init);
   TFitResultPtr r_1 = h_W2_gr_ps_cut_scaled ->Fit(fit_background_anticut_big,"RQ0");
   Int_t fitStatus_1 = r_1;
   if (fitStatus_1 !=0) {cout<<"fit error" <<endl;}

   TH1D *h_W2_gr_cut_subtracted = (TH1D*)(h_W2_gr_cut->Clone("h_W2_gr_cut_subtracted"));
   h_W2_gr_cut_subtracted ->Add(h_W2_gr_ps_cut_scaled, -1);

   TF1 *fit_W2_gr_cut_subtracted =  new TF1("fit_W2_gr_cut_subtracted", Gaus, 0.6,1.1,3);
   fit_W2_gr_cut_subtracted->SetParameters(600,0.9,1);
   TFitResultPtr r2_2 = h_W2_gr_cut_subtracted ->Fit(fit_W2_gr_cut_subtracted,"RQ0");
   Int_t fitStatus2_2 = r2_2;
   if (fitStatus2_2 !=0) { cout<<"fit error" <<endl; }


   TF1 *fit_W2_gr_cut_overall =  new TF1("fit_W2_gr_cut_overall", poly_Gaus, 0, 1.4,10);
   fit_W2_gr_cut_overall ->SetParameters(par0_cut_init, par1_cut_init, par2_cut_init, par3_cut_init, par4_cut_init, par5_cut_init, par6_cut_init, height_cut_init, mean_cut_init, sigma_cut_init);
   fit_W2_gr_cut_overall ->SetParNames("p0", "p1", "p2", "p3", "p4", "p5", "p6","height","mean","sigma");
   TFitResultPtr r3_2 = h_W2_gr_cut ->Fit(fit_W2_gr_cut_overall,"RQ0");
   Int_t fitStatus3_2 = r3_2;
   if (fitStatus3_2 !=0) { cout<<"fit error" <<endl; }
   
     Double_t par0_cut_overall = fit_W2_gr_cut_overall ->GetParameter(0);
     Double_t par1_cut_overall  = fit_W2_gr_cut_overall ->GetParameter(1);
     Double_t par2_cut_overall  = fit_W2_gr_cut_overall ->GetParameter(2);
     Double_t par3_cut_overall  = fit_W2_gr_cut_overall ->GetParameter(3);
     Double_t par4_cut_overall  = fit_W2_gr_cut_overall ->GetParameter(4);
     Double_t par5_cut_overall  = fit_W2_gr_cut_overall ->GetParameter(5);
     Double_t par6_cut_overall  = fit_W2_gr_cut_overall ->GetParameter(6);
     Double_t par0_cut_overall_err = fit_W2_gr_cut_overall ->GetParError(0);
     Double_t par1_cut_overall_err  = fit_W2_gr_cut_overall ->GetParError(1);
     Double_t par2_cut_overall_err  = fit_W2_gr_cut_overall ->GetParError(2);
     Double_t par3_cut_overall_err  = fit_W2_gr_cut_overall ->GetParError(3);
     Double_t par4_cut_overall_err  = fit_W2_gr_cut_overall ->GetParError(4);
     Double_t par5_cut_overall_err  = fit_W2_gr_cut_overall ->GetParError(5);
     Double_t par6_cut_overall_err  = fit_W2_gr_cut_overall ->GetParError(6);
     cout<<endl;
     cout<<"overall fit for CUT"<<endl;
     cout<< "par0: "<< par0_cut_overall <<" +/- "<<par0_cut_overall_err<<endl;
     cout<< "par1: "<< par1_cut_overall <<" +/- "<<par1_cut_overall_err<<endl;
     cout<< "par2: "<< par2_cut_overall <<" +/- "<<par2_cut_overall_err <<endl;
     cout<< "par3: "<< par3_cut_overall <<" +/- "<<par3_cut_overall_err <<endl;
     cout<< "par4: "<< par4_cut_overall <<" +/- "<<par4_cut_overall_err<<endl;
     cout<< "par5: "<< par5_cut_overall <<" +/- "<<par5_cut_overall_err<<endl;
     cout<< "par6: "<< par6_cut_overall <<" +/- "<<par6_cut_overall_err<<endl;
     Double_t height_cut_overall = fit_W2_gr_cut_overall ->GetParameter(7);
     Double_t mean_cut_overall  = fit_W2_gr_cut_overall ->GetParameter(8);
     Double_t sigma_cut_overall  =fit_W2_gr_cut_overall ->GetParameter(9);
      Double_t height_cut_overall_err = fit_W2_gr_cut_overall ->GetParError(7);
     Double_t mean_cut_overall_err  = fit_W2_gr_cut_overall ->GetParError(8);
     Double_t sigma_cut_overall_err  =fit_W2_gr_cut_overall ->GetParError(9);
     cout<<"height: "<<height_cut_overall <<" +/- "<<height_cut_overall_err<<endl;
     cout<<"mean: "<<mean_cut_overall <<" +/- "<<mean_cut_overall_err<<endl;
     cout<<"sigma: "<<sigma_cut_overall <<" +/- "<<sigma_cut_overall_err<<endl;   

     TF1 *gaus_result_cut =  new TF1("gaus_result_cut", Gaus, 0,1.4,3);
     gaus_result_cut ->SetParameters(height_cut_overall, mean_cut_overall, sigma_cut_overall);
     gaus_result_cut ->SetParNames("height", "mean","sigma");

     TF1 *poly_result_cut = new TF1("poly_result_cut",poly,0,1.4,7);
     poly_result_cut ->SetParameters(par0_cut_overall, par1_cut_overall, par2_cut_overall, par3_cut_overall, par4_cut_overall, par5_cut_overall,par6_cut_overall);

     c1->cd(1);
     
     h_W2_gr_cut ->Draw("hist");
     gaus_result_cut ->SetLineColor(kGreen -6);
     gaus_result_cut ->SetFillColor(kGreen - 6);
     gaus_result_cut ->SetFillStyle(3003);
     poly_result_cut ->SetLineColor(kMagenta-6);
     poly_result_cut ->SetFillColor(kMagenta-6);
     poly_result_cut ->SetFillStyle(3003);
     fit_W2_gr_cut_overall ->SetLineColor(kRed);
     fit_background_anticut_big ->SetLineColor(kViolet -6);
     fit_background_anticut_big ->SetLineStyle(2);

     poly_result_cut ->Draw("same");
     gaus_result_cut ->Draw("same");
     fit_W2_gr_cut_overall ->Draw("same");
     // fit_background_anticut_big ->Draw("same");

     Double_t integral_gaus_cut = gaus_result_cut ->Integral(0.55,1.25);
     Double_t counts_gaus_cut = integral_gaus_cut/width;
     cout<< "counts gaus for cut: "<< counts_gaus_cut <<endl;

     Double_t integral_poly_cut =  poly_result_cut ->Integral(0.55,1.25);
     Double_t counts_poly_cut = integral_poly_cut/width; 
     cout<< "counts poly: " << counts_poly_cut <<endl;
     //cout<< "sum: "<< integral_gaus/width + integral_poly/width <<endl;
 
     Double_t integral_overall_cut = fit_W2_gr_cut_overall ->Integral(0.55,1.25);
     //cout<< "integral of overall fit: "<< integral_overall<<endl;
     cout<<"counts overall: "<<integral_overall_cut/width<<endl;
   
   

   TF1 *fit_W2_gr_cut_overall_skew =  new TF1("fit_W2_gr_cut_overall_skew", poly_skewed_gaus, 0, 1.4,11);
   fit_W2_gr_cut_overall_skew ->SetParameters(par0_cut_init, par1_cut_init, par2_cut_init, par3_cut_init, par4_cut_init, par5_cut_init, par6_cut_init, amp_sk_cut_init, mean_sk_cut_init, sigma_sk_cut_init,alpha_cut_init);
    TH1D *h_W2_gr_cut_copy =  (TH1D*)(h_W2_gr_cut->Clone("h_W2_gr_cut_copy"));
    TFitResultPtr r3_3 = h_W2_gr_cut_copy ->Fit(fit_W2_gr_cut_overall_skew,"RQ0");
  Int_t fitStatus3_3 = r3_3;
   if (fitStatus3_3 !=0) { cout<<"fit error" <<endl; }
   
     Double_t par0_sk_cut_overall = fit_W2_gr_cut_overall_skew ->GetParameter(0);
     Double_t par1_sk_cut_overall  = fit_W2_gr_cut_overall_skew ->GetParameter(1);
     Double_t par2_sk_cut_overall  = fit_W2_gr_cut_overall_skew ->GetParameter(2);
     Double_t par3_sk_cut_overall  = fit_W2_gr_cut_overall_skew ->GetParameter(3);
     Double_t par4_sk_cut_overall  = fit_W2_gr_cut_overall_skew ->GetParameter(4);
     Double_t par5_sk_cut_overall  = fit_W2_gr_cut_overall_skew ->GetParameter(5);
     Double_t par6_sk_cut_overall  = fit_W2_gr_cut_overall_skew ->GetParameter(6);
     Double_t par0_sk_cut_overall_err = fit_W2_gr_cut_overall_skew ->GetParError(0);
     Double_t par1_sk_cut_overall_err  = fit_W2_gr_cut_overall_skew ->GetParError(1);
     Double_t par2_sk_cut_overall_err  = fit_W2_gr_cut_overall_skew->GetParError(2);
     Double_t par3_sk_cut_overall_err  = fit_W2_gr_cut_overall_skew ->GetParError(3);
     Double_t par4_sk_cut_overall_err  = fit_W2_gr_cut_overall_skew ->GetParError(4);
     Double_t par5_sk_cut_overall_err  = fit_W2_gr_cut_overall_skew ->GetParError(5);
     Double_t par6_sk_cut_overall_err  = fit_W2_gr_cut_overall_skew->GetParError(6);
     cout<<endl;
     cout<<"overall fit with skewed gaus for CUT"<<endl;
     cout<< "par0: "<< par0_sk_cut_overall <<" +/- "<<par0_sk_cut_overall_err<<endl;
     cout<< "par1: "<< par1_sk_cut_overall <<" +/- "<<par1_sk_cut_overall_err<<endl;
     cout<< "par2: "<< par2_sk_cut_overall <<" +/- "<<par2_sk_cut_overall_err <<endl;
     cout<< "par3: "<< par3_sk_cut_overall <<" +/- "<<par3_sk_cut_overall_err <<endl;
     cout<< "par4: "<< par4_sk_cut_overall <<" +/- "<<par4_sk_cut_overall_err<<endl;
     cout<< "par5: "<< par5_sk_cut_overall <<" +/- "<<par5_sk_cut_overall_err<<endl;
     cout<< "par6: "<< par6_sk_cut_overall <<" +/- "<<par6_sk_cut_overall_err<<endl;
     Double_t amp_sk_cut_overall = fit_W2_gr_cut_overall_skew ->GetParameter(7);
     Double_t mean_sk_cut_overall  = fit_W2_gr_cut_overall_skew  ->GetParameter(8);
     Double_t sigma_sk_cut_overall  =fit_W2_gr_cut_overall_skew  ->GetParameter(9);
     Double_t alpha_cut_overall  =fit_W2_gr_cut_overall_skew  ->GetParameter(10);
     Double_t amp_sk_cut_overall_err = fit_W2_gr_cut_overall_skew  ->GetParError(7);
     Double_t mean_sk_cut_overall_err  = fit_W2_gr_cut_overall_skew  ->GetParError(8);
     Double_t sigma_sk_cut_overall_err  =fit_W2_gr_cut_overall_skew  ->GetParError(9);
     Double_t alpha_cut_overall_err  =fit_W2_gr_cut_overall_skew  ->GetParError(10);
     cout<<"amp skew cut: "<<amp_sk_cut_overall <<" +/- "<<amp_sk_cut_overall_err<<endl;
     cout<<"mean skew: "<<mean_sk_cut_overall <<" +/- "<<mean_sk_cut_overall_err<<endl;
     cout<<"sigma skew: "<<sigma_sk_cut_overall <<" +/- "<<sigma_sk_cut_overall_err<<endl;  
     cout<<"alpha: "<<alpha_cut_overall <<" +/- "<<alpha_cut_overall_err<<endl;  

     TF1 *skewed_gaus_result_cut =  new TF1("skewed_gaus_result_cut", skewed_gaus, 0,1.4,4);
     skewed_gaus_result_cut ->SetParameters(amp_sk_cut_overall, mean_sk_cut_overall, sigma_sk_cut_overall,alpha_cut_overall);
     //gaus_result ->SetParNames("height", "mean","sigma");

     TF1 *skewed_poly_result_cut = new TF1("skewed_poly_result_cut",poly,0,1.4,7);
     skewed_poly_result_cut ->SetParameters(par0_sk_cut_overall, par1_sk_cut_overall, par2_sk_cut_overall, par3_sk_cut_overall, par4_sk_cut_overall, par5_sk_cut_overall,par6_sk_cut_overall);

     c2->cd(1);

     skewed_gaus_result_cut ->SetLineColor(kGreen -6);
     skewed_gaus_result_cut ->SetFillColor(kGreen - 6);
     skewed_gaus_result_cut ->SetFillStyle(3003);
     skewed_poly_result_cut ->SetLineColor(kMagenta-6);
     skewed_poly_result_cut ->SetFillColor(kMagenta-6);
     skewed_poly_result_cut ->SetFillStyle(3003);
     fit_W2_gr_cut_overall_skew ->SetLineColor(kRed);
     fit_background_anticut_big ->SetLineColor(kViolet -6);
     fit_background_anticut_big ->SetLineStyle(2);

     h_W2_gr_cut_copy ->Draw("hist");
     skewed_poly_result_cut ->Draw("same");
     skewed_gaus_result_cut ->Draw("same");
     fit_W2_gr_cut_overall_skew ->Draw("same");
     // fit_background_anticut_big ->Draw("same");
     
      Double_t integral_gaus_skew_cut = skewed_gaus_result_cut ->Integral(0.55,1.25); 
      Double_t counts_gaus_skew_cut =  integral_gaus_skew_cut/width;
     cout<< "counts skewed gaus cut: "<< counts_gaus_skew_cut <<endl;

     Double_t integral_poly_skew_cut =  skewed_poly_result_cut ->Integral(0.55,1.25);
     Double_t counts_poly_skew_cut = integral_poly_skew_cut/width;
     cout<< "counts poly skewed: " <<  counts_poly_skew_cut <<endl;
    
     Double_t integral_overall_skew_cut = fit_W2_gr_cut_overall_skew ->Integral(0.55,1.25);
     cout<<"counts overall skewed: "<<integral_overall_skew_cut/width<<endl;
  
   
     cout<<endl;
     Double_t gaus_eff = counts_gaus_cut / (counts_gaus_cut + counts_gaus);
     cout<<"Eff from normal gaus::  " << counts_gaus_cut <<" / ("<< counts_gaus_cut<< " + "<< counts_gaus<<") = " <<   gaus_eff<<endl;
     Double_t skew_eff = counts_gaus_skew_cut / (counts_gaus_skew_cut + counts_gaus_skew);
     cout<<"Eff from skewed gaus::  " << counts_gaus_skew_cut <<" / ("<< counts_gaus_skew_cut<< " + "<< counts_gaus_skew<<") = " <<  skew_eff<<endl;

     cout<<endl;
     Double_t gaus_iden_eff = (counts_gaus_cut + counts_poly) / (counts_gaus_cut + counts_poly +counts_poly_cut + counts_gaus);
     cout<<"Identification Eff gaus:: "<< "("<<counts_gaus_cut <<" + " << counts_poly <<") / ("<< counts_gaus_cut <<" + " <<counts_poly << " + " <<counts_poly_cut <<" + " << counts_gaus<<")= "<< gaus_iden_eff<<endl;

     cout<<endl;
     Double_t gaus_ineff = counts_poly_cut/ (counts_poly_cut + counts_poly);
     cout<<"fasle-positives evaluation from normal gaus::  " << counts_poly_cut<<" / ("<< counts_poly_cut<< " + "<< counts_poly<<") = " <<   gaus_ineff<<endl;
     Double_t skew_ineff = counts_poly_skew_cut / (counts_poly_skew_cut + counts_poly_skew);
     cout<<"false-positives evaluation from skewed gaus::  " << counts_poly_skew_cut <<" / ("<< counts_poly_skew_cut<< " + "<< counts_poly_skew<<") = " <<  skew_ineff<<endl;
    

     

  
}// end MAIN


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


Double_t poly_Gaus(Double_t *x, Double_t *par)
{
   Double_t height = par[7];
   Double_t mean = par[8];
   Double_t sigma = par[9];

   Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) + par[5] * pow(x[0],5) +  par[6] * pow(x[0],6)+  height * exp(-0.5 * pow((x[0] - mean)/sigma,2));
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
