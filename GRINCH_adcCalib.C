#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
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
#include<math.h>
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>

#include <iterator>
#include <istream>
#include <vector>
#include <string>


using namespace std;


TCanvas *c1;
TCanvas *allC;
//TCanvas *allC2;
TCanvas *allC3;
TCanvas *allCm;//maria
const Int_t nADC = 64;

TH1D *hvadc0[nADC];
TH1D *adc[nADC];
TH1D *hvadc3[nADC];
TH1D *adc_test[nADC];//maria 


TH1D *ped_w;
TH1D *spe;
TH1D *spe2;

string ped_chan[64];
string spe_chan[64];
string spe_chan_sub[64];

stringstream pc;
stringstream sc;
stringstream scs;

std::vector<int> v_adc[nADC];
Double_t ped[nADC];

TF1 *ped_fit;
TF1 *adc_fit;
TF1 *ped_adc_fit;

int i;

void style(int n)
{
  if(n == 0)
    {
      gStyle->SetOptStat(0);
    }
  else
    {
      gStyle->SetOptStat(1111111);
      gStyle->SetOptFit(1);
    }

}

void all_style()
{
      gStyle ->SetTitleSize(.15, " ");

}

//void GRINCH_adcCalib(TString basename="e1209019_trigtest_10460.evio.0") {
//void GRINCH_adcCalib(TString basename="e1209019_11186.evio.0") {
void GRINCH_adcCalib(TString basename="grinch_") {

  //RUNNUM
  // ******************
  Int_t runnum = 12338;
  // ******************


   TString inputroot;
   inputroot="../rootfiles/"+basename+runnum+".evio.0.root";
   // inputroot = "../../../Rootfiles/grinchled_13108_fev_0_nev_-1_fseg_0_lseg_4.root";


TFile *file = new TFile(inputroot); 
TTree *eventtree = (TTree*) file->Get("T");
 Int_t ntdc;
 Double_t grinch_tdc_pmt[512];
 Double_t grinch_tdc_le[512];
 Double_t grinch_tdc_te[512];
 Double_t grinch_tdc_mult[512];
 Double_t grinch_tdc_tot[512];

 Double_t grinch_adc_pmt[nADC];
 Double_t grinch_adc[nADC];
 Double_t adcGAtime[1000];

// Fill the arrays from ROOT branches (variables to use)
 eventtree->SetBranchAddress("Ndata.bb.grinch_tdc.tdc",&ntdc) ; // Number of TDC (LE) in event
 eventtree->SetBranchAddress("bb.grinch_tdc.tdcelemID",&grinch_tdc_pmt) ; // PMT number
 eventtree->SetBranchAddress("bb.grinch_tdc.tdc",&grinch_tdc_le) ; // TDC LE
 eventtree->SetBranchAddress("bb.grinch_tdc.tdc_te",&grinch_tdc_te) ; // TDC TE
 eventtree->SetBranchAddress("bb.grinch_tdc.tdc_mult",&grinch_tdc_mult) ; // TDC multiplicity
 eventtree->SetBranchAddress("bb.grinch_adc.a_time",&adcGAtime) ;
// eventtree->SetBranchAddress("bb.grinch_tdc.tot",&grinch_tdc_tot) ; // TDC time over threshold

 eventtree->SetBranchAddress("bb.grinch_adc.adcelemID",&grinch_adc_pmt) ; // ADC number
 eventtree->SetBranchAddress("bb.grinch_adc.a_p",&grinch_adc) ; // ADC counts 
 //changed above from .a_p // I set it back for a minute -B.Y.

// Histograms
 TH1F *GRINCH_ntdc = new TH1F("GRINCH_ntdc"," ; Nhit ; Counts",512,0,511); // 
 TH2F *GRINCH_tdc_LE_vs_pmt = new TH2F("GRINCH_tdc_LE_vs_pmt"," ; PMT # ; TDC leading edge (ns)", 512,0,511,100,1,100);
 TH2F *GRINCH_tdc_TE_vs_pmt = new TH2F("GRINCH_tdc_TE_vs_pmt"," ; PMT # ; TDC trailing edge (ns)", 512,0,511,100,1,100);
 TH2F *GRINCH_mult_vs_pmt = new TH2F("GRINCH_mult_vs_pmt"," ; PMT # ; PMT hit multiplicity", 512,0,511,20,1,20);
 TH2F *GRINCH_tdc_ToT_vs_pmt = new TH2F("GRINCH_tdc_ToT_vs_pmt"," ; PMT # ; Time over threshold (ns)", 512,0,511,200,1,200);

 TH2F *GRINCH_adc_chan_vs_num = new TH2F("GRINCH_adc_chan_vs_num "," ; ADC # ; ADC channel",64,0,63,311, -10, 300);

 TH1D *projY = new TH1D("projY", "projY",311,-10,300);

// ADC histograms to fit
for (int k = 0; k< nADC; k++)
{
   hvadc0[k] = new TH1D(Form("hvadc0[%d]",k),Form("hvadc0[%d]",k), 201, -0.5, 200.5);
   adc[k] = new TH1D(Form("adc[%d] %i",k,runnum),Form("adc[%d] %i",k,runnum), 201, -0.5, 200.5);
   //adc_test[k] = new TH1D(Form("adc_test[%d] %i",k,runnum),Form("adc_test[%d] %i",k,runnum), 201, -0.5, 200.5);
   hvadc3[k] = new TH1D(Form("hvadc3[%d]",k),Form("hvadc3[%d]",k), 201, -0.5, 200.5);
}

 Long64_t nentries = eventtree->GetEntries();
	for (int i = 0; i < nentries; i++) {
//      for (int i = 0; i < 1000; i++) { // for debugging
      		eventtree->GetEntry(i);		
                if (i%1000==0) cout << " Entry = " << i << endl;

//		GRINCH_ntdc->Fill(float(ntdc));

//                // Loop over TDC LE
//		for (Int_t n=0;n<512;n++) {
//		  GRINCH_tdc_LE_vs_pmt->Fill( grinch_tdc_pmt[n],grinch_tdc_le[n] );
//                  GRINCH_tdc_TE_vs_pmt->Fill( grinch_tdc_pmt[n],grinch_tdc_te[n] );
//                  GRINCH_mult_vs_pmt->Fill( grinch_tdc_pmt[n],grinch_tdc_mult[n] );
//                  GRINCH_tdc_ToT_vs_pmt->Fill( grinch_tdc_pmt[n],(grinch_tdc_te[n]-grinch_tdc_le[n]) );
                  //h_ToT_vs_pmt->Fill( grinch_tdc_pmt[n],grinch_tdc_tot[n] ); // The SBS-Offline variable isn't filled correctly...
//		}

                // Fill the "ADC channels vs. ADC #" array//probably wrong?
                for (Int_t n=0;n<nADC;n++) {
                   
                   GRINCH_adc_chan_vs_num->Fill( grinch_adc_pmt[n],grinch_adc[n] );

                   //hvadc0[n]->Fill(grinch_adc[n]);
                   //adc[n]->Fill(grinch_adc[n]);
                   //hvadc3[n]->Fill(grinch_adc[n]);
		   //adc_test[n]->Fill(grinch_adc[n]);
		   
                   hvadc0[n]=GRINCH_adc_chan_vs_num->ProjectionY(Form("hvadc0[%d]",n+1),n+1,n+1,"");
                   adc[n]=GRINCH_adc_chan_vs_num->ProjectionY(Form("adc[%d]",n+1),n+1,n+1,"");
                   hvadc3[n]=GRINCH_adc_chan_vs_num->ProjectionY(Form("hvadc3[%d]",n+1),n+1,n+1,"");
		   //adc_test[n]=GRINCH_adc_chan_vs_num->ProjectionY(Form("adc_test[%d]",n+1),n+1,n+1,"");

                }

	}  // event loop

	//projY = GRINCH_adc_chan_vs_num -> ProjectionY("py",1,1);
  // Draws the ADC vs. channel (SBS-Offline output)
  c1 = new TCanvas("c1", "ADC Channels vs. ADC #",10,10,1000,750);
  c1->cd()->SetLogz();
  gStyle->SetOptStat(11);
  gStyle->SetStatX(0.60);
  //c1->Modified(); c1->Update(); gSystem->ProcessEvents();
  GRINCH_adc_chan_vs_num->Draw("COLZ");

  // Draws the fit ADC peaks, and spe vs. channel
  allC =  new TCanvas("allC","All ADC Fits",800,800);
  allC -> SetCanvasSize(1800, 2000);
  allC -> SetWindowSize(1800, 800);

  allC -> Divide(8,8);

  // Maria's. Draws the fit ADC peaks, and spe vs. channel
  //allCm =  new TCanvas("allCm","All ADC Test",800,800);
  //allCm -> SetCanvasSize(1800, 2000);
  //allCm -> SetWindowSize(1800, 800);

  //allCm -> Divide(8,8);

//  allC2 =  new TCanvas("allC2","Pedestal Width vs. ADC Channel",800,800);
//  allC2 -> SetCanvasSize(1000, 400);
//  allC2 -> SetWindowSize(1000, 450);
//  allC2 -> Divide(8,8);

  allC3 =  new TCanvas("allC3","SPE Channel vs. ADC #",800,800);
  allC3 -> SetCanvasSize(1000, 400);
  allC3 -> SetWindowSize(1000, 450);

  allC3 -> Divide(8,8);

  style(1);

//  for (int i = 0; i< nADC; i++)
//    {
      //      delete adc[i];
      //adc[i] = new TH1D(Form("vADC%d",i),Form("vADC%d",i), 101, -0.5, 4000.5);
//      adc[i] = new TH1D(Form("ADC[%d]",i),Form("ADC[%d]",i), 501, -0.5, 500.5);
//      hvadc2[i] = new TH1D(Form("vADC%d",i),Form("vADC%d",i), 501, -0.5, 500.5); // fewer bins for lower statistics
//      hvadc3[i] = new TH1D(Form("vADC[%d]",i),Form("vADC[%d]",i), 501, -0.5, 500.5);
//      for (int j = 0; j < grinch_adc[i].size(); j++)
//      {
//         adc[i]->Fill(grinch_adc[i]);
//         hvadc2[i]->Fill(grinch_adc[i]);
//         hvadc3[i]->Fill(grinch_adc[i]);
//      }
//    }

  ped_w = new TH1D("ped_w","Pedestal mean vs. ADC #", 64, 0, 63);
  spe = new TH1D("spe","SPE channels (pedestal subtracted)", 64, 0, 63);
  spe2 = new TH1D("spe2","SPE channel vs. ADC #", 64, 0, 63);

  cout<<"Histo #: "<<i<<endl;

  //maria
  for (int i=0; i<nADC;i++)
    {
      //allCm->cd(i+1)->SetLogy();//maria
      //all_style();
      //adc_test[i]->GetXaxis()->SetRange(0,300);
      //adc_test[i]->Draw("same"); // Commented out because the original is fixed now -B.Y.

    }

  for (int i = 0; i< nADC; i++)
    {
     
      allC->cd(i+1)->SetLogy();
      all_style();
      //      cout<< adc[i]->GetMean()<<endl;

      //BG_sig = new TF1("ped_sig","[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2])) + x*log([3])-[3]-TMath::LnGamma(x+1)");
      //BG_sig -> SetParNames("ped_const","ped_mean","ped_sigma","","","");

      //ped = new TF1("ped","gaus");

      //adc[i]->Fit("gaus","","",0,150); // Fit the ADC pedestal

      ped_fit = new TF1("ped_fit","[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2]))");
      double peak_pos = hvadc0[i]->GetMaximumBin();
      //double peak_pos = grinch_adc[i].max();

      ped_fit->SetParameters(hvadc0[i]->GetBinContent(peak_pos),peak_pos,2);
      //ped_fit->SetParameter(1,peak_pos);

      //ped_adc_fit = new TF1("ped_adc_fit","[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2])) + [3]*(x-[4])*TMath::Exp(-(x-[4])/[5])");

      //ped_fit->SetLineColor(2);
      hvadc0[i]->Fit("ped_fit","","",peak_pos-20,peak_pos+20);
      //adc[i]->Fit("ped_fit");

      // "end of pedestal", defined to be "x-sigma from pedestal" (adjustable)
      double ped_end = ped_fit->GetParameter(1)+5*ped_fit->GetParameter(2);

//      adc_fit = new TF1("adc_fit","[3]*(x-[0])*TMath::Exp(-(x-[1])/[2])");
      //adc_fit->SetParameter(3,ped_fit->GetParameter(0)/10); // scale spe height param by pedestal height (adjustable)
//      adc_fit->SetParameter(0,ped_end);

      adc_fit = new TF1("adc_fit","[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2]))");
      //adc_fit->SetParameter(3,ped_fit->GetParameter(0)/10); // scale spe height param by pedestal height (adjustable)
      adc_fit->SetParameter(1,50);
      //adc_fit->SetParameter(2,4);



      // TRY THIS
      //ped_adc_fit = new TF1("ped_adc_fit","[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2])) + [3]*(x-[4])*TMath::Exp(-(x-[4])/[5])");


      //cout<< hvadc0[i]->GetBinContent(peak_pos)<<endl;
      


      //adc[i] -> Draw();

      //allC2->cd(i+1)->SetLogy();
      //double ped_end = ped_fit->GetParameter(1)+5*ped_fit->GetParameter(2);
      //hvadc0[i]->Clone(Form("adc[%d]",i));
      //adc[i]->GetXaxis()->SetRange(ped_end,500);
      //double spe_mean = adc[i]->GetMean();

      adc[i]->GetXaxis()->SetRange(ped_end,300);
//      adc_fit->SetParameter(1,adc[i]->GetMean());
//      adc_fit->SetParameter(2,adc[i]->GetRMS());
      adc_fit->SetLineColor(3);
      adc[i]->GetXaxis()->SetRange(0,300);

      //cout<< adc[i]->GetRMS()<<endl;

      //adc[i]->Fit("adc_fit","","",ped_end,300); // Fit the spe spectrum
//      adc[i]->Fit("adc_fit","","",ped_end,20);
      adc[i]->Fit("adc_fit","","",ped_end,100);
      //hvadc0[i] -> Draw();
      adc[i] -> Draw("same");


      ped_adc_fit = new TF1("ped_adc_fit","[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2])) + (x-[3])*TMath::Exp(-(x-[4])/[5])");
      ped_adc_fit -> FixParameter(0,ped_fit->GetParameter(0));
      ped_adc_fit -> FixParameter(1,ped_fit->GetParameter(1));
      ped_adc_fit -> FixParameter(2,ped_fit->GetParameter(2));
      ped_adc_fit -> FixParameter(3,adc_fit->GetParameter(0));
      ped_adc_fit -> SetParameter(4,adc_fit->GetParameter(1));
      ped_adc_fit -> SetParameter(5,adc_fit->GetParameter(2));
      //ped_adc_fit -> SetParameter(6,adc_fit->GetParameter(3));

      ped_adc_fit -> SetParNames("ped_const","ped_mean","ped_sigma","start","sig_MPV","sig_width");
      ped_adc_fit -> SetLineColor(94);

      hvadc3[i]->Fit("ped_adc_fit","","",peak_pos,300);

      //adc = new TF1("adc","[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2]))");
      //adc -> FixParameter(0,ped_adc_fit->GetParameter(3));
      //adc -> FixParameter(1,ped_adc_fit->GetParameter(4));
      //adc -> FixParameter(2,ped_adc_fit->GetParameter(5));

      //spe->Draw();
      adc[i] -> Draw();
      hvadc0[i] -> Draw("same");
  //    hvadc3[i] -> Draw("same");
      //adc -> SetLineColor(5);
      //adc_fit -> SetLineStyle(5);
      //adc -> Draw("same");

      ped_w->SetBinContent(i,ped_fit->GetParameter(1));
      //spe->SetBinContent(i,adc_fit->GetParameter(1));
      //ped_w->SetBinContent(i,ped_adc_fit->GetParameter(2));
      //spe->SetBinContent(i,ped_adc_fit->GetParameter(4));
      //adc[i]->SetRangeUser(ped_end,500);


      double peddouble = ped_fit->GetParameter(1);
      //double peddouble = peak_pos;
//      double spedouble = adc_fit->GetMaximumX(ped_end,50); // TODO: Currently uses "maximum X". Use fits instead.
      double spedouble = adc_fit->GetParameter(1);

      int pedint = round(peddouble);
      if( pedint<0 ){ pedint = 0; }
      int speint = round(spedouble);
      if( speint<0 ){ speint = 0; }
      int spepedint = speint-pedint;
      if(spepedint<0){ spepedint = 0; }

      //spe->SetBinContent(i,adc_fit->GetMaximumX(ped_end,500)-ped_fit->GetParameter(1));
      spe->SetBinContent(i,speint);
      //spe->SetBinContent(i,adc_fit->GetParameter(1));

      cout<<""<<endl;
      cout<<"****************************"<<endl;
      cout<<"ADC #"<<i<<endl;
      cout<<""<<endl;
      cout<<"pedestal mean channel: "<<pedint<<endl;
      //cout<<"ped_height (max content): "<<peak_pos<<endl;
      cout<<"pedestal end (+3 sigma): "<<ped_end<<endl;
      cout<<""<<endl;
      cout<<"spe channel: "<<speint<<endl;
      cout<<"counts at spe peak: "<<adc_fit->GetMaximum(ped_end,200)<<endl;
      cout<<""<<endl;
      cout<<"****************************"<<endl;
      cout<<""<<endl;

      // Make entries for .txt file
      pc.str("");
      pc.clear();
      pc << pedint;
      ped_chan[i] = pc.str();

      sc.str("");
      sc.clear();
      sc << speint;
      spe_chan[i] = sc.str();

      scs.str("");
      scs.clear();
      scs << spepedint;
      spe_chan_sub[i] = scs.str();

    }


  //allC2->cd()->SetLogy(0);
  //style(0);

  ped_w->SetXTitle("ADC #");
  ped_w->SetYTitle("pedestal mean");
  //ped_w->Draw();

  allC3->cd()->SetLogy(0);
  style(0);
  spe->GetYaxis()->SetRangeUser(0,200);
  spe->SetXTitle("ADC #");
  spe->SetYTitle("spe channel (pedestal subtracted)");
   spe->Draw();
   spe->Draw("text same");
  spe2->Draw("same");

  allC->Modified(); allC->Update(); gSystem->ProcessEvents();
  //allC2->Modified(); allC2->Update(); gSystem->ProcessEvents();
  allC3->Modified(); allC3->Update(); gSystem->ProcessEvents();
  //allCm->Modified(); allC->Update(); gSystem->ProcessEvents();//maria



///// OUTPUT HISTOGRAMS .png //////////////////////////////////////////////////////
//TCanvas *c1 = new TCanvas("c1", "c1",10,10,1000,750);
//c1->SetLogz();
//gStyle->SetOptStat(11);
//gStyle->SetStatX(0.60);

ofstream fw;
//pedestals to text file
fw.open("../grinch_adc_calib/ped.txt", std::ofstream::trunc);
if (fw.is_open())
{
   for (int i = 0; i < nADC; i++) {
   fw << ped_chan[i] << "\n";
}
   fw.close();
}
   else cout << "Problem with opening 'pedestals.txt'"<<endl;

//spe channel to text file
fw.open("../grinch_adc_calib/spe.txt", std::ofstream::trunc);
if (fw.is_open())
{
   for (int i = 0; i < nADC; i++) {
   fw << spe_chan[i] << "\n";
}
   fw.close();
}
   else cout << "Problem with opening 'spe.txt'"<<endl;

fw.open("../grinch_adc_calib/spe_sub.txt", std::ofstream::trunc);
//ofstream fw("../adc_calib/pedestals.txt", std::ofstream::out);
if (fw.is_open())
{
   //subtracted spe channel to text file
   for (int i = 0; i < nADC; i++) {
   fw << spe_chan_sub[i] << "\n";
}
   fw.close();
}
   else cout << "Problem with opening 'spe_sub.txt'"<<endl;


//char outputname[200];

//GRINCH_tdc_LE_vs_pmt->Draw("COLZ");
//strcpy(outputname,"../monitoringPlots/tdc_LE_vs_pmt");
//c1->SaveAs(strcat(outputname,".png"));

//GRINCH_tdc_TE_vs_pmt->Draw("COLZ");
//strcpy(outputname,"../monitoringPlots/tdc_TE_vs_pmt");
//c1->SaveAs(strcat(outputname,".png"));

//GRINCH_mult_vs_pmt->Draw("COLZ");
//strcpy(outputname,"../monitoringPlots/tdc_mult_vs_pmt");
//c1->SaveAs(strcat(outputname,".png"));

//GRINCH_tdc_ToT_vs_pmt->Draw("COLZ");
//strcpy(outputname,"../monitoringPlots/tdc_ToT_vs_pmt");
//c1->SaveAs(strcat(outputname,".png"));

//GRINCH_adc_vs_chan->Draw("COLZ");
//strcpy(outputname,"../adcCalibPlots/adc_vs_chan");
//c1->SaveAs(strcat(outputname,".png"));


}

