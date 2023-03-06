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
#include<math.h>
#include <stack>



using namespace std;


//global variables
stack <Int_t> short_row_on_left_search;
stack <Int_t> short_row_on_right_search;
stack <Int_t> long_row_on_left_search;
stack <Int_t> long_row_on_right_search;
stack <Int_t> middle_search;
stack <Int_t> top_row_middle_search;
stack <Int_t> top_row_left_search;
stack <Int_t> top_row_right_search;
stack <Int_t> bottom_row_middle_search;
stack <Int_t> bottom_row_left_search;
stack <Int_t> bottom_row_right_search;

const Int_t N_ROW=60;

Int_t row_array[510];
Int_t col_array[510];

//function declarations
void make_row_col_array();
void Fill_Search_Stacks();
Bool_t Is_NOT_In_Stack(stack <Int_t> checkstack, Int_t PMT);
void Print_Stack(stack <Int_t> inputstack);
stack <Int_t> Clear_Stack();



void make_hist_grinch(TString basename="",Int_t nrun=2043,TString configfilename="run_list_grinch.txt", Bool_t show_track = kFALSE){
  if (basename=="") {
    cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
    cin >> basename;
  }
  TString fullname = Form("%s_%d",basename.Data(),nrun);
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
  TString inputroot;
  inputroot="rootfiles/"+fullname+".root";
  TString outputhist;
  outputhist= "hist/"+fullname+"_grinch_hist.root";
  TObjArray HList(0);
  TString Chainroot;
  configfilename="RunLists/"+configfilename;
  ifstream configfile(configfilename.Data());
  TString currentline;
  TChain *fchain = new TChain("T");
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ) {
     if( !currentline.BeginsWith("#") ){
      cout << " add file : " << currentline << endl;
      fchain->Add(currentline);
    }   
  } 
  //  Int_t hcalNum;
  //fchain->SetBranchAddress("Ndata.sbs.hcal.tdcelemID",&hcalNum) ;
  
  Double_t ps_e; // raw adc
  fchain->SetBranchAddress("bb.ps.e",&ps_e) ;
  Double_t sh_e; // raw adc
  fchain->SetBranchAddress("bb.sh.e",&sh_e) ;
  Double_t pmom; // raw adc
  fchain->SetBranchAddress("BB.gold.p",&pmom) ;
   Double_t xptar; // raw adc
  fchain->SetBranchAddress("BB.gold.th",&xptar) ;
   Double_t yptar; // raw adc
  fchain->SetBranchAddress("BB.gold.ph",&yptar) ;
   Double_t ytar; // raw adc
  fchain->SetBranchAddress("BB.gold.y",&ytar) ;
  Int_t GrinchNum;
  fchain->SetBranchAddress("Ndata.bb.grinch_tdc.tdcelemID",&GrinchNum) ;
   Double_t tdcGID[1000]; // raw adc
  Double_t tdcGLe[1000]; // raw adc
  Double_t tdcGTot[1000]; // raw adc
  Double_t tdcGTe[1000]; // trailing edge 
  
  Double_t tdcGMult[1000]; // raw adc
  fchain->SetBranchAddress("bb.grinch_tdc.tdcelemID",&tdcGID) ;
  fchain->SetBranchAddress("bb.grinch_tdc.tdc",&tdcGLe) ;
  fchain -> SetBranchAddress("bb.grinch_tdc.tdc_te",&tdcGTe) ;
  fchain->SetBranchAddress("bb.grinch_tdc.tdc_mult",&tdcGMult) ;
 fchain->SetBranchAddress("bb.grinch_tdc.tdc_tot",&tdcGTot) ;
  Int_t GrinchADCNum;
   Double_t adcGID[1000]; // raw adc
   Double_t adcGAtime[1000]; // raw adc
   Double_t adcGAmp[1000]; // raw adc
   Double_t grinch_adc[1000];//maria
   Double_t adcGMult[1000]; // raw adc //maria

   UInt_t fTrigBits;
   UInt_t fEvtNum;
   
 fchain->SetBranchAddress("Ndata.bb.grinch_adc.adcelemID",&GrinchADCNum) ;
 fchain->SetBranchAddress("bb.grinch_adc.adcelemID",&adcGID) ;
 fchain->SetBranchAddress("bb.grinch_adc.a_time",&adcGAtime) ;
 fchain->SetBranchAddress("bb.grinch_adc.a_amp_p",&adcGAmp) ; //amp_p
 fchain->SetBranchAddress("bb.grinch_adc.a_mult",&adcGMult) ;//maria
 fchain->SetBranchAddress("bb.grinch_adc.a_p",&grinch_adc) ; // ADC counts //maria a_p

 fchain->SetBranchAddress("fEvtHdr.fTrigBits",&fTrigBits);
 fchain->SetBranchAddress("fEvtHdr.fEvtNum", &fEvtNum);
 
  TH1F* h_ps_e = new TH1F("h_ps_e"," ; Total PreShower E ",300,0.0,2.0);
      HList.Add(h_ps_e);
  TH1F* h_TBB_e = new TH1F("h_TBB_e"," ; Total PRe+Shower E ",300,0.0,5.0);
      HList.Add(h_TBB_e);
  TH1F* h_ratio_e = new TH1F("h_ratio_e"," ; Total PRe+Shower E/ track norm ",300,0.0,2.5);
      HList.Add(h_ratio_e);
 

  TH2F* h_grinch_atime_elem = new TH2F("h_grinch_atime_elem"," ; GRINCH ADC elemID ; GRINCH ADC time (ns) ",63,0,63,500,0,2000);
      HList.Add(h_grinch_atime_elem);
 TH2F* h_grinch_atime_elem_multcut = new TH2F("h_grinch_atime_elem_multcut"," ; GRINCH ADC elemID ; GRINCH ADC time (ns) multcut ",63,0,63,500,0,2000);
 HList.Add(h_grinch_atime_elem_multcut);
 TH2F* h_grinch_atime_elem_SHcut  = new TH2F("h_grinch_atime_elem_SHcut"," ; GRINCH ADC elemID ; GRINCH ADC time (ns) Shower Cut ",63,0,63,500,0,2000);
 HList.Add(h_grinch_atime_elem_SHcut);

  TH2F* h_grinch_amp_elem = new TH2F("h_grinch_amp_elem"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) ",63,0,63,500,0,500);
      HList.Add(h_grinch_amp_elem);     
  TH2F* h_grinch_amp_elem_multcut = new TH2F("h_grinch_amp_elem_multcut"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) multcut ",63,0,63,600,0,600);
  HList.Add(h_grinch_amp_elem_multcut);//maria
 TH2F* h_grinch_amp_elem_SHcut = new TH2F("h_grinch_amp_elem_SHcut"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) Shower Cut ",63,0,63,600,0,600);
  HList.Add(h_grinch_amp_elem_SHcut);//maria
  TH2F* h_grinch_amp_elem_multcut_timecut = new TH2F("h_grinch_amp_elem_multcut_timecut"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) multcut timecut ",63,0,63,600,0,600);
   HList.Add(h_grinch_amp_elem_multcut_timecut);//maria
  TH2F* h_grinch_adcMult_elem = new TH2F("h_grinch_adcMult_elem"," ; GRINCH ADC elem id ; GRINCH ADC Mult. ",63,0,63,10,0,10);
      HList.Add(h_grinch_adcMult_elem);//maria
  TH2F* h_grinch_adcMult_elem_SHcut = new TH2F("h_grinch_adcMult_elem_SHcut"," ; GRINCH ADC elem id ; GRINCH ADC Mult. Shower Cut ",63,0,63,10,0,10);
      HList.Add(h_grinch_adcMult_elem_SHcut);//maria
  TH2F* h_grinch_le_elem = new TH2F("h_grinch_le_elem"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,2700,-100,2600);
      HList.Add(h_grinch_le_elem);
  TH2F* h_grinch_le_elem_multcut = new TH2F("h_grinch_le_elem_multcut"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) multcut",510,0,510,2700,-100,2600);
      HList.Add(h_grinch_le_elem_multcut);//maria

      TH2F* h_grinch_le_elem_SHcut = new TH2F("h_grinch_le_elem_SHcut"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) Shower Cut",510,0,510,2700,-100,2600);
      HList.Add(h_grinch_le_elem_SHcut);//maria

  TH2F* h_grinch_te_elem_multcut = new TH2F("h_grinch_te_elem_multcut"," ; GRINCH TDC elemID ; GRINCH TDC TE (ns) multcut",510,0,510,2700,-100,2600);
      HList.Add(h_grinch_te_elem_multcut);//maria

      TH2F* h_grinch_te_elem = new TH2F("h_grinch_te_elem"," ; GRINCH TDC elemID ; GRINCH TDC TE (ns)",510,0,510,2700,-100,2600);
      HList.Add(h_grinch_te_elem);//maria

      TH2F* h_grinch_te_elem_SHcut = new TH2F("h_grinch_te_elem_SHcut"," ; GRINCH TDC elemID ; GRINCH TDC TE (ns) Shower Cut",510,0,510,2700,-100,2600);
      HList.Add(h_grinch_te_elem_SHcut);//maria






  TH1F* h_grinch_ghits = new TH1F("h_grinch_ghits"," ; GRINCH Good hits ; ",10,0,10);
      HList.Add(h_grinch_ghits);
  TH1F* h_grinch_elem = new TH1F("h_grinch_elem"," ; GRINCH TDC elemID ; ",510,0,510);
      HList.Add(h_grinch_elem);
  TH1F* h_grinch_tot_all = new TH1F("h_grinch_tot_all"," ; GRINCH TDC Tot (no cut) ; ",50,0,50);
      HList.Add(h_grinch_tot_all);
  TH1F* h_grinch_tot_all_lecut = new TH1F("h_grinch_tot_all_lecut"," ; GRINCH TDC Tot (TDC LE cut) ; ",50,0,50);
      HList.Add(h_grinch_tot_all_lecut);
  TH2F* h_grinch_tot_elem = new TH2F("h_grinch_tot_elem"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) ",510,0,510,110,-10,100);
      HList.Add(h_grinch_tot_elem);
  TH2F* h_grinch_tot_elem_multcut = new TH2F("h_grinch_tot_elem_multcut"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) multcut ",510,0,510,110,-10,100);
      HList.Add(h_grinch_tot_elem_multcut); //maria
  TH2F* h_grinch_tot_elem_SHcut = new TH2F("h_grinch_tot_elem_SHcut"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) Shower Cut ",510,0,510,110,-10,100);
      HList.Add(h_grinch_tot_elem_SHcut); //maria
  TH2F* h_grinch_mult_elem = new TH2F("h_grinch_mult_elem"," ; GRINCH TDC elemID ; GRINCH TDC Mult ",510,0,510,15,0,15);
      HList.Add(h_grinch_mult_elem);

  TH2F* h_grinch_mult_elem_SHcut = new TH2F("h_grinch_mult_elem_SHcut"," ; GRINCH TDC elemID ; GRINCH TDC Mult Shower Cut ",510,0,510,15,0,15);
      HList.Add(h_grinch_mult_elem_SHcut);


  TH2F* h_grinch_adc_chan_vs_num = new TH2F("h_grinch_adc_chan_vs_num","; GRINCH ADC elemID ; GRINCH ADC channel ",64,0,63,510,-10,500);//maria
      HList.Add(h_grinch_adc_chan_vs_num);
 TH2F* h_grinch_adc_chan_vs_num_SHcut = new TH2F("h_grinch_adc_chan_vs_num_SHcut","; GRINCH ADC elemID ; GRINCH ADC channel Shower Cut",64,0,63,510,-10,500);//maria
      HList.Add(h_grinch_adc_chan_vs_num_SHcut);
  TH2F* h_grinch_adc_chan_vs_num_multcut = new TH2F("h_grinch_adc_chan_vs_num_multcut","; GRINCH ADC elemID ; GRINCH ADC channel multcut ",64,0,64,510,-10,500);//ma1ria
      HList.Add(h_grinch_adc_chan_vs_num_multcut);

      TH1F* h_EvtHdr_TrigBits = new TH1F("h_EvtHdr_TrigBits", "; fTrigBits;",35,0,35);
      HList.Add(h_EvtHdr_TrigBits);


 //
 TH1F* h_grinch_le[511];
 TH1F* h_grinch_tot[511];
 TH1F* h_grinch_tot_lecut[511];
 for (Int_t ig=0;ig<511;ig++) {
   h_grinch_le[ig] = new TH1F(Form("h_grinch_le_%d",ig),Form(" ; GRINCH TDC LE (ns) PMT %d  ; ",ig),25,875,925);
      HList.Add(h_grinch_le[ig]);
   h_grinch_tot[ig] = new TH1F(Form("h_grinch_tot_%d",ig),Form(" ; GRINCH TDC Tot PMT %d (no cut) ; ",ig),50,0,50);
      HList.Add(h_grinch_tot[ig]);
      h_grinch_tot_lecut[ig] = new TH1F(Form("h_grinch_tot_lecut_%d",ig),Form(" ; GRINCH TDC Tot PMT %d (LE cut) ; ",ig),50,0,50);
      HList.Add(h_grinch_tot_lecut[ig]);   
 }

 
TH1F* h_grinch_adc[64];//maria 
 for (Int_t ig=0;ig<64;ig++){
  h_grinch_adc[ig] = new TH1F(Form("h_grinch_adc_%d ",ig),Form(" ; GRINCH ADC Channels ADC %d ; ",ig), 201, -0.5, 200.5); 
   HList.Add(h_grinch_adc[ig]);
 }

TH1F* h_grinch_amp_p[64];//maria 
 for (Int_t ig=0;ig<64;ig++){
  h_grinch_amp_p[ig] = new TH1F(Form("h_grinch_amp_p_%d ",ig),Form(" ; GRINCH ADC Amplitude (pedsub) %d ; ",ig), 501, -0.5, 500.5); 
   HList.Add(h_grinch_amp_p[ig]);
 }

 Int_t zero_count = 0; 

      //
  Long64_t nentries = fchain->GetEntries();
  for (int i = 0; i < nentries; i++) {//nentries
    fchain->GetEntry(i);
    if (i%1000==0) cout << " Entry = " << i << endl;
    Double_t tot_e = ps_e+sh_e;
    Double_t rat = tot_e/pmom;
    h_ps_e->Fill(ps_e);
    h_TBB_e->Fill(tot_e);
    h_ratio_e->Fill(rat);
    h_EvtHdr_TrigBits ->Fill(fTrigBits);
  
    for (Int_t ig=0;ig<GrinchADCNum;ig++) {
      Int_t gindex=adcGID[ig];
      h_grinch_atime_elem->Fill(adcGID[ig],adcGAtime[ig]);
      h_grinch_amp_elem->Fill(adcGID[ig],adcGAmp[ig]);
      h_grinch_adc_chan_vs_num ->Fill(adcGID[ig],grinch_adc[ig]);//maria
      h_grinch_adcMult_elem ->Fill(adcGID[ig],adcGMult[ig]);//maria

      
      if (ps_e> 0.2 && rat > 0.5) {//cut that mark put in	
	h_grinch_amp_elem_SHcut ->Fill(adcGID[ig],adcGAmp[ig]);
	h_grinch_atime_elem_SHcut->Fill(adcGID[ig],adcGAtime[ig]);
	h_grinch_adc_chan_vs_num_SHcut ->Fill(adcGID[ig],grinch_adc[ig]);
	h_grinch_adcMult_elem_SHcut ->Fill(adcGID[ig],adcGMult[ig]);
      }
      

      if (adcGMult[ig]==1){
	h_grinch_amp_elem_multcut ->Fill(adcGID[ig],adcGAmp[ig]);
	h_grinch_atime_elem_multcut->Fill(adcGID[ig],adcGAtime[ig]);
	h_grinch_adc_chan_vs_num_multcut ->Fill(adcGID[ig],grinch_adc[ig]);

	if(adcGAtime[ig] > 1230 && adcGAtime[ig]<1300){
	  h_grinch_amp_elem_multcut_timecut->Fill(adcGID[ig],adcGAmp[ig]);
	}

	if (gindex <64 && gindex >-1 && adcGMult[ig]==1) {
	  h_grinch_adc[gindex]->Fill(grinch_adc[ig]);
	  h_grinch_amp_p[gindex] ->Fill(adcGAmp[ig]);
	}      
      }
    }


    Double_t goodhit=0;
    for (Int_t ig=0;ig<GrinchNum;ig++) {
      Int_t gindex=tdcGID[ig];
      if (gindex <511 && gindex >-1) {
      h_grinch_elem->Fill(tdcGID[ig]);
      h_grinch_le_elem->Fill(tdcGID[ig],tdcGLe[ig]);
      h_grinch_te_elem->Fill(tdcGID[ig],tdcGTe[ig]);
      h_grinch_tot_elem->Fill(tdcGID[ig],tdcGTot[ig]);
      h_grinch_mult_elem->Fill(tdcGID[ig],tdcGMult[ig]);
      h_grinch_tot_all->Fill(tdcGTot[ig]);
      h_grinch_le[gindex]->Fill(tdcGLe[ig]);
      
      if (ps_e> 0.2 && rat > 0.5) {//cut that mark put in
	h_grinch_le_elem_SHcut->Fill(tdcGID[ig],tdcGLe[ig]);
	h_grinch_te_elem_SHcut->Fill(tdcGID[ig],tdcGTe[ig]);
	h_grinch_tot_elem_SHcut->Fill(tdcGID[ig],tdcGTot[ig]);
	h_grinch_mult_elem_SHcut->Fill(tdcGID[ig],tdcGMult[ig]);
      }
      
      if (tdcGMult[ig]==1){//maria
	 h_grinch_le_elem_multcut->Fill(tdcGID[ig],tdcGLe[ig]);
	 h_grinch_te_elem_multcut->Fill(tdcGID[ig],tdcGTe[ig]);
	 h_grinch_tot_elem_multcut->Fill(tdcGID[ig],tdcGTot[ig]);
      }
     
	 if (abs(tdcGLe[ig]-900)<20) {
	   h_grinch_tot_lecut[gindex]->Fill(tdcGTot[ig]);
	 }
       
      h_grinch_tot[gindex]->Fill(tdcGTot[ig]);

      if (abs(tdcGLe[ig]-900)<20 && (ps_e+sh_e) >1.3){
	   goodhit++;
      }
      if (abs(tdcGLe[ig]-900)<20) {
	h_grinch_tot_all_lecut->Fill(tdcGTot[ig]);
      }
      
      }
      h_grinch_ghits->Fill(goodhit);
    }
  }
  TFile hsimc(outputhist,"recreate");
  HList.Write();

}// end main


void Fill_Search_Stacks()
{
  //// pre-filling stacks to make searching for neighbors and cluster finding cleaner in the other functions. 
  ////
  //// if [0] is the PMT we want to check the neighbors of, this is the PMT numbers of it's neighbors. 
  ////
  ////   [-9] [-8]
  //// [-1] [0] [+1]
  ////   [+8] [+9]
  ////
  ////
  //// So if we want to search around PMT 20, it would look like this:
  ////    [11] [12]
  ////  [19] [20] [21]
  ////    [28] [29]
  //// 
  //// If we are on the edge of the detector, we need to adjust slighty. 
  ////
  //// short_row_on_left : exclude [-1]
  //// short_row_on_right: exclude [+1]
  //// top_row_middle : exclue [-8], [-9]
  //// top_row_left : exclude [-1], [-8], [-9] this is only for PMT 0
  //// top_row_right: exclude [-9] ,[-8], [+1] this is only for PMT 7
  //// bottom_row_middle : exclude [+8], [+9] 
  //// bottom_row_left : exclude [-9], [-1], [+8], [+9] this is only for PMT 501
  //// bottom_row_right: exclude [-8], [+1], [+8], [+9] this is only for PMT 509
  ////

  Int_t search[6]= {-9,-8,9,8};
  for (Int_t counter = 0; counter < 4; counter ++)
    {
      short_row_on_left_search.push(search[counter]);
      short_row_on_right_search.push(search[counter]);
      middle_search.push(search[counter]);
    }
  short_row_on_left_search.push(1);
  middle_search.push(1);
  short_row_on_right_search.push(-1);
  middle_search.push(-1);
  
  long_row_on_left_search.push(-8);
  long_row_on_left_search.push(1);
  long_row_on_left_search.push(9);
  long_row_on_right_search.push(-9);
  long_row_on_right_search.push(-1);
  long_row_on_right_search.push(8);

  top_row_middle_search.push(-1);
  top_row_middle_search.push(1);
  top_row_middle_search.push(8);
  top_row_middle_search.push(9);

  top_row_left_search.push(8);
  top_row_left_search.push(1);
  top_row_left_search.push(9);

  top_row_right_search.push(-1);
  top_row_right_search.push(8);
  top_row_right_search.push(9);

  bottom_row_middle_search.push(-1);
  bottom_row_middle_search.push(1);
  bottom_row_middle_search.push(-8);
  bottom_row_middle_search.push(-9);  

  bottom_row_left_search.push(-8);
  bottom_row_left_search.push(1);

  bottom_row_right_search.push(-1);
  bottom_row_right_search.push(-9); 
}// end Fill_Search_Stacks

void make_row_col_array()// making it easy to convert from PMT num to row, col
{
  Int_t pmt = 0;
  Int_t N_COL;
  for (Int_t i = 0;i<N_ROW;i++)
    {
      N_COL = 8 + i%2;
      for (Int_t j=0 ; j< N_COL ; j++)
	{
	  row_array[pmt] = i;
	  col_array[pmt] = j;
	  pmt ++;
	}
    }
  /*
    for (Int_t r = 0; r< 510; r++)
    {
    cout<< "PMT "<<r<<": ROW "<< row_array[r] <<", COL"<< col_array[r]<<endl;
    }
  */
}// end make_row_col_array

Bool_t Is_NOT_In_Stack(stack <Int_t> checkstack, Int_t PMT)
{
  while (!checkstack.empty())
    {
      if (PMT == checkstack.top())
	{
	  return kFALSE; // the PMT is already in the stack
	}
      checkstack.pop();
    }
  return kTRUE; // the PMT is not already in the stack
}// end Is_NOT_In_Stack


void Print_Stack(stack <Int_t> inputstack)
{
  // cout<< "~~~~~~~"<<endl;
  cout<<endl;
  while (!inputstack.empty())
    {
      cout<<"| " <<inputstack.top()<< " |" <<endl;
      inputstack.pop();
    }
  
  cout<< "|_____|"<<endl;
}// end Print_Stack


stack <Int_t> Clear_Stack()
{
  stack <Int_t> clearstack; //default constructor is an empty stack
  return clearstack;
}// end Clear_Stack
