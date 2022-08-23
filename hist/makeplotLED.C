void makeplotLED(int run)
{
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
  TFile *_file0 = TFile::Open(Form("grinchled_%d_grinch_hist.root",run));
  cout<<gFile->GetName()<<endl;

  TLine *line= new TLine(16,0,16,500);
  line -> SetLineColor(kBlack);
  line ->SetLineStyle(9);
  line ->SetLineWidth(2);
  TLine *line2= new TLine(32,0,32,500);
  line2 -> SetLineColor(kBlack);
  line2 -> SetLineStyle(9);
  line2 ->SetLineWidth(2);
  TLine *line3= new TLine(48,0,48,500);
  line3 -> SetLineColor(kBlack);
  line3 -> SetLineStyle(9);
  line3 ->SetLineWidth(2);
  TLine *linetdc =  new TLine(192,0,192,2500);
  linetdc -> SetLineColor(kBlack);
  linetdc ->SetLineStyle(9);
  linetdc ->SetLineWidth(2);
  TLine *linetdc2 =  new TLine(255,0,255,2500);
  linetdc2 -> SetLineColor(kBlack);
  linetdc2 ->SetLineStyle(9);
  linetdc2 ->SetLineWidth(2);

  TLatex * title = new TLatex(0.15,0.94,Form("Run %d",run));
  title->SetNDC(kTRUE);


  /*
  TCanvas * c5 = new TCanvas;
  c5 -> SetCanvasSize(1800,2000);
  c5 -> SetWindowSize(1800,800);
  c5 -> Divide(8,8);
  //TH1F* temphisto = (TH1F*) gFile->Get("h_grinch_adc_0");
  for (int i=0;i<64;i++){
    TH1F* temphisto = (TH1F*) gFile->Get(Form("h_grinch_adc_%d",i));
    c5 -> cd(i+1) ->SetLogy();  
    temphisto ->GetXaxis()->SetRange(0,300);
    temphisto -> Draw("same");
  }
  */

  

  TCanvas * cADC =  new TCanvas;
  cADC -> SetCanvasSize(1300,1050);
  cADC -> SetWindowSize(1350,1100);
  cADC -> Divide(2,2);
  TH2F * adcMult = (TH2F*) gFile->Get("h_grinch_adcMult_elem");
  TH2F* amp = (TH2F*) gFile->Get("h_grinch_amp_elem");
  TH2F* adcInt = (TH2F*) gFile->Get("h_grinch_adc_chan_vs_num");
  TH2F* atime =  (TH2F*) gFile->Get("h_grinch_atime_elem");
  cADC -> cd(1) ->SetLogz(1);
  amp -> Draw("colz");
  title->Draw();  
  line -> Draw();
  line2 -> Draw();
  line3 -> Draw();
  cADC -> cd(2) ->SetLogz(1);
  adcInt -> Draw("colz"); 
  title->Draw();  
  line -> Draw();
  line2 -> Draw();
  line3 -> Draw();
  cADC -> cd(3) ->SetLogz(1);
  atime -> Draw("colz");
  //atime ->GetXaxis() ->SetRange(0,63);
  title->Draw();  
  line -> Draw();
  line2 -> Draw();
  line3 -> Draw();
  cADC -> cd(4)-> SetLogz(1);
  //adcMult ->GetXaxis() ->SetRange(0,63);
  adcMult -> Draw("colz");
  title->Draw();  
  line -> Draw();
  line2 -> Draw();
  line3 -> Draw();


  TCanvas * cADC_multcut =  new TCanvas;
  cADC_multcut -> SetCanvasSize(1300,1050);
  cADC_multcut -> SetWindowSize(1350,1100);
  cADC_multcut -> Divide(2,2);
  TH2F* amp_multcut = (TH2F*) gFile->Get("h_grinch_amp_elem_multcut");
  TH2F* adcInt_multcut = (TH2F*) gFile->Get("h_grinch_adc_chan_vs_num_multcut");
  TH2F* atime_multcut =  (TH2F*) gFile->Get("h_grinch_atime_elem_multcut");
  cADC_multcut -> cd(1) ->SetLogz(1);
  amp_multcut -> Draw("colz");
  title->Draw();  
  line -> Draw();
  line2 -> Draw();
  line3 -> Draw();
  cADC_multcut -> cd(2) ->SetLogz(1);
  adcInt_multcut -> Draw("colz"); 
  title->Draw();  
  line -> Draw();
  line2 -> Draw();
  line3 -> Draw();
  cADC_multcut -> cd(3) ->SetLogz(1);
  atime_multcut -> Draw("colz");
  title->Draw();  
  line -> Draw();
  line2 -> Draw();
  line3 -> Draw();
  cADC_multcut -> cd(4)-> SetLogz(1);
  adcMult -> Draw("colz");
  title->Draw();  
  line -> Draw();
  line2 -> Draw();
  line3 -> Draw();

 



  
  TCanvas *cTDC =  new TCanvas;
  cTDC -> SetCanvasSize(1300,1050);
  cTDC -> SetWindowSize(1350,1100);
  cTDC -> Divide(2,2);
  TH2F* tot = (TH2F*) gFile->Get("h_grinch_tot_elem");
  TH2F* histoLE = (TH2F*) gFile->Get("h_grinch_le_elem");
  TH2F* histoTE = (TH2F*) gFile->Get("h_grinch_te_elem");
  TH2F* tdcMult = (TH2F*) gFile->Get("h_grinch_mult_elem");
  cTDC -> cd(1)->SetLogz(1);
  histoLE -> Draw("colz");
  title->Draw(); 
  linetdc -> Draw();
  linetdc2 -> Draw();
  cTDC -> cd(2)->SetLogz(1);
  histoTE -> Draw("colz");
  title->Draw(); 
  linetdc -> Draw();
  linetdc2 -> Draw();
  cTDC -> cd(3)-> SetLogz(1);
  tot -> Draw("colz");
  title->Draw(); 
  linetdc -> Draw();
  linetdc2 -> Draw();
  cTDC -> cd(4) ->SetLogz(1);
  tdcMult -> Draw("colz");
  title->Draw();
  linetdc -> Draw();
  linetdc2 -> Draw();
  
  
  
  TCanvas *cTDC_mult =  new TCanvas;
  cTDC_mult -> SetCanvasSize(1300,1050);
  cTDC_mult -> SetWindowSize(1350,1100);
  cTDC_mult -> Divide(2,2);
  TH2F* tot_mult = (TH2F*) gFile->Get("h_grinch_tot_elem_multcut");
  TH2F* histoLE_mult = (TH2F*) gFile->Get("h_grinch_le_elem_multcut");
  TH2F* histoTE_mult = (TH2F*) gFile->Get("h_grinch_te_elem_multcut");
  cTDC_mult -> cd(1)->SetLogz(1);
  histoLE_mult -> Draw("colz");
  title->Draw(); 
  linetdc -> Draw();
  linetdc2 -> Draw();
  cTDC_mult -> cd(2)->SetLogz(1);
  histoTE_mult -> Draw("colz");
  title->Draw(); 
  linetdc -> Draw();
  linetdc2 -> Draw();
  cTDC_mult -> Draw();
  linetdc2 -> Draw();
  cTDC_mult -> cd(3)-> SetLogz(1);
  tot_mult -> Draw("colz");
  title->Draw(); 
  linetdc -> Draw();
  linetdc2 -> Draw();
  cTDC_mult -> cd(4) ->SetLogz(1);
  tdcMult -> Draw("colz");
  title->Draw();
  linetdc -> Draw();
  linetdc2 -> Draw();

   
 
}
