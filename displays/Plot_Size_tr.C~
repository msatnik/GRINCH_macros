#include <TH1D.h>
#include <TAttFill.h>

void Plot_Size(){
  
  TCanvas *c1 = new TCanvas("c1","GRINCH cluster size",1000,500);
  c1->Divide(2,1);
  c1->SetLogy();

  //for( Int_t i=0; i<nkine; i++ ){

    c1->cd(1);

    gStyle->SetOptStat(0);
    //Get pre-ecal energy dist
    TFile *f1 = TFile::Open("../output/sbs8_13487_ele.root");
    TH1D *h1= (TH1D*)f1->Get("h_BBgr_clus_size");
    TH1D *h1_norm =  (TH1D*)(h1 ->Clone("h1_norm"));
    h1_norm ->Scale(1./h1_norm->Integral(),"width");
    //h1_norm->GetXaxis()->SetRange(0,20);
    h1_norm->SetLineWidth(3);
    h1_norm ->SetFillColor(30);
    h1_norm ->SetFillStyle(3003);
    h1_norm->SetLineColor(30);
    //h1->GetXaxis()->SetRange(0,20);
    h1 ->SetLineWidth(3);
    h1->SetFillColor(30);
    h1 ->SetFillStyle(3003);
    h1 ->SetLineColor(30);
    
    
   
    TFile *f2 = TFile::Open("../output/sbs8_13487_pions.root");
    TH1D *h2 = (TH1D*)f2->Get("h_BBgr_clus_size");
    TH1D *h2_norm =  (TH1D*)(h2 ->Clone("h2_norm"));
    h2_norm ->Scale(1./h2_norm->Integral(),"width");
    // h2_norm->GetXaxis()->SetRange(0,20);
    h2_norm->SetLineWidth(2);
    h2_norm->SetLineColor(kMagenta - 6);
    h2_norm ->SetFillColor(kMagenta - 6);
    h2_norm ->SetFillStyle(3003);
    //h2 ->GetXaxis()->SetRange(0,20);
    h2 ->SetLineWidth(2);
    h2 ->SetLineColor(kMagenta - 6);
    h2 ->SetFillColor(kMagenta - 6);
    h2  ->SetFillStyle(3003);


			     
    //h2->Draw("hist same");

    // TF1 *fit1;
    // h1->Fit("gaus","","",fit1l[i],fit1u[i]);
    // fit1 = h1->GetFunction("gaus");
    // Double_t fit1m = fit1->GetParameter(1);

    // TF1 *fit2;
    // h2->Fit("gaus","","",fit2l[i],fit2u[i]);
    // fit2 = h2->GetFunction("gaus");
    // Double_t fit2m = fit2->GetParameter(1);
    // Double_t fit2s = fit2->GetParameter(2);
    
    h1->Draw("hist");
    h2->Draw("hist same");
   

    //Add a legend
    auto legend = new TLegend(0.43,0.7,0.89,0.89);
    legend->SetTextSize(0.035);
    legend->SetHeader("GRINCH cluster size");
    legend->AddEntry(h1,"ps.e > 0.2","l");
    legend->AddEntry(h2,"ps.e < 0.2","l");
    // legend->AddEntry(h3_norm,"No Offsets","l");
    legend->Draw();

    c1->cd(2);
     
    h1_norm->Draw("hist");
     h2_norm->Draw("hist same");
     // h3_2_norm ->Draw("hist same");
     auto legend2 = new TLegend(0.43,0.7,0.89,0.89);
     legend2->SetTextSize(0.035);
     legend2->SetHeader("NORMALIZED");
     legend2->AddEntry(h1_norm,"ps.e > 0.2","l");
     legend2->AddEntry(h2_norm,"ps.e < 0.2","l");
     legend2->Draw();


    
    //}

  // c1->SaveAs("/u/home/sseeds/Plots/ecalE_allkine.png");

}
