void Plot_Variance(){
  
  TCanvas *c1 = new TCanvas("c1","GRINCH cluster variance",1600,1200);
  c1->Divide(2,1);

  //for( Int_t i=0; i<nkine; i++ ){

    c1->cd(2);

    gStyle->SetOptStat(0);
    //Get pre-ecal energy dist
    TFile *f1 = TFile::Open("grinch_3038_corrected.root");
    TH1D *h1= (TH1D*)f1->Get("h_grinch_cluster_variance");
    TH1D *h1_2 = (TH1D*)f1->Get("h_grinch_cluster_tmean");
    TH1D *h1_norm =  (TH1D*)(h1 ->Clone("h1_norm"));
    h1_norm ->Scale(1./h1_norm->Integral(),"width");
    TH1D *h1_2_norm =  (TH1D*)(h1_2 ->Clone("h1_2_norm"));
    h1_2_norm ->Scale(1./h1_2_norm->Integral(),"width");
    h1_norm->GetXaxis()->SetRange(0,24);
    h1_norm->SetLineWidth(3);
    h1_norm->SetLineColor(kBlack);
    h1_2_norm->SetLineWidth(3);
    h1_2_norm->SetLineColor(kBlack);
    
    
    
    //h1->Draw("hist");
    //Get post-ecal energy dist
    TFile *f2 = TFile::Open("grinch_3038_uncorrected.root");
    TH1D *h2 = (TH1D*)f2->Get("h_grinch_cluster_variance");
    TH1D *h2_2 = (TH1D*)f2->Get("h_grinch_cluster_tmean");
    TH1D *h2_norm =  (TH1D*)(h2 ->Clone("h2_norm"));
    h2_norm ->Scale(1./h2_norm->Integral(),"width");
    TH1D *h2_2_norm =  (TH1D*)(h2_2 ->Clone("h2_2_norm"));
    h2_2_norm ->Scale(1./h2_2_norm->Integral(),"width");
    h2_norm->GetXaxis()->SetRange(0,24);
    h2_norm->SetLineWidth(2);
    h2_norm->SetLineColor(kRed);
    h2_2_norm->SetLineWidth(2);
    h2_2_norm->SetLineColor(kRed);


    TFile *f3 = TFile::Open("grinch_3036.root");
    TH1D *h3 = (TH1D*)f3->Get("h_grinch_cluster_variance");
    TH1D *h3_2 = (TH1D*)f3->Get("h_grinch_cluster_tmean");
    TH1D *h3_norm =  (TH1D*)(h3 ->Clone("h3_norm"));
    h3_norm ->Scale(1./h3_norm->Integral(),"width");
    TH1D *h3_2_norm =  (TH1D*)(h3_2 ->Clone("h3_2_norm"));
    h3_2_norm ->Scale(1./h3_2_norm->Integral(),"width");
    h3_norm->GetXaxis()->SetRange(0,24);
    h3_norm->SetLineWidth(2);
    h3_norm->SetLineColor(kBlue);
    h3_2_norm->SetLineWidth(2);
    h3_2_norm->SetLineColor(kBlue);
 
 

			     
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
    
    h1_norm->Draw("hist");
    h2_norm->Draw("hist same");
    //h3_norm->Draw("hist same");

    // Double_t res = fit2s/fit2m*100;

    //Add a legend
    auto legend = new TLegend(0.43,0.7,0.89,0.89);
    legend->SetTextSize(0.02);
    legend->SetHeader("GRINCH cluster variance");
    legend->AddEntry(h1_norm,"Additional Offsets","l");
    legend->AddEntry(h2_norm,"LED Offsets","l");
    // legend->AddEntry(h3_norm,"No Offsets","l");
    legend->Draw();

    c1->cd(1);
     
     h1_2_norm->Draw("hist");
     h2_2_norm->Draw("hist same");
     // h3_2_norm ->Draw("hist same");
     auto legend2 = new TLegend(0.43,0.7,0.89,0.89);
     legend2->SetTextSize(0.02);
     legend2->SetHeader("GRINCH cluster mean time");
     legend2->AddEntry(h1_2_norm,"Additional Offsets","l");
     legend2->AddEntry(h2_2_norm,"LED Offsets","l");
     // legend2->AddEntry(h3_2_norm,"No Offsets","l");
     legend2->Draw();


    
    //}

  // c1->SaveAs("/u/home/sseeds/Plots/ecalE_allkine.png");

}
