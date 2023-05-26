void Plot_W2(){
  
  TCanvas *c1 = new TCanvas("c1","GRINCH",1000,500);
  c1->Divide(2,2);
  c1->SetLogy();

  //for( Int_t i=0; i<nkine; i++ ){

    c1->cd(1);

    gStyle->SetOptStat(0);
    //Get pre-ecal energy dist
    TFile *f1 = TFile::Open("../output/sbs9_big.root");
    TH1D *h1= (TH1D*)f1->Get("h_W2_gr_anticut");
    //TH1D *h1_norm =  (TH1D*)(h1 ->Clone("h1_norm"));
    //h1_norm ->Scale(1./h1_norm->Integral(),"width");
    //h1_norm->GetXaxis()->SetRange(1,200);
    //h1_norm->SetLineWidth(3);
    //h1_norm ->SetFillColor(30);
    //h1_norm ->SetFillStyle(3003);
    //h1_norm->SetLineColor(30);
    h1 ->SetLineWidth(3);
    h1 ->SetLineColor(30);
    // h1->SetFillColor(30);
    //h1 ->SetFillStyle(3003);
   
    
    
    TH1D *h2 = (TH1D*)f1->Get("h_W2_ps_anticut");
    //TH1D *h2_norm =  (TH1D*)(h2 ->Clone("h2_norm"));
    //h2_norm ->Scale(1./h2_norm->Integral(),"width");
    //h2_norm->GetXaxis()->SetRange(1,200);
    //h2_norm->SetLineWidth(2);
    //h2_norm->SetLineColor(kMagenta - 6);
    //h2_norm ->SetFillColor(kMagenta - 6);
    //h2_norm ->SetFillStyle(3003);

    //h2 ->GetXaxis()->SetRange(1,200);
    h2 ->SetLineWidth(2);
    h2 ->SetLineColor(kMagenta - 6);
    //h2 ->SetFillColor(kMagenta - 6);
    //h2 ->SetFillStyle(3003);
    


    TH1D *h3 = (TH1D*)f1->Get("h_W2_gr_ps_anticut");
    h3 ->SetLineWidth(2);
    h3 ->SetLineColor(kBlue - 6);
    //h3 ->SetFillColor(kBlue - 6);
    //h3 ->SetFillStyle(3003);
    TH1D *h3_scaled =  (TH1D*)(h3 ->Clone("h3_scaled"));
    h3_scaled ->SetLineWidth(2);
    h3_scaled->SetLineColor(kViolet - 6);
    //h3_scaled ->SetFillColor(kViolet - 6);
    //h3_scaled ->SetFillStyle(3003);
    h3_scaled ->Scale(3.4);



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
    h3->Draw("hist same");
    h3_scaled ->Draw("hist same");
   

    //Add a legend
    auto legend = new TLegend(0.43,0.7,0.89,0.89);
    legend->SetTextSize(0.035);
    legend->SetHeader("W2 with sideband cuts");
    legend->AddEntry(h1,"W2 w/ cut gr_ToT_sum < 30","l");
    legend->AddEntry(h2,"W2 w/ cut ps.e < 0.2","l");
    legend->AddEntry(h3,"W2 w/ both cuts","l");
    legend->Draw();

    c1->cd(2);

    TH1D *h4 =  (TH1D*)(h1 ->Clone("h4"));
    h4 ->Add(h3_scaled,-1);
    h4 ->SetLineWidth(3);
    h4 ->SetLineColor(30);
    h4 ->Draw("hist");

   
     auto legend2 = new TLegend(0.43,0.7,0.89,0.89);
     legend2->SetTextSize(0.035);
     legend2->SetHeader("subtracted");
     legend2->AddEntry(h4,"gr cut - background","l");
     // legend2->Draw();


     TH1D *h5 = (TH1D*)f1->Get("h_W2");
      h5 ->SetLineWidth(2);
      h5 ->SetLineColor(kBlue - 6);
      TH1D *h3_scaled_big =  (TH1D*)(h3 ->Clone("h3_scaled_big"));
      h3_scaled_big -> Scale(14);
      h3_scaled_big ->SetLineWidth(2);
      h3_scaled_big ->SetLineColor(kViolet - 6);

      c1 ->cd(3);
      h5->Draw("hist");
      h3_scaled_big ->Draw("hist same");

      auto legend3 = new TLegend(0.43,0.7,0.89,0.89);
     legend3->SetTextSize(0.035);
     legend3->SetHeader("W2");
     legend3->AddEntry(h5,"W2","l");
      legend3->AddEntry(h3_scaled_big,"scaled background","l");
     legend3->Draw();


     c1 ->cd(4);      
     TH1D *h6 =  (TH1D*)(h5 ->Clone("h6"));
     h6 ->Add(h3_scaled_big,-1);
     h6 ->SetLineWidth(3);
     h6 ->SetLineColor(30);
     h6 ->Draw("hist");
     
      auto legend4= new TLegend(0.43,0.7,0.89,0.89);
     legend4->SetTextSize(0.035);
     legend4->SetHeader("subtracted");
     legend4->AddEntry(h4,"W2 - background","l");
     // legend4->Draw();


    
    //}

  // c1->SaveAs("/u/home/sseeds/Plots/ecalE_allkine.png");

}
