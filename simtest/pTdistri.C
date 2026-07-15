for (int R : {2,4,6})
{
    TCanvas *c = new TCanvas(Form("c_R%d",R),
                             Form("R = %d",R),
                             800,600);

    TH1F *hReco = new TH1F(Form("hReco_R%d",R),
                           Form("R = %d;Jet p_{T} (GeV/c);Entries",R),
                           220,0,220);

    TH1F *hTruth = new TH1F(Form("hTruth_R%d",R),
                            Form("R = %d;Jet p_{T} (GeV/c);Entries",R),
                            220,0,220);

    jetTree->Draw(Form("jetpT>>hReco_R%d",R),
                  Form("jetR==%d",R),
                  "goff");

    truthjetTree->Draw(Form("truthjetpT>>hTruth_R%d",R),
                       Form("truthjetR==%d",R),
                       "goff");

    hReco->SetLineColor(kBlue);
    hReco->SetLineWidth(2);

    hTruth->SetLineColor(kRed);
    hTruth->SetLineWidth(2);

    hReco->Draw("hist");
    hTruth->Draw("hist same");

    gPad->SetLogy();

    auto leg = new TLegend(0.60,0.75,0.88,0.88);
    leg->AddEntry(hReco,"Reco jets","l");
    leg->AddEntry(hTruth,"Truth jets","l");
    leg->Draw();

    c->SaveAs(Form("RecoVsTruthJetPt_R%d.pdf",R));
}
