#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TROOT.h>

void pTdistri(const char *filename="MergedAnalysisJets.root")
{
    TFile *file = TFile::Open(filename);

    if(!file || file->IsZombie())
    {
        printf("Cannot open %s\n", filename);
        return;
    }

    TTree *jetTree = (TTree*)file->Get("jetTree");
    TTree *truthjetTree = (TTree*)file->Get("truthjetTree");

    if(!jetTree || !truthjetTree)
    {
        printf("Cannot find jetTree or truthjetTree\n");
        return;
    }

    int radii[] = {2,4,6};

    for(int i=0;i<3;i++)
    {
        // pT distribution for reconstruct and truth jets
        int R = radii[i];

        TCanvas *c = new TCanvas(Form("c_R%d",R),
                                 Form("R=%d",R),
                                 800,600);

        TH1F *hReco = new TH1F(Form("hReco_R%d",R),
                               Form("R=%d;Jet p_{T} (GeV/c);Entries",R),
                               220,0,220);

        TH1F *hTruth = new TH1F(Form("hTruth_R%d",R),
                                Form("R=%d;Jet p_{T} (GeV/c);Entries",R),
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

        TLegend *leg = new TLegend(0.60,0.75,0.88,0.88);
        leg->AddEntry(hReco,"Reco jets","l");
        leg->AddEntry(hTruth,"Truth jets","l");
        leg->Draw();

        c->SaveAs(Form("RecoVsTruthJetPt_R%d.pdf",R));

        // eta distribution for reconstruct and truth jets
        TCanvas *cEta = new TCanvas(Form("cEta_R%d",R),
                            Form("Eta R=%d",R),
                            800,600);

        TH1F *hRecoEta = new TH1F(Form("hRecoEta_R%d",R),
                          Form("R=%d;Jet #eta;Entries",R),
                          100,-1.5,1.5);

        TH1F *hTruthEta = new TH1F(Form("hTruthEta_R%d",R),
                           Form("R=%d;Jet #eta;Entries",R),
                           100,-1.5,1.5);

        jetTree->Draw(Form("jetEta>>hRecoEta_R%d",R),
              Form("jetR==%d",R),
              "goff");

        truthjetTree->Draw(Form("truthjetEta>>hTruthEta_R%d",R),
                   Form("truthjetR==%d",R),
                   "goff");

        hRecoEta->SetLineColor(kBlue);
        hRecoEta->SetLineWidth(2);

        hTruthEta->SetLineColor(kRed);
        hTruthEta->SetLineWidth(2);

        hRecoEta->Draw("hist");
        hTruthEta->Draw("hist same");

        TLegend *legEta = new TLegend(0.60,0.75,0.88,0.88);
        legEta->AddEntry(hRecoEta,"Reco jets","l");
        legEta->AddEntry(hTruthEta,"Truth jets","l");
        legEta->Draw();

        cEta->SaveAs(Form("RecoVsTruthJetEta_R%d.pdf",R));

        // phi distribution for reconstruct and truth jets
        TCanvas *cPhi = new TCanvas(Form("cPhi_R%d",R),
                            Form("Phi R=%d",R),
                            800,600);

        TH1F *hRecoPhi = new TH1F(Form("hRecoPhi_R%d",R),
                          Form("R=%d;Jet #phi (rad);Entries",R),
                          64,-TMath::Pi(),TMath::Pi());

        TH1F *hTruthPhi = new TH1F(Form("hTruthPhi_R%d",R),
                           Form("R=%d;Jet #phi (rad);Entries",R),
                           64,-TMath::Pi(),TMath::Pi());

        jetTree->Draw(Form("jetPhi>>hRecoPhi_R%d",R),
              Form("jetR==%d",R),
              "goff");

        truthjetTree->Draw(Form("truthjetPhi>>hTruthPhi_R%d",R),
                   Form("truthjetR==%d",R),
                   "goff");

        hRecoPhi->SetLineColor(kBlue);
        hRecoPhi->SetLineWidth(2);

        hTruthPhi->SetLineColor(kRed);
        hTruthPhi->SetLineWidth(2);

        hRecoPhi->Draw("hist");
        hTruthPhi->Draw("hist same");

        TLegend *legPhi = new TLegend(0.60,0.75,0.88,0.88);
        legPhi->AddEntry(hRecoPhi,"Reco jets","l");
        legPhi->AddEntry(hTruthPhi,"Truth jets","l");
        legPhi->Draw();

        cPhi->SaveAs(Form("RecoVsTruthJetPhi_R%d.pdf",R));
    }

    file->Close();
}
