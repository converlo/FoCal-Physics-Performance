#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>

void pTdistri(const char *filename="MergedAnalysisJets.root")
{
    TFile *file = TFile::Open(filename);

    gSystem->mkdir("JetDistributions", kTRUE);

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
        int R = radii[i];

        // pT distribution for reconstruct and truth jets
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

        c->SaveAs(Form("JetDistributions/RecoVsTruthJetPt_R%d.pdf",R));

        // eta distribution for reconstruct and truth jets
        TCanvas *cEta = new TCanvas(Form("cEta_R%d",R),
                            Form("Eta R=%d",R),
                            800,600);

        TH1F *hRecoEta = new TH1F(Form("hRecoEta_R%d",R),
                          Form("R=%d;Jet #eta;Entries",R),
                          100,-5.5,5.5);

        TH1F *hTruthEta = new TH1F(Form("hTruthEta_R%d",R),
                           Form("R=%d;Jet #eta;Entries",R),
                           100,-5.5,5.5);

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

        cEta->SaveAs(Form("JetDistributions/RecoVsTruthJetEta_R%d.pdf",R));

        // energy distribution for reconstruct and truth jets

        TCanvas *cEnergy = new TCanvas(Form("cEnergy_R%d",R),
                         Form("R=%d",R),
                         800,600);

        TH1F *hRecoEnergy = new TH1F(Form("hRecoEnergy_R%d",R),
                       Form("R=%d;Jet Energy (GeV);Entries",R),
                       220,0,220);

        TH1F *hTruthEnergy = new TH1F(Form("hTruthEnergy_R%d",R),
                        Form("R=%d;Jet Energy (GeV);Entries",R),
                        220,0,220);

        jetTree->Draw(Form("jetE>>hRecoEnergy_R%d",R),
              Form("jetR==%d",R),
              "goff");

        truthjetTree->Draw(Form("truthjetE>>hTruthEnergy_R%d",R),
                   Form("truthjetR==%d",R),
                   "goff");

        hRecoEnergy->SetLineColor(kBlue);
        hRecoEnergy->SetLineWidth(2);

        hTruthEnergy->SetLineColor(kRed);
        hTruthEnergy->SetLineWidth(2);

        hRecoEnergy->Draw("hist");
        hTruthEnergy->Draw("hist same");

        gPad->SetLogy();

        TLegend *legEnergy = new TLegend(0.60,0.75,0.88,0.88);
        legEnergy->AddEntry(hRecoEnergy,"Reco jets","l");
        legEnergy->AddEntry(hTruthEnergy,"Truth jets","l");
        legEnergy->Draw();

        cEnergy->SaveAs(Form("JetDistributions/RecoVsTruthJetEnergy_R%d.pdf",R));

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

        cPhi->SaveAs(Form("JetDistributions/RecoVsTruthJetPhi_R%d.pdf",R));
    }

    file->Close();
}
