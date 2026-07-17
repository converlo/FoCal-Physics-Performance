#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TSystem.h>
#include <vector>
#include <iostream>
#include <functional>

using namespace std;

// Delta-phi brut, sans repliement, entre 0 et 2pi (comme demande)
double DeltaPhiRaw(double phi1, double phi2)
{
    return fabs(phi1 - phi2);
}

// Repliement dans [0, pi], utilise uniquement en interne pour calculer une
// distance geometrique DeltaR correcte (independant du choix d'affichage de Delta-phi)
double DeltaPhiFold(double phi1, double phi2)
{
    double dphi = fabs(phi1 - phi2);
    if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
    return dphi;
}

double DeltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = DeltaPhiFold(phi1, phi2);
    return TMath::Sqrt(deta*deta + dphi*dphi);
}

void AnalyzeDijetDeltaPhi(const char *filename="MergedAnalysisJets.root",
                           const int R=2,
                           const int pTcut=10,
                           const int pT2cut=0,        // (2) coupure sur le subleading jet
                           const double dRcut=0.0,     // (3) separation geometrique min. leading-subleading (0 = desactive)
                           const bool require2jets=false) // (4) ne garder que les evenements a exactement 2 jets reconstruits
{
    //==========================================================
    // Open input file
    //==========================================================

    TFile *file = TFile::Open(filename);

    gSystem->mkdir("DijetAnalysis", kTRUE);

    if(!file || file->IsZombie())
    {
        cout << "Cannot open " << filename << endl;
        return;
    }

    TTree *jetTree = (TTree*)file->Get("jetTree");

    if(!jetTree)
    {
        cout << "jetTree not found." << endl;
        return;
    }

    //==========================================================
    // Branches
    //==========================================================

    Int_t ievt;
    Float_t jetpT;
    Float_t jetPhi;
    Float_t jetEta;
    Int_t jetR;

    jetTree->SetBranchAddress("ievt",&ievt);
    jetTree->SetBranchAddress("jetpT",&jetpT);
    jetTree->SetBranchAddress("jetPhi",&jetPhi);
    jetTree->SetBranchAddress("jetEta",&jetEta);
    jetTree->SetBranchAddress("jetR",&jetR);

    //==========================================================
    // Histograms
    //==========================================================

    TH1D *hDeltaPhi = new TH1D("hDeltaPhi", ";#Delta#phi (rad);Events", 64, 0, 2*TMath::Pi());

    TH1D *hpT1 = new TH1D("hpT1", ";Leading jet p_{T} (GeV/c);Events", 100,0,300);
    TH1D *hpT2 = new TH1D("hpT2", ";Subleading jet p_{T} (GeV/c);Events", 100,0,300);
    TH1D *hDeltapT = new TH1D("hDeltapT", ";p_{T1}-p_{T2} (GeV/c);Events", 100,0,200);
    TH1D *hDeltaEta = new TH1D("hDeltaEta", ";#Delta#eta;Events", 100,0,3);

    TH2D *hpT1pT2 = new TH2D("hpT1pT2",
        ";Leading jet p_{T} (GeV/c);Subleading jet p_{T} (GeV/c)",
        100, 0, 300, 100, 0, 300);

    TH2D *hDeltaEtaDeltaPhi = new TH2D("hDeltaEtaDeltaPhi",
        ";#Delta#eta;#Delta#phi (rad)",
        100, 0, 3, 64, 0, 2*TMath::Pi());

    // (1) plots de correlation phi1-phi2 et eta1-eta2
    TH2D *hPhiCorr = new TH2D("hPhiCorr",
        ";Leading jet #phi (rad);Subleading jet #phi (rad)",
        64, 0, 2*TMath::Pi(), 64, 0, 2*TMath::Pi());

    TH2D *hEtaCorr = new TH2D("hEtaCorr",
        ";Leading jet #eta;Subleading jet #eta",
        100, -7, 7, 100, -7, 7);

    // (3) ΔR entre leading et subleading, utile pour verifier l'effet de dRcut
    TH1D *hDeltaR = new TH1D("hDeltaR", ";#DeltaR(leading,subleading);Events", 100, 0, 8);

    // (4) nombre de jets reconstruits par evenement (avant toute selection dijet)
    TH1D *hNjets = new TH1D("hNjets", ";N_{jets}/event;Events", 21, -0.5, 20.5);

    //==========================================================
    // Fonction de traitement d'un evenement complet
    //==========================================================

    Long64_t nEventsSeen = 0;
    Long64_t nEventsKept  = 0;

    auto processEvent = [&](const vector<float>& pT,
                             const vector<float>& phi,
                             const vector<float>& eta,
                             int evtId)
    {
        if(pT.empty()) return;

        nEventsSeen++;
        hNjets->Fill((double)pT.size());

        // (4) selection stricte: on ne garde que les evenements a exactement 2 jets reconstruits
        if(require2jets && pT.size() != 2) return;

        if(pT.size() < 2) return;

        // Tri des indices par pT decroissant (permet de "sauter" un jet, cf point 3)
        vector<int> idx(pT.size());
        for(unsigned int j=0;j<pT.size();j++) idx[j]=j;
        sort(idx.begin(), idx.end(), [&](int a, int b){ return pT[a] > pT[b]; });

        int leadIdx = idx[0];
        if(pT[leadIdx] < pTcut) return;

        // (3) selection geometrique back-to-back: on prend le premier jet, dans l'ordre
        // decroissant de pT, qui est separe du leading par au moins dRcut. Si dRcut==0,
        // ce filtre est desactive et on prend simplement le 2e jet le plus energetique.
        int subIdx = -1;
        for(unsigned int k=1;k<idx.size();k++)
        {
            int cand = idx[k];
            double dR = DeltaR(eta[leadIdx], phi[leadIdx], eta[cand], phi[cand]);
            if(dR >= dRcut)
            {
                subIdx = cand;
                break;
            }
        }
        if(subIdx < 0) return; // aucun jet suffisamment separe du leading

        double pt1 = pT[leadIdx];
        double pt2 = pT[subIdx];

        // (2) coupure sur le subleading jet
        if(pt2 < pT2cut) return;

        double dphi = DeltaPhiRaw(phi[leadIdx], phi[subIdx]);
        double deta = fabs(eta[leadIdx] - eta[subIdx]);
        double dR   = DeltaR(eta[leadIdx], phi[leadIdx], eta[subIdx], phi[subIdx]);

        if(evtId < 10)
        {
            cout << "Event " << evtId
                 << "  Njets = " << pT.size()
                 << "  lead pT = " << pt1
                 << "  sublead pT = " << pt2
                 << "  phi1 = " << phi[leadIdx]
                 << "  phi2 = " << phi[subIdx]
                 << "  dphi = " << dphi
                 << "  dR = " << dR
                 << endl;
        }

        nEventsKept++;

        hDeltaPhi->Fill(dphi);
        hpT1->Fill(pt1);
        hpT2->Fill(pt2);
        hDeltapT->Fill(pt1 - pt2);
        hDeltaEta->Fill(deta);
        hDeltaR->Fill(dR);
        hpT1pT2->Fill(pt1, pt2);
        hDeltaEtaDeltaPhi->Fill(deta, dphi);
        hPhiCorr->Fill(phi[leadIdx], phi[subIdx]);
        hEtaCorr->Fill(eta[leadIdx], eta[subIdx]);
    };

    //==========================================================
    // Variables / boucle sur l'arbre
    //==========================================================

    vector<float> eventpT, eventPhi, eventEta;

    Long64_t nEntries = jetTree->GetEntries();
    cout << "Number of jets = " << nEntries << endl;

    Int_t currentEvent = -1;

    // Cherche le premier evenement ayant le rayon demande
    for(Long64_t i=0; i<nEntries; i++)
    {
        jetTree->GetEntry(i);
        if(jetR == R) { currentEvent = ievt; break; }
    }

    if(currentEvent < 0)
    {
        cout << "No jets found for R = " << R << endl;
        return;
    }

    for(Long64_t i=0;i<nEntries;i++)
    {
        jetTree->GetEntry(i);
        if(jetR != R) continue;

        if(ievt == currentEvent)
        {
            eventpT.push_back(jetpT);
            eventPhi.push_back(jetPhi);
            eventEta.push_back(jetEta);
        }
        else
        {
            processEvent(eventpT, eventPhi, eventEta, currentEvent);

            eventpT.clear(); eventPhi.clear(); eventEta.clear();
            currentEvent = ievt;
            eventpT.push_back(jetpT);
            eventPhi.push_back(jetPhi);
            eventEta.push_back(jetEta);
        }
    }

    // Dernier evenement
    processEvent(eventpT, eventPhi, eventEta, currentEvent);

    cout << "Evenements avec >=1 jet (R=" << R << ") : " << nEventsSeen << endl;
    cout << "Evenements retenus dans les plots dijet    : " << nEventsKept << endl;

    //==========================================================
    // Save in root file
    //==========================================================

    TFile *out = new TFile(Form("DijetDistributions_R%d.root",R),"RECREATE");
    hDeltaPhi->Write();
    hpT1->Write();
    hpT2->Write();
    hDeltapT->Write();
    hDeltaEta->Write();
    hDeltaR->Write();
    hNjets->Write();
    hDeltaEtaDeltaPhi->Write();
    hpT1pT2->Write();
    hPhiCorr->Write();
    hEtaCorr->Write();
    out->Close();

    //==========================================================
    // Save figures
    //==========================================================

    TString tag = Form("R%d_cut%d_cut2_%d_dR%.1f_2j%d", R, pTcut, pT2cut, dRcut, (int)require2jets);

    TCanvas *c1 = new TCanvas("c1","DeltaPhi",800,600);
    hDeltaPhi->Draw();
    c1->SaveAs(Form("DijetAnalysis/DeltaPhi_%s.pdf", tag.Data()));

    TCanvas *c2 = new TCanvas("c2","Leading jet pT",800,600);
    hpT1->Draw();
    c2->SaveAs(Form("DijetAnalysis/pT1_%s.pdf", tag.Data()));

    TCanvas *c3 = new TCanvas("c3","Subleading jet pT",800,600);
    hpT2->Draw();
    c3->SaveAs(Form("DijetAnalysis/pT2_%s.pdf", tag.Data()));

    TCanvas *c4 = new TCanvas("c4","Delta pT",800,600);
    hDeltapT->Draw();
    c4->SaveAs(Form("DijetAnalysis/DeltapT_%s.pdf", tag.Data()));

    TCanvas *c5 = new TCanvas("c5","Delta eta",800,600);
    hDeltaEta->Draw();
    c5->SaveAs(Form("DijetAnalysis/DeltaEta_%s.pdf", tag.Data()));

    TCanvas *c6 = new TCanvas("c6","Delta eta vs Delta phi",800,600);
    hDeltaEtaDeltaPhi->Draw("COLZ");
    c6->SaveAs(Form("DijetAnalysis/DeltaEtaDeltaPhi_%s.pdf", tag.Data()));

    TCanvas *c7 = new TCanvas("c7","Leading vs Subleading jet pT",800,600);
    hpT1pT2->Draw("COLZ");
    c7->SaveAs(Form("DijetAnalysis/pT1pT2_%s.pdf", tag.Data()));

    // (1) plots de correlation
    TCanvas *c8 = new TCanvas("c8","Phi1 vs Phi2",800,600);
    hPhiCorr->Draw("COLZ");
    c8->SaveAs(Form("DijetAnalysis/PhiCorr_%s.pdf", tag.Data()));

    TCanvas *c9 = new TCanvas("c9","Eta1 vs Eta2",800,600);
    hEtaCorr->Draw("COLZ");
    c9->SaveAs(Form("DijetAnalysis/EtaCorr_%s.pdf", tag.Data()));

    // (3) DeltaR leading-subleading
    TCanvas *c10 = new TCanvas("c10","DeltaR leading-subleading",800,600);
    hDeltaR->Draw();
    c10->SaveAs(Form("DijetAnalysis/DeltaR_%s.pdf", tag.Data()));

    // (4) multiplicite de jets par evenement
    TCanvas *c11 = new TCanvas("c11","Njets per event",800,600);
    hNjets->Draw();
    c11->SaveAs(Form("DijetAnalysis/Njets_R%d.pdf", R));

    cout << "Done." << endl;
}
