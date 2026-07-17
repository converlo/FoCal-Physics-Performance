#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TSystem.h>
#include <vector>
#include <iostream>

using namespace std;

void AnalyzeDijetDeltaPhi(const char *filename="MergedAnalysisJets.root", const int R=2, const int pTcut=10)
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
    // Histogram
    //==========================================================

    TH1D *hDeltaPhi = new TH1D(
        "hDeltaPhi",
        ";#Delta#phi (rad);Events",
        64,
        0,
        2*TMath::Pi()
    );

    TH1D *hpT1 = new TH1D(
        "hpT1",
        ";Leading jet p_{T} (GeV/c);Events",
        100,0,300
    );

    TH1D *hpT2 = new TH1D(
        "hpT2",
        ";Subleading jet p_{T} (GeV/c);Events",
        100,0,300
    );

    TH1D *hDeltapT = new TH1D(
        "hDeltapT",
        ";p_{T1}-p_{T2} (GeV/c);Events",
        100,0,200
    );

    TH1D *hDeltaEta = new TH1D(
        "hDeltaEta",
        ";#Delta#eta;Events",
        100,0,3
    );

    TH2D *hpT1pT2 = new TH2D(
        "hpT1pT2",
        ";Leading jet p_{T} (GeV/c);Subleading jet p_{T} (GeV/c)",
        100, 0, 300,
        100, 0, 300
    );

    TH2D *hDeltaEtaDeltaPhi = new TH2D(
        "hDeltaEtaDeltaPhi",
        ";#Delta#eta;#Delta#phi (rad)",
        100, 0, 3,
        64, 0, 2*TMath::Pi()
    );

    //==========================================================
    // Variables
    //==========================================================

    vector<float> eventpT;
    vector<float> eventPhi;
    vector<float> eventEta;

    Long64_t nEntries = jetTree->GetEntries();

    cout << "Number of jets = " << nEntries << endl;


    Int_t currentEvent = -1;

    // Cherche le premier événement ayant le rayon demandé
    for(Long64_t i=0; i<nEntries; i++)
    {
        jetTree->GetEntry(i);

        if(jetR == R)
        {
            currentEvent = ievt;
            break;
        }
    }

    if(currentEvent < 0)
    {
        cout << "No jets found for R = " << R << endl;
        return;
    }
    //==========================================================
    // Loop over all jets
    //==========================================================

    for(Long64_t i=0;i<nEntries;i++)
    {
        jetTree->GetEntry(i);

        if(jetR != R)
            continue;
            
        if(ievt == currentEvent)
        {
            eventpT.push_back(jetpT);
            eventPhi.push_back(jetPhi);
            eventEta.push_back(jetEta);
        }
        else
        {
            if(eventpT.size() >= 2)
            {
                int leading = -1;
                int subleading = -1;

                float pt1 = -1.;
                float pt2 = -1.;

                for(unsigned int j=0;j<eventpT.size();j++)
                {
                    if(eventpT[j] > pt1)
                    {
                        pt2 = pt1;
                        subleading = leading;

                        pt1 = eventpT[j];
                        leading = j;
                    }
                    else if(eventpT[j] > pt2)
                    {
                        pt2 = eventpT[j];
                        subleading = j;
                    }
                }

                if(leading >= 0 && subleading >= 0)
                {

                    // cut on leading jet pT
                    if(pt1 >= pTcut){

                    double dphi = fabs(eventPhi[leading]-eventPhi[subleading]);

                    if(currentEvent < 10)
                    {
                        cout << "Event " << currentEvent
                             << "  Njets = " << eventpT.size()
                             << "  lead pT = " << pt1
                             << "  sublead pT = " << pt2
                             << "  phi1 = " << eventPhi[leading]
                             << "  phi2 = " << eventPhi[subleading]
                             << "  dphi = " << dphi
                             << endl;
                    }

                    hDeltaPhi->Fill(dphi);
                    hpT1->Fill(pt1);
                    hpT2->Fill(pt2);
                    hDeltapT->Fill(pt1 - pt2);
                    double deta = fabs(eventEta[leading] - eventEta[subleading]);
                    hDeltaEta->Fill(deta);
                    hpT1pT2->Fill(pt1, pt2);
                    hDeltaEtaDeltaPhi->Fill(deta, dphi);
                    }
                }
            }

            eventpT.clear();
            eventPhi.clear();
            eventEta.clear();
            currentEvent = ievt;
            eventpT.push_back(jetpT);
            eventPhi.push_back(jetPhi);
            eventEta.push_back(jetEta);
        }
    }

    //==========================================================
    // Last event
    //==========================================================

    if(eventpT.size() >= 2)
    {
        int leading = -1;
        int subleading = -1;

        float pt1 = -1.;
        float pt2 = -1.;

        for(unsigned int j=0;j<eventpT.size();j++)
        {
            if(eventpT[j] > pt1)
            {
                pt2 = pt1;
                subleading = leading;

                pt1 = eventpT[j];
                leading = j;
            }
            else if(eventpT[j] > pt2)
            {
                pt2 = eventpT[j];
                subleading = j;
            }
        }

        if(leading >= 0 && subleading >= 0)
        {
            // cut on leading jet pT
            if(pt1 >= pTcut)
            {

            double dphi = fabs(eventPhi[leading]-eventPhi[subleading]);

            double deta = fabs(eventEta[leading] - eventEta[subleading]);
            /*if(currentEvent < 10)
            {
                cout << "Event " << currentEvent
                     << "  Njets = " << eventpT.size()
                     << "  lead pT = " << pt1
                     << "  sublead pT = " << pt2
                     << "  phi1 = " << eventPhi[leading]
                     << "  phi2 = " << eventPhi[subleading]
                     << "  dphi = " << dphi
                     << endl;
            }*/

            hDeltaPhi->Fill(dphi);
            hpT1->Fill(pt1);
            hpT2->Fill(pt2);
            hDeltapT->Fill(pt1 - pt2);
            hDeltaEta->Fill(deta);
            hpT1pT2->Fill(pt1, pt2);
            hDeltaEtaDeltaPhi->Fill(deta, dphi);
            }
        }
    }

    // Save in root file

    TFile *out = new TFile(Form("DijetDistributions_R%d.root",R),"RECREATE");

    hDeltaPhi->Write();
    hpT1->Write();
    hpT2->Write();
    hDeltapT->Write();
    hDeltaEta->Write();
    hDeltaEtaDeltaPhi->Write();
    hpT1pT2->Write();
    out->Close();

    // Save figures
    // Delta Phi distribution
    TCanvas *c1 = new TCanvas("c1","DeltaPhi",800,600);
    hDeltaPhi->Draw();
    c1->SaveAs(Form("DijetAnalysis/DeltaPhi_R%d_cut%d.pdf", R, pTcut));

    // Leading jet pT distribution
    TCanvas *c2 = new TCanvas("c2","Leading jet pT",800,600);
    hpT1->Draw();
    c2->SaveAs(Form("DijetAnalysis/pT1_R%d_cut%d.pdf", R, pTcut));

    // Subleading jet pT distribution
    TCanvas *c3 = new TCanvas("c3","Subleading jet pT",800,600);
    hpT2->Draw();
    c3->SaveAs(Form("DijetAnalysis/pT2_R%d_cut%d.pdf", R, pTcut));

    // Delta pT distribution
    TCanvas *c4 = new TCanvas("c4","Delta pT",800,600);
    hDeltapT->Draw();
    c4->SaveAs(Form("DijetAnalysis/DeltapT_R%d_cut%d.pdf", R, pTcut));

    // Delta eta distribution
    TCanvas *c5 = new TCanvas("c5","Delta eta",800,600);
    hDeltaEta->Draw();
    c5->SaveAs(Form("DijetAnalysis/DeltaEta_R%d_cut%d.pdf", R, pTcut));

    TCanvas *c6 = new TCanvas("c6","Delta eta vs Delta phi",800,600);
    hDeltaEtaDeltaPhi->Draw("COLZ");
    c6->SaveAs(Form("DijetAnalysis/DeltaEtaDeltaPhi_R%d_cut%d.pdf", R, pTcut));

    TCanvas *c7 = new TCanvas("c7","Leading vs Subleading jet pT",800,600);
    hpT1pT2->Draw("COLZ");
    c7->SaveAs(Form("DijetAnalysis/pT1pT2_R%d_cut%d.pdf", R, pTcut));

    cout << "Done." << endl;
}
