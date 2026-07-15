#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>

#include <vector>
#include <iostream>

using namespace std;

void AnalyzeDijetDeltaPhi(const char *filename="MergedAnalysisJets.root", const int R=2)
{
    //==========================================================
    // Open input file
    //==========================================================

    TFile *file = TFile::Open(filename);

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
        TMath::Pi()
    );

    TH1D *hPt1 = new TH1D(
        "hPt1",
        ";Leading jet p_{T} (GeV/c);Events",
        100,0,300
    );

    TH1D *hPt2 = new TH1D(
        "hPt2",
        ";Subleading jet p_{T} (GeV/c);Events",
        100,0,300
    );

    TH1D *hDeltaPt = new TH1D(
        "hDeltaPt",
        ";p_{T1}-p_{T2} (GeV/c);Events",
        100,0,200
    );

    TH1D *hDeltaEta = new TH1D(
        "hDeltaEta",
        ";#Delta#eta;Events",
        100,0,3
    );

    TH2D *hPt1Pt2 = new TH2D(
        "hPt1Pt2",
        ";Leading jet p_{T} (GeV/c);Subleading jet p_{T} (GeV/c)",
        100, 0, 300,
        100, 0, 300
    );

    TH2D *hDeltaEtaDeltaPhi = new TH2D(
        "hDeltaEtaDeltaPhi",
        ";#Delta#eta;#Delta#phi (rad)",
        100, 0, 3,
        64, 0, TMath::Pi()
    );

    //==========================================================
    // Variables
    //==========================================================

    vector<float> eventPt;
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
            eventPt.push_back(jetpT);
            eventPhi.push_back(jetPhi);
            eventEta.push_back(jetEta);
        }
        else
        {
            if(eventPt.size() >= 2)
            {
                int leading = -1;
                int subleading = -1;

                float pt1 = -1.;
                float pt2 = -1.;

                for(unsigned int j=0;j<eventPt.size();j++)
                {
                    if(eventPt[j] > pt1)
                    {
                        pt2 = pt1;
                        subleading = leading;

                        pt1 = eventPt[j];
                        leading = j;
                    }
                    else if(eventPt[j] > pt2)
                    {
                        pt2 = eventPt[j];
                        subleading = j;
                    }
                }

                if(leading >= 0 && subleading >= 0)
                {

                    // cut on leading jet pT
                    if(pt1 >= 100. && pt1 < 200.){

                    double dphi = fabs(eventPhi[leading]-eventPhi[subleading]);

                    if(dphi > TMath::Pi())
                        dphi = 2.*TMath::Pi()-dphi;

                    /*if(currentEvent < 10)
                    {
                        cout << "Event " << currentEvent
                             << "  Njets = " << eventPt.size()
                             << "  lead pT = " << pt1
                             << "  sublead pT = " << pt2
                             << "  phi1 = " << eventPhi[leading]
                             << "  phi2 = " << eventPhi[subleading]
                             << "  dphi = " << dphi
                             << endl;
                    }*/

                    hDeltaPhi->Fill(dphi);
                    hPt1->Fill(pt1);
                    hPt2->Fill(pt2);
                    hDeltaPt->Fill(pt1 - pt2);
                    double deta = fabs(eventEta[leading] - eventEta[subleading]);
                    hDeltaEta->Fill(deta);
                    hPt1Pt2->Fill(pt1, pt2);
                    hDeltaEtaDeltaPhi->Fill(deta, dphi);
                    }
                }
            }

            eventPt.clear();
            eventPhi.clear();
            eventEta.clear();
            currentEvent = ievt;
            eventPt.push_back(jetpT);
            eventPhi.push_back(jetPhi);
            eventEta.push_back(jetEta);
        }
    }

    //==========================================================
    // Last event
    //==========================================================

    if(eventPt.size() >= 2)
    {
        int leading = -1;
        int subleading = -1;

        float pt1 = -1.;
        float pt2 = -1.;

        for(unsigned int j=0;j<eventPt.size();j++)
        {
            if(eventPt[j] > pt1)
            {
                pt2 = pt1;
                subleading = leading;

                pt1 = eventPt[j];
                leading = j;
            }
            else if(eventPt[j] > pt2)
            {
                pt2 = eventPt[j];
                subleading = j;
            }
        }

        if(leading >= 0 && subleading >= 0)
        {
            // cut on leading jet pT
            if(pt1 >= 100. && pt1 < 200.)
            {

            double dphi = fabs(eventPhi[leading]-eventPhi[subleading]);

            if(dphi > TMath::Pi())
                dphi = 2.*TMath::Pi()-dphi;

            double deta = fabs(eventEta[leading] - eventEta[subleading]);
            /*if(currentEvent < 10)
            {
                cout << "Event " << currentEvent
                     << "  Njets = " << eventPt.size()
                     << "  lead pT = " << pt1
                     << "  sublead pT = " << pt2
                     << "  phi1 = " << eventPhi[leading]
                     << "  phi2 = " << eventPhi[subleading]
                     << "  dphi = " << dphi
                     << endl;
            }*/

            hDeltaPhi->Fill(dphi);
            hPt1->Fill(pt1);
            hPt2->Fill(pt2);
            hDeltaPt->Fill(pt1 - pt2);
            hDeltaEta->Fill(deta);
            hPt1Pt2->Fill(pt1, pt2);
            hDeltaEtaDeltaPhi->Fill(deta, dphi);
            }
        }
    }

    // Save in root file

    TFile *out = new TFile(Form("DijetDistributions_R%d.root",R),"RECREATE");

    hDeltaPhi->Write();
    hPt1->Write();
    hPt2->Write();
    hDeltaPt->Write();
    hDeltaEta->Write();
    hDeltaEtaDeltaPhi->Write();
    hPt1Pt2->Write();
    out->Close();

    // Save figures
    // Delta Phi distribution
    TCanvas *c1 = new TCanvas("c1","DeltaPhi",800,600);
    hDeltaPhi->Draw();
    c1->SaveAs(Form("DeltaPhi_R%d.pdf", R));

    // Leading jet pT distribution
    TCanvas *c2 = new TCanvas("c2","Leading jet pT",800,600);
    hPt1->Draw();
    c2->SaveAs(Form("Pt1_R%d.pdf", R));

    // Subleading jet pT distribution
    TCanvas *c3 = new TCanvas("c3","Subleading jet pT",800,600);
    hPt2->Draw();
    c3->SaveAs(Form("Pt2_R%d.pdf", R));

    // Delta pT distribution
    TCanvas *c4 = new TCanvas("c4","Delta pT",800,600);
    hDeltaPt->Draw();
    c4->SaveAs(Form("DeltaPt_R%d.pdf", R));

    // Delta eta distribution
    TCanvas *c5 = new TCanvas("c5","Delta eta",800,600);
    hDeltaEta->Draw();
    c5->SaveAs(Form("DeltaEta_R%d.pdf", R));

    TCanvas *c6 = new TCanvas("c6","Delta eta vs Delta phi",800,600);
    hDeltaEtaDeltaPhi->Draw("COLZ");
    c6->SaveAs(Form("DeltaEtaDeltaPhi_R%d.pdf", R));

    TCanvas *c7 = new TCanvas("c7","Leading vs Subleading jet pT",800,600);
    hPt1Pt2->Draw("COLZ");
    c7->SaveAs(Form("Pt1Pt2_R%d.pdf", R));

    cout << "Done." << endl;
}