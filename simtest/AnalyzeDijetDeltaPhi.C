#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
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

    //==========================================================
    // Variables
    //==========================================================

    vector<float> eventPt;
    vector<float> eventPhi;

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
                }
            }

            eventPt.clear();
            eventPhi.clear();

            currentEvent = ievt;
            eventPt.push_back(jetpT);
            eventPhi.push_back(jetPhi);
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
        }
    }

    // Save output

    TFile *out = new TFile(Form("DijetDeltaPhi_R%d.root", R), "RECREATE");
    hDeltaPhi->Write();
    out->Close();

    // Draw

    TCanvas *c = new TCanvas("c","DeltaPhi",800,600);

    hDeltaPhi->Draw();

    c->SaveAs(Form("DeltaPhi_R%d.pdf", R));

    cout << "Done." << endl;
}