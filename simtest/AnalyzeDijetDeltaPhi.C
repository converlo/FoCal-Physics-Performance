#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>

#include <vector>
#include <iostream>

using namespace std;

void AnalyzeDijetDeltaPhi(const char *filename="MergedAnalysisJets.root")
{

    // Open the input file

    TFile *file = TFile::Open(filename);

    // JetTree

    TTree *jetTree = (TTree*)file->Get("jetTree");

    if (!jetTree) {
        cout << "jetTree not found." << endl;
        return;
    }

    // Branch

    Int_t ievt;
    Float_t jetpT;
    Float_t jetPhi;
    Float_t jetEta;

    jetTree->SetBranchAddress("ievt",&ievt);
    jetTree->SetBranchAddress("jetpT",&jetpT);
    jetTree->SetBranchAddress("jetPhi",&jetPhi);
    jetTree->SetBranchAddress("jetEta",&jetEta);

    // Histogram

    TH1D *hDeltaPhi = new TH1D(
        "hDeltaPhi",
        ";#Delta#phi (rad);Events",
        64,
        0,
        TMath::Pi()
    );

    // Variables describing one event

    vector<float> eventPt;
    vector<float> eventPhi;

    Long64_t nEntries = jetTree->GetEntries();

    cout << "Number of jets = " << nEntries << endl;

    jetTree->GetEntry(0);

    Int_t currentEvent = ievt;

    // Loop over all jets

    for(Long64_t i=0;i<nEntries;i++)
    {
        jetTree->GetEntry(i);

        //------------------------------------------------------
        // Same event
        //------------------------------------------------------

        if(ievt == currentEvent)
        {
            eventPt.push_back(jetpT);
            eventPhi.push_back(jetPhi);
        }

        //------------------------------------------------------
        // New event
        //------------------------------------------------------

        else
        {
            //--------------------------------------------------
            // Analyse previous event
            //--------------------------------------------------

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

                if(leading >=0 && subleading >=0)
                {
                    double dphi = fabs(
                        eventPhi[leading]
                        -
                        eventPhi[subleading]
                    );

                    if(dphi > TMath::Pi())
                        dphi = 2.*TMath::Pi()-dphi;

                    hDeltaPhi->Fill(dphi);
                }
            }

            //--------------------------------------------------
            // Start new event
            //--------------------------------------------------

            eventPt.clear();
            eventPhi.clear();

            currentEvent = ievt;

            eventPt.push_back(jetpT);
            eventPhi.push_back(jetPhi);
        }
    }

    //==========================================================
    // Don't forget the last event
    //==========================================================

    if(eventPt.size() >= 2)
    {
        int leading=-1;
        int subleading=-1;

        float pt1=-1.;
        float pt2=-1.;

        for(unsigned int j=0;j<eventPt.size();j++)
        {
            if(eventPt[j] > pt1)
            {
                pt2=pt1;
                subleading=leading;

                pt1=eventPt[j];
                leading=j;
            }
            else if(eventPt[j]>pt2)
            {
                pt2=eventPt[j];
                subleading=j;
            }
        }

        if(leading>=0 && subleading>=0)
        {
            double dphi=fabs(
                eventPhi[leading]
                -
                eventPhi[subleading]
            );

            if(dphi>TMath::Pi())
                dphi=2.*TMath::Pi()-dphi;

            hDeltaPhi->Fill(dphi);
        }
    }

    //==========================================================
    // Save output
    //==========================================================



    TFile *out = new TFile("DijetDeltaPhi.root","RECREATE");

    if (currentEvent <10){
        cout << "Event " << currentEvent
        << "  lead pt = " << pt1
        << "  sublead pt = " << pt2
        << "  dphi = " << dphi
        << endl;
    }

    hDeltaPhi->Write();

    out->Close();

    //==========================================================
    // Draw
    //==========================================================

    TCanvas *c = new TCanvas("c","DeltaPhi",800,600);

    hDeltaPhi->Draw();

    c->SaveAs("DeltaPhi.png");

    cout << "Done." << endl;
}