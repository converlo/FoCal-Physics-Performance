/* Macro for plotting some histograms from jet trees
    2.21064
    0.0669805
    0.00182628
    0.000139462
    0.0000225822
*/
#include "iostream"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TF1.h"

double GetMedian(TH1D *h);

void JetPlottingMergedFile(int Rvalue = 2)
{
    const int nR = 3;                        // 5
    const float Rvals[nR] = {0.2, 0.4, 0.6}; // Cone radii

    const int nCol = 10;
    const int gcolors[nCol] = {1, 2, 6, 4, 7, 1, 2, 4, 6, 7};
    const int gmarkers[nCol] = {4, 8, 25, 21, 8, 21, 25, 4, 8, 21};

    const float etaMin = 3.4;
    const float etaMax = 5.5;

    const int nPtBins = 14;
    const int nEBins = 15;
    // const double JetPtBorders[nPtBins] = {5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 400.0}; // nPtBins = 10
    const double JetPtBorders[nPtBins] = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 150.0};                         // nPtBins = 14
    const double JetEBorders[nEBins] = {100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 3000.0}; // const int nEBins  = 15;
    double medianspT[nPtBins - 1];
    double meanspT[nPtBins - 1];
    double SDpT[nPtBins - 1];
    double mediansE[nEBins - 1];
    double meansE[nEBins - 1];
    double SDE[nEBins - 1];

    double gausMeanspT[nPtBins - 1];
    double gausMeansE[nEBins - 1];
    double gausSDpT[nPtBins - 1];
    double gausSDE[nEBins - 1];
    // eta binned
    const int nEtaBins = 3; // 11;
    const float EtaBinBorders[nEtaBins] = {3.8, 4.5, 5.1};

    double EtamedianspT[nEtaBins - 1][nPtBins - 1];
    double EtameanspT[nEtaBins - 1][nPtBins - 1];
    double EtaSDpT[nEtaBins - 1][nPtBins - 1];

    double EtamediansE[nEtaBins - 1][nEBins - 1];
    double EtameansE[nEtaBins - 1][nEBins - 1];
    double EtaSDE[nEtaBins - 1][nEBins - 1];

    TFile *jetFile = TFile::Open(Form("Data20230816/JES/EnMergedR%d_pT_0Mass.root", int(Rvals[Rvalue] * 10)));

    TTree *jetTree = (TTree *)jetFile->Get("jetTree");
    TTree *TruthjetTree = (TTree *)jetFile->Get("truthjetTree");
    // TTree *PDGTree = (TTree *)jetFile->Get("inPDGTree");

    TFile *fout = new TFile(Form("Data20230816/JES/EneMerged_OutputR%d_pT_0Mass.root", int(Rvals[Rvalue] * 10)), "RECREATE");

    // histograms

    TH2D *hRespMatrix_pT_Eta[nEtaBins - 1];
    TH2D *hRespMatrix_E_Eta[nEtaBins - 1];
    TH1D *hjetRatiopT_Eta[nEtaBins - 1][nPtBins - 1];
    TH1D *hjetRatioE_Eta[nEtaBins - 1][nEBins - 1];
    TH1D *hEtaMedianpT[nEtaBins - 1], *hEtaMedianE[nEtaBins - 1];
    TH1D *hEtaMeanpT[nEtaBins - 1], *hEtaMeanE[nEtaBins - 1];
    TH1D *hEtaSDpT[nEtaBins - 1], *hEtaSDE[nEtaBins - 1];

    // Not-eta-binned histograms
    TH1F *hjetpT[nR], *hjetE[nR];
    TH1F *hTruthjetpT[nR], *hTruthjetE[nR];
    TH1F *hTruthjetpT_Matched[nR], *hTruthjetE_Matched[nR];

    TH1D *hjetRatiopT[nPtBins - 1], *hjetRatioE[nEBins - 1];
    TH1D *hMedianpT, *hMedianE;
    TH1D *hMeanpT, *hMeanE;
    TH1D *hSDpT, *hSDE;
    TH1D *hgausMeanspT, *hgausMeansE;
    TH1D *hgausSDpT, *hgausSDE;

    hMedianpT = new TH1D("hMedianpT", "Mean and median of #Deltap_{T} distribution;p_{T}^{part} (GeV/c);Mean or median", nPtBins - 1, JetPtBorders);
    hMedianE = new TH1D("hMedianE", "Mean and median of #DeltaE distribution;p_{T}^{part} (GeV/c);Mean or median", nEBins - 1, JetEBorders);
    hMeanpT = new TH1D("hMeanpT", "Mean and median of #Deltap_{T} distribution;p_{T}^{part} (GeV/c);Standard deviation", nPtBins - 1, JetPtBorders);
    hMeanE = new TH1D("hMeanE", "Mean and median of #DeltaE distribution;E^{part} (GeV);Mean or median", nEBins - 1, JetEBorders);
    hSDpT = new TH1D("hSDpT", "Standard deviation of #Deltap_{T} distribution;E^{part} (GeV);Mean or median", nPtBins - 1, JetPtBorders);
    hSDE = new TH1D("hSDE", "Standard deviation of #DeltaE distribution;E^{part} (GeV);Standard deviation", nEBins - 1, JetEBorders);

    hgausMeanspT = new TH1D("hgausMeanspT", "Mean and median of #Deltap_{T} distribution", nPtBins - 1, JetPtBorders);
    hgausMeansE = new TH1D("hgausMeansE", "Mean and median of #DeltaE distribution", nEBins - 1, JetEBorders);
    hgausSDpT = new TH1D("hgausSDpT", "Standard deviation of #Deltap_{T} distribution", nPtBins - 1, JetPtBorders);
    hgausSDE = new TH1D("hgausSDE", "Standard deviation of #DeltaE distribution", nEBins - 1, JetEBorders);

    for (int iE = 0; iE < nEtaBins - 1; ++iE)
    {
        // hRespMatrix_pT_Eta[iE] = new TH2D(Form("hRespMatrix_pT_Eta_%d", iE), Form("Detector response matrix, eta %0.1f - %0.1f, R=%0.1f;#it{p}_{T}^{det} (GeV/c);#it{p}_{T}^{part} (GeV/c)", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]), 1000, 0, 1000, 1000, 0, 1000);
        // hRespMatrix_E_Eta[iE] = new TH2D(Form("hRespMatrix_E_Eta_%d", iE), Form("Detector response matrix, eta %0.1f - %0.1f, R=%0.1f;#it{E}_{jet}^{det} (GeV);#it{E}_{jet}^{part} (GeV)", EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]), 1000, 0, 1000, 1000, 0, 1000);

        hRespMatrix_pT_Eta[iE] = (TH2D *)((TH2D *)jetFile->Get(Form("hRespMatrix_pT_Eta_%d", iE)))->Clone();

        hRespMatrix_E_Eta[iE] = (TH2D *)((TH2D *)jetFile->Get(Form("hRespMatrix_E_Eta_%d", iE)))->Clone();

        hEtaMedianpT[iE] = new TH1D(Form("hEtaMedianpT_%d", iE), Form("Mean and median of #Deltap_{T} distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE + 1], Rvals[Rvalue]), nPtBins - 1, JetPtBorders);
        hEtaMedianE[iE] = new TH1D(Form("hEtaMedianE_%d", iE), Form("Mean and median of #DeltaE distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE + 1], Rvals[Rvalue]), nEBins - 1, JetEBorders);
        hEtaMeanpT[iE] = new TH1D(Form("hEtaMeanpT_%d", iE), Form("Mean and median of #Deltap_{T} distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE + 1], Rvals[Rvalue]), nPtBins - 1, JetPtBorders);
        hEtaMeanE[iE] = new TH1D(Form("hEtaMeanE_%d", iE), Form("Mean and median of #DeltaE distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE + 1], Rvals[Rvalue]), nEBins - 1, JetEBorders);
        hEtaSDpT[iE] = new TH1D(Form("hEtaSDpT_%d", iE), Form("Standard deviation of #Deltap_{T} distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE + 1], Rvals[Rvalue]), nPtBins - 1, JetPtBorders);
        hEtaSDE[iE] = new TH1D(Form("hEtaSDE_%d", iE), Form("Standard deviation of #DeltaE distribution, eta %0.1f - %0.1f, R=%0.1f", EtaBinBorders[iE], EtaBinBorders[iE + 1], Rvals[Rvalue]), nEBins - 1, JetEBorders);
    }

    // TH2D *hRespMatrix_pT = new TH2D("hRespMatrix_pT", Form("Detector response matrix, R=%0.1f;#it{p}_{T}^{det} (GeV/c);#it{p}_{T}^{part} (GeV/c)", Rvals[Rvalue]), 1000, 0, 1000, 1000, 0, 1000);
    // TH2D *hRespMatrix_E = new TH2D("hRespMatrix_E", Form("Detector response matrix, R=%0.1f;#it{E}_{jet}^{det} (GeV);#it{E}_{jet}^{part} (GeV)", Rvals[Rvalue]), 1000, 0, 1000, 1000, 0, 1000);

    TH2D *hRespMatrix_pT = (TH2D *)((TH2D *)jetFile->Get("hRespMatrix_pT"))->Clone();
    TH2D *hRespMatrix_E = (TH2D *)((TH2D *)jetFile->Get("hRespMatrix_E"))->Clone();

    TH1D *hJetRecEfficiency_pT[nR];
    TH1D *hJetRecEfficiency_E[nR];

    for (int iR = 0; iR < nR; iR++)
    {
        hjetpT[iR] = (TH1F *)((TH1F *)jetFile->Get(Form("hjetpT_%d", iR)))->Clone();
        hjetE[iR] = (TH1F *)((TH1F *)jetFile->Get(Form("hjetE_%d", iR)))->Clone();
        hTruthjetpT_Matched[iR] = (TH1F *)((TH1F *)jetFile->Get(Form("hTruthjetpT_Matched_%d", iR)))->Clone();
        hTruthjetE_Matched[iR] = (TH1F *)((TH1F *)jetFile->Get(Form("hTruthjetE_Matched_%d", iR)))->Clone();
        hTruthjetpT[iR] = (TH1F *)((TH1F *)jetFile->Get(Form("hTruthjetpT_%d", iR)))->Clone();
        hTruthjetE[iR] = (TH1F *)((TH1F *)jetFile->Get(Form("hTruthjetE_%d", iR)))->Clone();

        hJetRecEfficiency_pT[iR] = (TH1D *)hTruthjetpT_Matched[iR]->Clone();
        hJetRecEfficiency_pT[iR] = (TH1D *)hJetRecEfficiency_pT[iR]->Rebin(nPtBins - 1, Form("hJetRecEfficiency_pT_%d", iR), JetPtBorders);
        hJetRecEfficiency_pT[iR]->Divide((TH1D *)hTruthjetpT[iR]->Rebin(nPtBins - 1, "", JetPtBorders));
        hJetRecEfficiency_E[iR] = (TH1D *)hTruthjetE_Matched[iR]->Clone();
        hJetRecEfficiency_E[iR] = (TH1D *)hJetRecEfficiency_E[iR]->Rebin(nEBins - 1, Form("hJetRecEfficiency_E_%d", iR), JetEBorders);
        hJetRecEfficiency_E[iR]->Divide((TH1D *)hTruthjetE[iR]->Rebin(nEBins - 1, "", JetEBorders));
    }

    // pt loop
    for (int ipt = 0; ipt < nPtBins - 1; ++ipt)
    {

        // TF1 *f1 = new TF1("f1","gaus",-1,0.6);

        // hjetRatiopT[ipt] = new TH1D(Form("hjetRatiopT_%d", ipt), Form("Jet-by-jet #Deltap_{T} distribution, p_{T}: %d - %d GeV/c", int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1])), 50, -1.0, 1.0);
        // hjetRatioE[ipt] = new TH1D(Form("hjetRatioE_%d", ipt), Form("Jet-by-jet #Delta E distribution, p_{T}: %d - %d GeV/c", int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1])), 50, -1.0, 1.0);

        // Fill deltapt histos, no eta cut
        // Add here filling from th2d
        // jetTree->Draw(Form("(jetpT-jetpT_match)/jetpT_match>>hjetRatiopT_%d", ipt), Form("jetR==%d && jetpT_match>=%d  && jetpT_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f", int(Rvals[Rvalue] * 10), int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]),Rvals[Rvalue],Rvals[Rvalue],Rvals[Rvalue],Rvals[Rvalue], Rvals[Rvalue]*0.6), "goff");
        // hjetRatiopT[ipt]= hRespMatrix_pT->ProjectionX(Form("hjetRatiopT_%d", ipt), int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]), "e");
        hjetRatiopT[ipt] = (TH1D *)((TH1D *)jetFile->Get(Form("hjetRatiopT_%d", ipt)))->Clone();
        // hjetRatiopT[ipt]->Fit("f1");

        medianspT[ipt] = GetMedian(hjetRatiopT[ipt]);
        meanspT[ipt] = hjetRatiopT[ipt]->GetMean();
        SDpT[ipt] = hjetRatiopT[ipt]->GetStdDev();

        // medianspT[ipt]=f1->GetParameter(1);
        // meanspT[ipt]=f1->GetParameter(1);
        // SDpT[ipt]=f1->GetParameter(2);

        hMedianpT->Fill(JetPtBorders[ipt + 1] - ((JetPtBorders[ipt + 1] - JetPtBorders[ipt]) / 2.0), medianspT[ipt]);
        hMedianpT->SetBinError(ipt + 1, hjetRatiopT[ipt]->GetMeanError());

        hMeanpT->Fill(JetPtBorders[ipt + 1] - ((JetPtBorders[ipt + 1] - JetPtBorders[ipt]) / 2.0), meanspT[ipt]);
        hMeanpT->SetBinError(ipt + 1, hjetRatiopT[ipt]->GetMeanError());
        hSDpT->Fill(JetPtBorders[ipt + 1] - ((JetPtBorders[ipt + 1] - JetPtBorders[ipt]) / 2.0), SDpT[ipt]);
        hSDpT->SetBinError(ipt + 1, hjetRatiopT[ipt]->GetStdDevError());

        // SCALE?
        // hjetRatiopT[ipt]->Scale(1., "width");
        hjetRatiopT[ipt]->Scale(1. / hjetRatiopT[ipt]->Integral(), "");
        // hjetRatiopT[ipt]->Scale(1./hjetRatiopT[ipt]->GetEntries());
        std::cout << "In jet pT range " << JetPtBorders[ipt] << " - " << JetPtBorders[ipt + 1] << "GeV/c, number of jets accepted = " << hjetRatiopT[ipt]->GetEntries() << std::endl;

        hjetRatiopT[ipt]->Fit("gaus");
        gausMeanspT[ipt] = hjetRatiopT[ipt]->GetFunction("gaus")->GetParameter(1);
        gausSDpT[ipt] = hjetRatiopT[ipt]->GetFunction("gaus")->GetParameter(2);

        hgausMeanspT->Fill(JetPtBorders[ipt + 1] - ((JetPtBorders[ipt + 1] - JetPtBorders[ipt]) / 2.0), gausMeanspT[ipt]);
        hgausMeanspT->SetBinError(ipt + 1, hjetRatiopT[ipt]->GetFunction("gaus")->GetParError(1));
        hgausSDpT->Fill(JetPtBorders[ipt + 1] - ((JetPtBorders[ipt + 1] - JetPtBorders[ipt]) / 2.0), gausSDpT[ipt]);
        hgausSDpT->SetBinError(ipt + 1, hjetRatiopT[ipt]->GetFunction("gaus")->GetParError(2));

        // same for eta cut
        for (int iE = 0; iE < nEtaBins - 1; ++iE)
        {
            // hjetRatiopT_Eta[iE][ipt] = new TH1D(Form("hjetRatiopT_Eta_%d_%d", iE, ipt), Form("Jet-by-jet #Deltap_{T} distribution, #eta: %0.1f - %0.1f (det jet), p_{T}: %d - %d GeV/c", EtaBinBorders[iE], EtaBinBorders[iE+1], int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1])), 50, -1.0, 1.0);
            // jetTree->Draw(Form("(jetpT-jetpT_match)/jetpT_match>>hjetRatiopT_Eta_%d_%d", iE, ipt), Form("jetR==%d && jetpT_match>=%d && jetpT_match<%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f", int(Rvals[Rvalue] * 10), int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]),EtaBinBorders[iE], EtaBinBorders[iE+1],EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]*0.6), "goff");
            hjetRatiopT_Eta[iE][ipt] = (TH1D *)((TH1D *)jetFile->Get(Form("hjetRatiopT_Eta_%d_%d", iE, ipt)))->Clone();
            // hjetRatiopT_Eta[iE][ipt]->Fit("f1");

            EtamedianspT[iE][ipt] = GetMedian(hjetRatiopT_Eta[iE][ipt]);
            EtameanspT[iE][ipt] = hjetRatiopT_Eta[iE][ipt]->GetMean();
            EtaSDpT[iE][ipt] = hjetRatiopT_Eta[iE][ipt]->GetStdDev();
            // EtamedianspT[iE][ipt]=f1->GetParameter(1);
            // EtameanspT[iE][ipt]=f1->GetParameter(1);
            // EtaSDpT[iE][ipt]=f1->GetParameter(2);

            hEtaMedianpT[iE]->Fill(JetPtBorders[ipt + 1] - ((JetPtBorders[ipt + 1] - JetPtBorders[ipt]) / 2.0), EtamedianspT[iE][ipt]);
            hEtaMedianpT[iE]->SetBinError(ipt + 1, hjetRatiopT_Eta[iE][ipt]->GetMeanError());

            hEtaMeanpT[iE]->Fill(JetPtBorders[ipt + 1] - ((JetPtBorders[ipt + 1] - JetPtBorders[ipt]) / 2.0), EtameanspT[iE][ipt]);
            hEtaMeanpT[iE]->SetBinError(ipt + 1, hjetRatiopT_Eta[iE][ipt]->GetMeanError());
            hEtaSDpT[iE]->Fill(JetPtBorders[ipt + 1] - ((JetPtBorders[ipt + 1] - JetPtBorders[ipt]) / 2.0), EtaSDpT[iE][ipt]);
            hEtaSDpT[iE]->SetBinError(ipt + 1, hjetRatiopT_Eta[iE][ipt]->GetStdDevError());

            //            hjetRatiopT_Eta[iE][ipt]->Scale(1., "width");
            hjetRatiopT_Eta[iE][ipt]->Scale(1. / hjetRatiopT_Eta[iE][ipt]->Integral(), "");
            // hjetRatiopT_Eta[iE][ipt]->Scale(1./hjetRatiopT_Eta[iE][ipt]->GetEntries());
            // cout << "In jet pT range " << JetPtBorders[ipt] << " - " << JetPtBorders[ipt+1] << "GeV/c, number of jets accepted = " << hjetRatiopT_Eta[iE][ipt]->GetEntries() << endl;
        }
    }

    // E loop
    for (int ipt = 0; ipt < nEBins - 1; ++ipt)
    {

        // TF1 *f1 = new TF1("f1","gaus",-1,0.6);

        // Fill deltaE histos, no eta cut
        // jetTree->Draw(Form("(jetE-jetE_match)/jetE_match>>hjetRatioE_%d", ipt), Form("jetR==%d && jetpT_match>=%d  && jetpT_match<%d  && jetEta>3.4+%f && jetEta<5.5-%f && jetEta_match>3.4+%f && jetEta_match<5.5-%f && jet_distmatch<%f", int(Rvals[Rvalue] * 10), int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]),Rvals[Rvalue],Rvals[Rvalue],Rvals[Rvalue],Rvals[Rvalue], Rvals[Rvalue]*0.6), "goff");
        // hjetRatioE[ipt]= hRespMatrix_E->ProjectionX(Form("hjetRatioE_%d", ipt), int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]), "e");

        hjetRatioE[ipt] = (TH1D *)((TH1D *)jetFile->Get(Form("hjetRatioE_%d", ipt)))->Clone();
        // hjetRatioE[ipt]->Fit("f1");

        mediansE[ipt] = GetMedian(hjetRatioE[ipt]);
        meansE[ipt] = hjetRatioE[ipt]->GetMean();
        SDE[ipt] = hjetRatioE[ipt]->GetStdDev();
        // mediansE[ipt]=f1->GetParameter(1);
        // meansE[ipt]=f1->GetParameter(1);
        // SDE[ipt]=f1->GetParameter(2);

        hMedianE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), mediansE[ipt]);
        hMedianE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetMeanError());
        hMeanE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), meansE[ipt]);
        hMeanE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetMeanError());
        hSDE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2), SDE[ipt]);
        hSDE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetStdDevError());

        // SCALE?
        // hjetRatioE[ipt]->Scale(1., "width");
        hjetRatioE[ipt]->Scale(1 / hjetRatioE[ipt]->Integral(), "");
        // hjetRatioE[ipt]->Scale(1/hjetRatioE[ipt]->GetEntries());

        hjetRatioE[ipt]->Fit("gaus");
        gausMeansE[ipt] = hjetRatioE[ipt]->GetFunction("gaus")->GetParameter(1);
        gausSDE[ipt] = hjetRatioE[ipt]->GetFunction("gaus")->GetParameter(2);
        hgausMeansE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2.0), gausMeansE[ipt]);
        hgausMeansE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetFunction("gaus")->GetParError(1));
        hgausSDE->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2.0), gausSDE[ipt]);
        hgausSDE->SetBinError(ipt + 1, hjetRatioE[ipt]->GetFunction("gaus")->GetParError(2));

        hjetRatioE[ipt]->GetXaxis()->SetTitle("(E^{det}-E^{part})/E^{part}");
        hjetRatioE[ipt]->GetYaxis()->SetTitle("probability");

        // same for eta cut
        for (int iE = 0; iE < nEtaBins - 1; ++iE)
        {
            // hjetRatioE_Eta[iE][ipt] = new TH1D(Form("hjetRatioE_Eta_%d_%d", iE, ipt), Form("Jet-by-jet #Deltap_{T} distribution, #eta: %0.1f - %0.1f (det jet), p_{T}: %d - %d GeV/c", EtaBinBorders[iE], EtaBinBorders[iE+1], int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1])), 50, -1.0, 1.0);
            // jetTree->Draw(Form("(jetpT-jetpT_match)/jetpT_match>>hjetRatioE_Eta_%d_%d", iE, ipt), Form("jetR==%d && jetpT_match>=%d && jetpT_match<%d && jetEta>%f && jetEta<%f && jetEta_match>%f && jetEta_match<%f && jet_distmatch<%f", int(Rvals[Rvalue] * 10), int(JetPtBorders[ipt]), int(JetPtBorders[ipt+1]),EtaBinBorders[iE], EtaBinBorders[iE+1],EtaBinBorders[iE], EtaBinBorders[iE+1], Rvals[Rvalue]*0.6), "goff");
            hjetRatioE_Eta[iE][ipt] = (TH1D *)((TH1D *)jetFile->Get(Form("hjetRatioE_Eta_%d_%d", iE, ipt)))->Clone();
            // hjetRatioE_Eta[iE][ipt]->Fit("f1");

            EtamediansE[iE][ipt] = GetMedian(hjetRatioE_Eta[iE][ipt]);
            EtameansE[iE][ipt] = hjetRatioE_Eta[iE][ipt]->GetMean();
            EtaSDE[iE][ipt] = hjetRatioE_Eta[iE][ipt]->GetStdDev();
            // EtamediansE[iE][ipt]=f1->GetParameter(1);
            // EtameansE[iE][ipt]=f1->GetParameter(1);
            // EtaSDE[iE][ipt]=f1->GetParameter(2);

            hEtaMedianE[iE]->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2.0), EtamediansE[iE][ipt]);
            hEtaMedianE[iE]->SetBinError(ipt + 1, hjetRatioE_Eta[iE][ipt]->GetMeanError());

            hEtaMeanE[iE]->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2.0), EtameansE[iE][ipt]);
            hEtaMeanE[iE]->SetBinError(ipt + 1, hjetRatioE_Eta[iE][ipt]->GetMeanError());
            hEtaSDE[iE]->Fill(JetEBorders[ipt + 1] - ((JetEBorders[ipt + 1] - JetEBorders[ipt]) / 2.0), EtaSDE[iE][ipt]);
            hEtaSDE[iE]->SetBinError(ipt + 1, hjetRatioE_Eta[iE][ipt]->GetStdDevError());

            // hjetRatioE_Eta[iE][ipt]->Scale(1., "width");
            hjetRatioE_Eta[iE][ipt]->Scale(1. / hjetRatioE_Eta[iE][ipt]->Integral(), "");
            // hjetRatioE_Eta[iE][ipt]->Scale(1./hjetRatioE_Eta[iE][ipt]->GetEntries());
            // cout << "In jet pT range " << JetEBorders[ipt] << " - " << JetEBorders[ipt+1] << "GeV/c, number of jets accepted = " << hjetRatioE_Eta[iE][ipt]->GetEntries() << endl;
        }
    }

    fout->Write();
    fout->Close();
}

double GetMedian(TH1D *h)
{

    double x, q;
    q = 0.5;              // 0.5 for "median"
    h->ComputeIntegral(); // just a precaution
    h->GetQuantiles(1, &x, &q);
    std::cout << "median = " << x << std::endl;
    return x;
}
