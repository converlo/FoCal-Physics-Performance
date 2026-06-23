#include <TROOT.h>
#include <TStyle.h>
#include <TColor.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TSpline.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "paperPlotsHeader.h"

double textSize = 0.05; // 5;

// Macro to plot the pT binned JES and JER as well as response matrices and projections
// Base macro provided to me by F. Bock

void focalJetResolutions(
    TString inputfileR02 = "Data20230816/JES/EneMerged_OutputR6_pT_0Mass.root",
    TString inputfileR04 = "Data20230816/JES/EneMerged_OutputR6_pT_0Mass.root")
{
  StyleSettingsPaper();
  TGaxis::SetMaxDigits(4);

  int RMax = 1;      // how many R values you want to draw
  int JESMinX = 5;   // xmin for the jes and jer plots, also used for other pt min
  int JESMaxX = 102; // xmin for the jes and jer plots, also used for other pT max

  Color_t colors[3] = {kBlue + 2, kRed + 1, kBlack};   // kBlack,  kRed+1, kBlue+2
  Color_t colorsMC[3] = {kBlack, kRed + 1, kBlue + 2}; //{kBlue+2,  kRed+1, kBlack}; //{kGray+1, kRed+3, kBlue+3};
  Style_t marker[3] = {20, 21, 33};
  Style_t markerMC[3] = {33, 20, 21}; // 24, 25, 27
  Style_t style[3] = {1, 5, 7};
  Size_t markerS[3] = {2, 2, 2}; // 2.4

  TString radiusOut[2] = {"R06", "R04"};
  TString radiusLabel[2] = {"#it{R} = 0.6", "#it{R} = 0.4"};
  TString etaOut[3] = {"Full", "40to45", "45to49"};
  TString etaRange[3] = {"4.0 < #it{#eta}_{jet} < 4.9", "4.0 < #it{#eta}_{jet} < 4.5", "4.5 < #it{#eta}_{jet} < 4.9"}; //{3.8, 4.5, 5.1}  = {"3.4+R < #it{#eta}_{jet} < 5.5-R", "4.0 < #it{#eta}_{jet} < 4.5", "4.5 < #it{#eta}_{jet} < 4.9"};
  double rangeJES[2][2] = {{-0.35, -0.05}, {-0.35, -0.05}};
  const int maxNPtbins = 13;
  int exampleBins[3] = {0, 2, 5}; // 0,2,5 5,6,7  3,5,6
  // double binningPt[10]  = {2., 5., 10., 15., 20., 25., 30., 35., 40., 70.};
  // double binningPt[10]  = {2.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 400.0}; //for the larger pt bin
  // double binningPt[14]= {2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0};// 400.0};
  double binningPt[14] = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 150.0}; // nPtBins = 14

  TH2D *histResponseMat_pT[2][3];
  TH1D *histMean_pT[2][3];
  TH1D *histMedian_pT[2][3];
  TH1D *histSigma_pT[2][3];
  TH1D *histDeltaPt_bins[2][3][maxNPtbins];

  TH1D *hgausMeanspT;
  TH1D *hgausSDpT;

  TFile *fileR02 = new TFile(inputfileR02.Data());
  for (int i = 0; i < 3; i++)
  {
    if (i == 0)
    {
      histResponseMat_pT[0][i] = (TH2D *)fileR02->Get("hRespMatrix_pT");
      histMean_pT[0][i] = (TH1D *)fileR02->Get("hMeanpT");
      histMedian_pT[0][i] = (TH1D *)fileR02->Get("hMedianpT");
      histSigma_pT[0][i] = (TH1D *)fileR02->Get("hSDpT");
      hgausMeanspT = (TH1D *)fileR02->Get("hgausMeanspT");
      hgausSDpT = (TH1D *)fileR02->Get("hgausSDpT");
      for (int p = 0; p < maxNPtbins; p++)
      {
        histDeltaPt_bins[0][i][p] = (TH1D *)fileR02->Get(Form("hjetRatiopT_%d", p));
        histDeltaPt_bins[0][i][p]->Scale(histDeltaPt_bins[0][i][p]->GetBinWidth(1));
        histDeltaPt_bins[0][i][p]->Scale(1. / histDeltaPt_bins[0][i][p]->Integral());
      }
    }
    else
    {
      histResponseMat_pT[0][i] = (TH2D *)fileR02->Get(Form("hRespMatrix_pT_Eta_%i", i - 1));
      histMean_pT[0][i] = (TH1D *)fileR02->Get(Form("hEtaMeanpT_%i", i - 1));
      histMedian_pT[0][i] = (TH1D *)fileR02->Get(Form("hEtaMedianpT_%i", i - 1));
      histSigma_pT[0][i] = (TH1D *)fileR02->Get(Form("hEtaSDpT_%i", i - 1));
      for (int p = 0; p < maxNPtbins; p++)
      {
        histDeltaPt_bins[0][i][p] = (TH1D *)fileR02->Get(Form("hjetRatiopT_Eta_%i_%i", i - 1, p));
        histDeltaPt_bins[0][i][p]->Scale(histDeltaPt_bins[0][i][p]->GetBinWidth(1));
        histDeltaPt_bins[0][i][p]->Scale(1. / histDeltaPt_bins[0][i][p]->Integral());
      }
    }
    std::cout << "R = 0.2 \t" << histMean_pT[0][i] << "\t" << histSigma_pT[0][i] << "\t" << histResponseMat_pT[0][i] << std::endl;
    // histResponseMat_pT[0][i]->Rebin2D(1,2);
    for (int enty = 1; enty <= histResponseMat_pT[0][i]->GetNbinsY(); enty++)
    {
      double sumX = 0.0;
      for (int entx = 1; entx <= histResponseMat_pT[0][i]->GetNbinsX(); entx++)
      {
        // sum all x
        sumX += histResponseMat_pT[0][i]->GetBinContent(histResponseMat_pT[0][i]->GetBin(entx, enty));
      }
      for (int entx = 1; entx <= histResponseMat_pT[0][i]->GetNbinsX(); entx++)
      {
        // normalize with sum
        double oldBin = histResponseMat_pT[0][i]->GetBinContent(histResponseMat_pT[0][i]->GetBin(entx, enty));
        double newBin = 0.0;
        if (sumX != 0.0)
        {
          double newBin = oldBin / sumX;
          histResponseMat_pT[0][i]->SetBinContent(entx, enty, newBin);
        }
      }
    }
  }

  TFile *fileR04 = new TFile(inputfileR04.Data());
  for (int i = 0; i < 3; i++)
  {
    if (i == 0)
    {
      histResponseMat_pT[1][i] = (TH2D *)fileR04->Get("hRespMatrix_pT");
      histMean_pT[1][i] = (TH1D *)fileR04->Get("hMeanpT");
      histMedian_pT[1][i] = (TH1D *)fileR04->Get("hMedianpT");
      histSigma_pT[1][i] = (TH1D *)fileR04->Get("hSDpT");
      for (int p = 0; p < maxNPtbins; p++)
      {
        histDeltaPt_bins[1][i][p] = (TH1D *)fileR04->Get(Form("hjetRatiopT_%d", p));
        histDeltaPt_bins[1][i][p]->Scale(histDeltaPt_bins[1][i][p]->GetBinWidth(1));
        histDeltaPt_bins[1][i][p]->Scale(1. / histDeltaPt_bins[1][i][p]->Integral());
      }
    }
    else
    {
      histResponseMat_pT[1][i] = (TH2D *)fileR04->Get(Form("hRespMatrix_pT_Eta_%i", i - 1));
      histMean_pT[1][i] = (TH1D *)fileR04->Get(Form("hEtaMeanpT_%i", i - 1));
      histMedian_pT[1][i] = (TH1D *)fileR04->Get(Form("hEtaMedianpT_%i", i - 1));
      histSigma_pT[1][i] = (TH1D *)fileR04->Get(Form("hEtaSDpT_%i", i - 1));
      for (int p = 0; p < maxNPtbins; p++)
      {
        histDeltaPt_bins[1][i][p] = (TH1D *)fileR04->Get(Form("hjetRatiopT_Eta_%i_%i", i - 1, p));
        histDeltaPt_bins[1][i][p]->Scale(histDeltaPt_bins[1][i][p]->GetBinWidth(1));
        histDeltaPt_bins[1][i][p]->Scale(1. / histDeltaPt_bins[1][i][p]->Integral());
      }
    }
    std::cout << "R = 0.4 \t" << histMean_pT[1][i] << "\t" << histSigma_pT[1][i] << "\t" << histResponseMat_pT[1][i] << std::endl;
  }

  TCanvas *cJES = new TCanvas("cJES", "cJES", 2 * 600, 2 * 400);
  DrawPaperCanvasSettings(cJES, 0.081, 0.01, 0.013, 0.11);
  cJES->cd();

  TH2D *jesframe = new TH2D("jesframe", "", 5000, JESMinX, JESMaxX, 2000, -0.7, 0.2); // 72 -0.7
  SetStyleHistoTH2ForGraphs(jesframe, "#it{p}_{T,part} (GeV/#it{c})", "JES", 0.85 * textSize, textSize, 0.85 * textSize, textSize, 0.92, 0.85, 510, 505, 42, 62);
  jesframe->GetXaxis()->SetLabelOffset(0.004);
  jesframe->GetYaxis()->SetLabelOffset(0.006);

  TCanvas *cJER = new TCanvas("cJER", "cJER", 2 * 600, 2 * 400);
  DrawPaperCanvasSettings(cJER, 0.081, 0.01, 0.013, 0.11);
  cJER->cd();

  TH2D *jerframe = new TH2D("jerframe", "", 5000, JESMinX, JESMaxX, 2000, 0.045, 0.22); // 72 0.05
  SetStyleHistoTH2ForGraphs(jerframe, "#it{p}_{T,part} (GeV/#it{c})", "JER", 0.85 * textSize, textSize, 0.85 * textSize, textSize, 0.92, 0.85, 510, 505, 42, 62);
  jerframe->GetXaxis()->SetLabelOffset(0.004);
  jerframe->GetYaxis()->SetLabelOffset(0.006);

  for (int r = 0; r < RMax; r++)
  {
    cJES->Clear();
    jesframe->GetYaxis()->SetRangeUser(rangeJES[r][0], rangeJES[r][1]);
    jesframe->Draw("axis");

    TLegend *leg2 = GetAndSetLegend2(0.2, 0.93 - 3 * 1.0 * textSize, 0.3, 0.93, textSize, 1, "", 42, 0.35);
    for (int e = 0; e < 3; e++)
    {
      DrawSetMarker(histMean_pT[r][e], marker[e], markerS[e], colors[e], colors[e]);
      histMean_pT[r][e]->Draw("same");
      leg2->AddEntry(histMean_pT[r][e], etaRange[e].Data(), "p");
    }
    leg2->Draw("same");
    drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV", 0.95, 0.965 - 1 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd("FoCal upgrade", 0.95, 0.965 - 2 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s", radiusLabel[r].Data()), 0.95, 0.965 - 3 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    jesframe->Draw("axis,same");

    cJES->cd();
    cJES->Update();
    cJES->SaveAs(Form("figs/PT/JES_%s.pdf", radiusOut[r].Data()));

    cJER->Clear();
    jerframe->Draw("axis");
    TLegend *leg3 = GetAndSetLegend2(0.74, 0.94 - 6 * 1.1 * textSize, 0.96, 0.94 - 3 * 1.1 * textSize, textSize, 1, "", 42, 0.25);
    for (int e = 0; e < 3; e++)
    {
      DrawSetMarker(histSigma_pT[r][e], marker[e], markerS[e], colors[e], colors[e]);
      histSigma_pT[r][e]->Draw("same");
      leg3->AddEntry(histSigma_pT[r][e], etaRange[e].Data(), "p");
    }
    leg3->Draw("same");
    drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV", 0.95, 0.965 - 1 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd("FoCal upgrade", 0.95, 0.965 - 2 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s", radiusLabel[r].Data()), 0.95, 0.965 - 3 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    jerframe->Draw("axis,same");

    cJER->cd();
    cJER->Update();
    cJER->SaveAs(Form("figs/PT/JER_%s.pdf", radiusOut[r].Data()));
  }

  for (int e = 0; e < 3; e++)
  {
    cJES->Clear();
    jesframe->GetYaxis()->SetRangeUser(rangeJES[0][0], rangeJES[1][1]);
    jesframe->Draw("axis");

    TLegend *leg2 = GetAndSetLegend2(0.2, 0.8 - 3 * 1.0 * textSize, 0.3, 0.82, textSize, 1, "", 42, 0.35);
    for (int r = 0; r < RMax; r++)
    {
      if (r == 1)
        DrawSetMarker(histMean_pT[r][e], marker[e], markerS[e], colors[e], colors[e]);
      else
      {
        DrawSetMarker(histMean_pT[r][e], markerMC[e], markerS[e] * 1.2, colorsMC[e], colorsMC[e]);
        DrawSetMarker(histMedian_pT[r][e], marker[e], markerS[e] * 1.2, colors[e], colors[e]);
        DrawSetMarker(hgausMeanspT, marker[e + 1], markerS[e] * 1.2, colors[e + 1], colors[e + 1]);
      }
      histMean_pT[r][e]->Draw("same");
      histMedian_pT[r][e]->Draw("same");
      hgausMeanspT->Draw("same");
      leg2->AddEntry(histMean_pT[r][e], "mean", "p"); // leg2->AddEntry(histMean_pT[r][e], radiusLabel[r].Data(), "p");
      leg2->AddEntry(histMedian_pT[r][e], "median", "p");
      leg2->AddEntry(hgausMeanspT, "mean from Gaussian fit", "p");
    }
    leg2->Draw("same");
    drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV", 0.95, 0.965 - 1 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd("FoCal upgrade", 0.95, 0.965 - 2 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd(Form("jets, anti-#it{k}_{T}, #it{R} = 0.6, %s", etaRange[e].Data()), 0.95, 0.965 - 3 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    // drawLatexAdd("#it{R}=0.6",0.95,0.965-4.2*1.1*textSize, textSize,kFALSE, kFALSE, kTRUE);
    drawLatexAdd("#Delta#it{p}_{T} = (#it{p}_{T,det} #font[122]{-} #it{p}_{T,part})/#it{p}_{T,part} ", 0.95, 0.965 - 4.5 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    jesframe->Draw("axis,same");

    cJES->cd();
    cJES->Update();
    // TLine *line = new TLine(gPad->GetUxmin(), 1, gPad->GetUxmax(), 1);
    // line->Draw("same");
    auto fline = new TF1("fline", "0.0", 5, 102);
    fline->SetLineColor(kBlack);
    fline->SetLineStyle(2);
    fline->Draw("same");
    cJES->Update();
    cJES->SaveAs(Form("figs/PT/gausJES_%s.pdf", etaOut[e].Data()));

    cJER->Clear();
    jerframe->Draw("axis");
    TLegend *leg3 = GetAndSetLegend2(0.25, 0.78 - 4 * 1.1 * textSize, 0.35, 0.78 - 1 * 1.1 * textSize, textSize, 1, "", 42, 0.3);
    for (int r = 0; r < RMax; r++)
    {
      if (r == 1)
        DrawSetMarker(histSigma_pT[r][e], marker[e], markerS[e], colors[e], colors[e]);
      else
      {
        DrawSetMarker(histSigma_pT[r][e], markerMC[e], markerS[e] * 1.2, colorsMC[e], colorsMC[e]);
        DrawSetMarker(hgausSDpT, marker[e + 1], markerS[e] * 1.2, colors[e + 1], colors[e + 1]);
      }
      histSigma_pT[r][e]->Draw("same");
      hgausSDpT->Draw("same");
      // leg3->AddEntry((TObject*)0, radiusLabel[r].Data(), "");
      leg3->AddEntry(histSigma_pT[r][e], "standard deviation", "p");
      leg3->AddEntry(hgausSDpT, "standard deviation from Gaussian fit", "p");
    }
    leg3->Draw("same");
    drawLatexAdd("ALICE simulation, pp #sqrt{#it{s}} = 14 TeV", 0.95, 0.965 - 1 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd("FoCal upgrade", 0.95, 0.965 - 2 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd(Form("jets, anti-#it{k}_{T}, #it{R} = 0.6, %s", etaRange[e].Data()), 0.95, 0.965 - 3 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    drawLatexAdd("#Delta#it{p}_{T} = (#it{p}_{T,det} #font[122]{-} #it{p}_{T,part})/#it{p}_{T,part} ", 0.95, 0.965 - 4.5 * 1.1 * textSize, textSize, kFALSE, kFALSE, kTRUE);
    jerframe->Draw("axis,same");

    cJER->cd();
    cJER->Update();
    cJER->SaveAs(Form("figs/PT/gausJER_%s.pdf", etaOut[e].Data()));
  }

  // Style projection
  TCanvas *cCorr = new TCanvas("cCorr", "cCorr", 2 * 900, 2 * 860);
  // DrawPaperCanvasSettings(cCorr, 0.095, 0.11, 0.02, 0.1);1386
  DrawPaperCanvasSettings(cCorr, 0.12, 0.121, 0.092, 0.11);
  cCorr->SetCanvasSize(2 * 900, 2 * 860);
  cCorr->SetLogz(1);

  for (int e = 0; e < 3; e++)
  {
    for (int r = 0; r < RMax; r++)
    {
      cCorr->Clear();
      SetStyleHistoTH2ForGraphs(histResponseMat_pT[r][e], "#it{p}_{T, det} (GeV/#it{c})", "#it{p}_{T, part} (GeV/#it{c})", 0.75 * textSize * 1.1, textSize * 1.1, 0.75 * textSize * 1.1, textSize * 1.1, 0.87, 0.95, 510, 510, 42, 62);
      histResponseMat_pT[r][e]->GetZaxis()->SetLabelSize(0.85 * textSize);
      histResponseMat_pT[r][e]->GetZaxis()->SetLabelOffset(- 0.005);
      // histResponseMat_pT[r][e]->GetZaxis()->SetRangeUser(0.00001, 220);
      histResponseMat_pT[r][e]->GetXaxis()->SetRangeUser(JESMinX, JESMaxX);
      histResponseMat_pT[r][e]->GetYaxis()->SetRangeUser(JESMinX, JESMaxX);
      histResponseMat_pT[r][e]->Draw("colz");
      auto fline = new TF1("fline", "x", JESMinX, JESMaxX);
      fline->SetLineColor(kGray);
      fline->SetLineStyle(2);
      fline->Draw("same");

      drawLatexAdd("#bf{ALICE simulation, pp #sqrt{#it{s}} = 14 TeV}", 0.15, 1.01 - 1 * 0.9 * textSize, 0.85 * textSize, kFALSE, kFALSE, kFALSE);
      drawLatexAdd("#bf{FoCal upgrade}", 0.15, 1.01 - 2 * 0.9 * textSize, 0.85 * textSize, kFALSE, kFALSE, kFALSE);
      drawLatexAdd(Form("#bf{jets, anti-#it{k}_{T}, %s}", radiusLabel[r].Data()), 0.84, 0.16 + 0.9 * textSize, 0.85 * textSize, kFALSE, kFALSE, kTRUE);
      drawLatexAdd(Form("#bf{%s}", etaRange[e].Data()), 0.84, 0.16, 0.85 * textSize, kFALSE, kFALSE, kTRUE);

      cCorr->Update();
      cCorr->SaveAs(Form("figs/PT/ResponseMatrix_%s_%s.pdf", etaOut[e].Data(), radiusOut[r].Data()));
    }
  }

  double maxY = 0.25; // 0.151;
  int yDivs = 505;
  double pixelX = 1400 * 2;
  double pixelY = 500 * 2;
  double arrayBoundariesX1_4[4];
  double arrayBoundariesY1_4[2];
  double relativeMarginsX[3];
  double relativeMarginsY[3];
  int textSizeLabelsPixel = (int)(textSize * pixelY);
  ReturnCorrectValuesForCanvasScaling(pixelX, pixelY, 3, 1, 0.048, 0.005, 0.01, 0.1, arrayBoundariesX1_4, arrayBoundariesY1_4, relativeMarginsX, relativeMarginsY);
  double textSizeLabels[3] = {0.045, 0.045, 0.045};
  double textSizeFac[3] = {0.045, 0.045, 0.045};
  double margin = relativeMarginsX[0] * pixelX;

  TCanvas *cSlicesJES = new TCanvas("JESslice", "JESslice", pixelX, pixelY);
  DrawPaperCanvasSettings(cSlicesJES, 0, 0, 0, 0);

  std::cout << arrayBoundariesX1_4[0] << "\t" << arrayBoundariesY1_4[0] << "\t" << arrayBoundariesX1_4[1] << "\t" << arrayBoundariesY1_4[1] << std::endl;
  TPad *padRight = new TPad("padRight", "padRight", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0]);
  DrawPaperPadSettings(padRight, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[2]);
  ReturnCorrectValuesTextSize(padRight, textSizeLabels[0], textSizeFac[0], textSizeLabelsPixel, 0.08);

  std::cout << arrayBoundariesX1_4[1] << "\t" << arrayBoundariesY1_4[0] << "\t" << arrayBoundariesX1_4[2] << "\t" << arrayBoundariesY1_4[1] << std::endl;
  TPad *padMiddle = new TPad("padMiddle", "padMiddle", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0]);
  DrawPaperPadSettings(padMiddle, relativeMarginsX[1], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[2]);
  ReturnCorrectValuesTextSize(padMiddle, textSizeLabels[1], textSizeFac[1], textSizeLabelsPixel, 0.08);

  std::cout << arrayBoundariesX1_4[2] << "\t" << arrayBoundariesY1_4[0] << "\t" << arrayBoundariesX1_4[3] << "\t" << arrayBoundariesY1_4[1] << std::endl;
  TPad *padLeft = new TPad("padLeft", "padLeft", arrayBoundariesX1_4[2], arrayBoundariesY1_4[1], arrayBoundariesX1_4[3], arrayBoundariesY1_4[0]);
  DrawPaperPadSettings(padLeft, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[2]);
  ReturnCorrectValuesTextSize(padLeft, textSizeLabels[2], textSizeFac[2], textSizeLabelsPixel, 0.08);

  for (int r = 0; r < RMax; r++)
  {
    padRight->Draw();
    padMiddle->Draw();
    padLeft->Draw();

    padRight->cd();
    TH2D *esframeR = new TH2D("effframe", "", 5000, -1.01, 0.5, 2000, 0., maxY);
    SetStyleHistoTH2ForGraphs(esframeR, "#Delta#it{p}_{T} = (#it{p}_{T,det} #font[122]{-} #it{p}_{T,part})/#it{p}_{T,part} ", "1/#it{k} d#it{N}/d#Delta#it{p}_{T}", 0.85 * textSizeLabels[0], textSizeLabels[0], 0.85 * textSizeLabels[0], textSizeLabels[0], 0.955, 1.45, 510, yDivs, 42, 62);
    esframeR->GetYaxis()->SetLabelOffset(0.01);
    esframeR->GetYaxis()->SetTickLength(0.035);
    esframeR->Draw("axis");

    for (int e = 0; e < 3; e++)
    {
      DrawSetMarker(histDeltaPt_bins[r][e][exampleBins[0]], marker[e], markerS[e] * 1.5, colors[e], colors[e]);
      histDeltaPt_bins[r][e][exampleBins[0]]->Draw("same");
    }

    drawLatexAdd("ALICE simulation", 0.17, 0.965 - 1 * 1.1 * textSizeLabels[0], textSizeLabels[0], kFALSE, kFALSE, kFALSE);
    drawLatexAdd("FoCal upgrade", 0.17, 0.965 - 2 * 1.1 * textSizeLabels[0], textSizeLabels[0], kFALSE, kFALSE, kFALSE);
    drawLatexAdd("pp #sqrt{#it{s}} = 14 TeV", 0.17, 0.965 - 3 * 1.1 * textSizeLabels[0], textSizeLabels[0], kFALSE, kFALSE, kFALSE);
    drawLatexAdd(Form("%.0f< #it{p}_{T,part}<%.0f GeV/#it{c}", binningPt[exampleBins[0]], binningPt[exampleBins[0] + 1]), 0.96, 0.97 - 1 * 1.1 * textSizeLabels[0], textSizeLabels[0], kFALSE, kFALSE, kTRUE);

    DrawLines(0, 0, 0, 0.05, 2, kGray + 2, 7);

    esframeR->Draw("axis,same");
    padMiddle->cd();
    TH2D *esframeM = new TH2D("effframe", "", 5000, -1.01, 0.5, 2000, 0., maxY);
    SetStyleHistoTH2ForGraphs(esframeM, "#Delta#it{p}_{T} = (#it{p}_{T,det} #font[122]{-} #it{p}_{T,part})/#it{p}_{T,part} ", "1/k d#it{N}/d#Delta#it{p}_{T}", 0.85 * textSizeLabels[1], textSizeLabels[1], 0.85 * textSizeLabels[1], textSizeLabels[1], 0.84, 0.95, 510, yDivs, 42, 62);
    esframeM->GetYaxis()->SetTickLength(0.04);
    esframeM->GetXaxis()->SetLabelOffset(-0.001);
    esframeM->Draw("axis");

    TLegend *leg2 = GetAndSetLegend2(0.03, 0.945 - 3 * 1. * textSizeLabels[1], 0.4, 0.945, textSizeLabels[1], 1, "", 42, 0.35);
    for (int e = 0; e < 3; e++)
    {
      DrawSetMarker(histDeltaPt_bins[r][e][exampleBins[1]], marker[e], markerS[e] * 1.5, colors[e], colors[e]);
      histDeltaPt_bins[r][e][exampleBins[1]]->Draw("same");
      leg2->AddEntry(histDeltaPt_bins[r][e][exampleBins[1]], etaRange[e].Data(), "p");
    }
    leg2->Draw("same");
    drawLatexAdd(Form("%.0f< #it{p}_{T,part}<%.0f GeV/#it{c}", binningPt[exampleBins[1]], binningPt[exampleBins[1] + 1]), 0.95, 0.97 - 1 * 1.1 * textSizeLabels[1], textSizeLabels[1], kFALSE, kFALSE, kTRUE);

    DrawLines(0, 0, 0, 0.05, 2, kGray + 2, 7);
    esframeM->Draw("axis,same");

    padLeft->cd();
    TH2D *esframeL = new TH2D("esframeL", "", 5000, -1.01, 0.5, 2000, 0., maxY);
    SetStyleHistoTH2ForGraphs(esframeL, "#Delta#it{p}_{T} = (#it{p}_{T,det} #font[122]{-} #it{p}_{T,part})/#it{p}_{T,part} ", "#it{p}", 0.85 * textSizeLabels[1], textSizeLabels[1], 0.85 * textSizeLabels[1], textSizeLabels[1], 0.83, 0.95, 510, yDivs, 42, 62);
    esframeL->GetYaxis()->SetTickLength(0.038);
    esframeL->GetXaxis()->SetLabelOffset(-0.001);
    esframeL->Draw("axis");

    for (int e = 0; e < 3; e++)
    {
      DrawSetMarker(histDeltaPt_bins[r][e][exampleBins[2]], marker[e], markerS[e] * 1.5, colors[e], colors[e]);
      histDeltaPt_bins[r][e][exampleBins[2]]->Draw("same");
    }
    drawLatexAdd(Form("jets, anti-#it{k}_{T}, %s", radiusLabel[r].Data()), 0.03, 0.965 - 1 * 1.1 * textSizeLabels[2], textSizeLabels[2], kFALSE, kFALSE, kFALSE);
    drawLatexAdd(Form("%.0f< #it{p}_{T,part}<%.0f GeV/#it{c}", binningPt[exampleBins[2]], binningPt[exampleBins[2] + 1]), 0.93, 0.965 - 1 * 1.1 * textSizeLabels[2], textSizeLabels[2], kFALSE, kFALSE, kTRUE);

    DrawLines(0, 0, 0, 0.05, 2, kGray + 2, 7);
    esframeL->Draw("axis,same");

    cSlicesJES->cd();
    cSlicesJES->Update();
    cSlicesJES->SaveAs(Form("figs/PT/JetEscaleProj_%s.pdf", radiusOut[r].Data()));
  }

  for (int e = 0; e < 3; e++)
  {
    padRight->Draw();
    padMiddle->Draw();
    padLeft->Draw();

    padRight->cd();
    TH2D *esframeR = new TH2D("effframe", "", 5000, -1.01, 0.5, 2000, 0., maxY);
    SetStyleHistoTH2ForGraphs(esframeR, "#Delta#it{p}_{T} = (#it{p}_{T,det} #font[122]{-} #it{p}_{T,part})/#it{p}_{T,part} ", "1/#it{k} d#it{N}/d#Delta#it{p}_{T}", 0.85 * textSizeLabels[0], textSizeLabels[0], 0.85 * textSizeLabels[0], textSizeLabels[0], 0.955, 1.45, 510, yDivs, 42, 62);
    esframeR->GetYaxis()->SetLabelOffset(0.01);
    esframeR->GetYaxis()->SetTickLength(0.035);
    esframeR->Draw("axis");

    for (int r = 0; r < RMax; r++)
    {
      if (r == 1)
        DrawSetMarker(histDeltaPt_bins[r][e][exampleBins[0]], marker[e], markerS[e] * 1.5, colors[e], colors[e]);
      else
        DrawSetMarker(histDeltaPt_bins[r][e][exampleBins[0]], markerMC[e], markerS[e] * 1.9, colorsMC[e], colorsMC[e]);
      histDeltaPt_bins[r][e][exampleBins[0]]->Draw("same");
    }

    drawLatexAdd("ALICE simulation", 0.17, 0.965 - 1 * 1.1 * textSizeLabels[0], textSizeLabels[0], kFALSE, kFALSE, kFALSE);
    drawLatexAdd("FoCal upgrade", 0.17, 0.965 - 2 * 1.1 * textSizeLabels[0], textSizeLabels[0], kFALSE, kFALSE, kFALSE);
    drawLatexAdd("pp #sqrt{#it{s}} = 14 TeV", 0.17, 0.965 - 3 * 1.1 * textSizeLabels[0], textSizeLabels[0], kFALSE, kFALSE, kFALSE);
    drawLatexAdd(Form("%.0f< #it{p}_{T,part}<%.0f GeV/#it{c}", binningPt[exampleBins[0]], binningPt[exampleBins[0] + 1]), 0.96, 0.97 - 1 * 1.1 * textSizeLabels[0], textSizeLabels[0], kFALSE, kFALSE, kTRUE);

    DrawLines(0, 0, 0, 0.05, 2, kGray + 2, 7);

    esframeR->Draw("axis,same");
    padMiddle->cd();
    TH2D *esframeM = new TH2D("effframe", "", 5000, -1.01, 0.5, 2000, 0., maxY);
    SetStyleHistoTH2ForGraphs(esframeM, "#Delta#it{p}_{T} = (#it{p}_{T,det} #font[122]{-} #it{p}_{T,part})/#it{p}_{T,part} ", "1/k d#it{N}/d#Delta#it{p}_{T}", 0.85 * textSizeLabels[1], textSizeLabels[1], 0.85 * textSizeLabels[1], textSizeLabels[1], 0.84, 0.95, 510, yDivs, 42, 62);
    esframeM->GetYaxis()->SetTickLength(0.04);
    esframeM->GetXaxis()->SetLabelOffset(-0.001);
    esframeM->Draw("axis");

    TLegend *leg2 = GetAndSetLegend2(0.03, 0.945 - 2 * 1. * textSizeLabels[1], 0.4, 0.945, textSizeLabels[1], 1, "", 42, 0.35);
    for (int r = 0; r < RMax; r++)
    {
      if (r == 1)
        DrawSetMarker(histDeltaPt_bins[r][e][exampleBins[1]], marker[e], markerS[e] * 1.5, colors[e], colors[e]);
      else
        DrawSetMarker(histDeltaPt_bins[r][e][exampleBins[1]], markerMC[e], markerS[e] * 1.9, colorsMC[e], colorsMC[e]);
      histDeltaPt_bins[r][e][exampleBins[1]]->Draw("same");
      leg2->AddEntry(histDeltaPt_bins[r][e][exampleBins[1]], radiusLabel[r].Data(), "p");
    }
    leg2->Draw("same");
    drawLatexAdd(Form("%.0f< #it{p}_{T,part}<%.0f GeV/#it{c}", binningPt[exampleBins[1]], binningPt[exampleBins[1] + 1]), 0.95, 0.97 - 1 * 1.1 * textSizeLabels[1], textSizeLabels[1], kFALSE, kFALSE, kTRUE);

    DrawLines(0, 0, 0, 0.05, 2, kGray + 2, 7);
    esframeM->Draw("axis,same");

    padLeft->cd();
    TH2D *esframeL = new TH2D("esframeL", "", 5000, -1.01, 0.5, 2000, 0., maxY);
    SetStyleHistoTH2ForGraphs(esframeL, "#Delta#it{p}_{T} = (#it{p}_{T,det} #font[122]{-} #it{p}_{T,part})/#it{p}_{T,part} ", "#it{p}", 0.85 * textSizeLabels[1], textSizeLabels[1], 0.85 * textSizeLabels[1], textSizeLabels[1], 0.83, 0.95, 510, yDivs, 42, 62);
    esframeL->GetYaxis()->SetTickLength(0.038);
    esframeL->GetXaxis()->SetLabelOffset(-0.001);
    esframeL->Draw("axis");

    for (int r = 0; r < RMax; r++)
    {
      if (r == 1)
        DrawSetMarker(histDeltaPt_bins[r][e][exampleBins[2]], marker[e], markerS[e] * 1.5, colors[e], colors[e]);
      else
        DrawSetMarker(histDeltaPt_bins[r][e][exampleBins[2]], markerMC[e], markerS[e] * 1.9, colorsMC[e], colorsMC[e]);

      histDeltaPt_bins[r][e][exampleBins[2]]->Draw("same");
    }
    drawLatexAdd("jets, anti-#it{k}_{T}", 0.03, 0.965 - 1 * 1.1 * textSizeLabels[2], textSizeLabels[2], kFALSE, kFALSE, kFALSE);
    drawLatexAdd(Form("%s", etaRange[e].Data()), 0.03, 0.965 - 2 * 1.1 * textSizeLabels[2], textSizeLabels[2], kFALSE, kFALSE, kFALSE);
    drawLatexAdd(Form("%.0f< #it{p}_{T,part}<%.0f GeV/#it{c}", binningPt[exampleBins[2]], binningPt[exampleBins[2] + 1]), 0.93, 0.965 - 1 * 1.1 * textSizeLabels[2], textSizeLabels[2], kFALSE, kFALSE, kTRUE);

    DrawLines(0, 0, 0, 0.05, 2, kGray + 2, 7);
    esframeL->Draw("axis,same");

    cSlicesJES->cd();
    cSlicesJES->Update();
    cSlicesJES->SaveAs(Form("figs/PT/JetEscaleProj_%s.pdf", etaOut[e].Data()));
  }
}
