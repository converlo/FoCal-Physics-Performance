//__________________________________________________________________________________________________________

void SetPlotStyle() {
  const Int_t nRGBs = 5;
  const Int_t nCont = 255;

  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.31, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1., 0.12, 0.00, 0.00};

  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
  gStyle->SetNumberContours(nCont);
}

void StyleSettingsPaper( TString format = ""){
    //gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);   //show day and time
    gStyle->SetOptStat(0);  //show statistic
    gStyle->SetPalette(1,0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTextSize(0.5);
    gStyle->SetLabelSize(0.03,"xyz");
    gStyle->SetLabelOffset(0.002,"xyz");
    gStyle->SetTitleFontSize(0.04);
    gStyle->SetTitleOffset(1,"y");
    gStyle->SetTitleOffset(0.7,"x");
    gStyle->SetCanvasColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(1);

    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.09);
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadLeftMargin(0.13);


    TGaxis::SetMaxDigits(3);
    gErrorIgnoreLevel=kError;

    if (format.CompareTo("eps") == 0 ||format.CompareTo("pdf") == 0  ) gStyle->SetLineScalePS(1);
    SetPlotStyle();
}
//
//  TCanvas *canvasExample  = new TCanvas("canvasExample","",0,0,1300,850);
//  Example:     DrawPaperCanvasSettings( canvasExample, 0.085, 0.01, 0.01, 0.105);
//__________________________________________________________________________________________________________
void DrawPaperCanvasSettings(
    TCanvas* c1,
    Double_t leftMargin,
    Double_t rightMargin,
    Double_t topMargin,
    Double_t bottomMargin
){
    c1->SetTickx();
    c1->SetTicky();
    c1->SetGridx(0);
    c1->SetGridy(0);
    c1->SetLogy(0);
    c1->SetLeftMargin(leftMargin);
    c1->SetRightMargin(rightMargin);
    c1->SetTopMargin(topMargin);
    c1->SetBottomMargin(bottomMargin);
    c1->SetFillColor(0);
}

//
//  TPad *padExample  = new TPad("padExample","",0,0,1300,850);
//  Example:     DrawPaperPadSettings( padExample, 0.085, 0.01, 0.01, 0.105);
//__________________________________________________________________________________________________________
void DrawPaperPadSettings(
    TPad* c1,
    Double_t leftMargin,
    Double_t rightMargin,
    Double_t topMargin,
    Double_t bottomMargin
){
    c1->SetTickx();
    c1->SetTicky();
    c1->SetGridx(0);
    c1->SetGridy(0);
    c1->SetLogy(0);
    c1->SetLeftMargin(leftMargin);
    c1->SetRightMargin(rightMargin);
    c1->SetTopMargin(topMargin);
    c1->SetBottomMargin(bottomMargin);
    c1->SetFillColor(0);
}

//
//     Example
//     Double_t minY                                 = 0.91;
//     Double_t maxY                                 = 1.039;
//     Double_t minX                                 = 0.27;
//     Double_t maxX                                 = 2.99e2;
//     Double_t textSizeSinglePad                    = 0.05;
//     Double_t textSizeLabelsPixel                  = 35;
//     Double_t textSizeLabelsRel                    = 35./canvasheight;
//
// it is often best to use a dummy histogram for the style settings
//     TH2F * histExampleDummy    = new TH2F("histExampleDummy","histExampleDummy",1000,minX, maxX,1000,minY, maxY);
// set the style of the dummy (dummy, x title, y title, x label size, x title size, y label size, y title size, x title offset, y title offset, ndivisions x, ndivisions y)
//     SetStyleHistoTH2ForGraphs(histExampleDummy, "#it{E}_{rec} (GeV)","#it{E}_{rec}/#it{E}_{in}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81, 510, 510);
//
//   histExampleDummy->GetXaxis()->SetNoExponent(); // when SetLogx() is used, one can think about adding more labels
//   histExampleDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
//   histExampleDummy->DrawCopy();
//
//__________________________________________________________________________________________________________
void SetStyleHistoTH2ForGraphs(
    TH2* histo,
    TString XTitle,
    TString YTitle,
    Size_t xLableSize,
    Size_t xTitleSize,
    Size_t yLableSize,
    Size_t yTitleSize,
    Float_t xTitleOffset    = 1,
    Float_t yTitleOffset    = 1,
    Int_t xNDivisions       = 510,
    Int_t yNDivisions       = 510,
    Font_t textFontLabel    = 42,
    Font_t textFontTitle    = 62,
    TString ZTitle          = "",
    Size_t zLableSize       = 0.03,
    Size_t zTitleSize       = 0.03,
    Float_t zTitleOffset    = 1,
    Int_t zNDivisions       = 510
    
){
    histo->SetXTitle(XTitle);
    histo->SetYTitle(YTitle);
    histo->SetZTitle(ZTitle);
    histo->SetTitle("");

    histo->GetXaxis()->SetLabelFont(textFontLabel);
    histo->GetYaxis()->SetLabelFont(textFontLabel);
    histo->GetZaxis()->SetLabelFont(textFontLabel);
    histo->GetXaxis()->SetTitleFont(textFontTitle);
    histo->GetYaxis()->SetTitleFont(textFontTitle);
    histo->GetZaxis()->SetTitleFont(textFontTitle);

    histo->GetXaxis()->SetLabelSize(xLableSize);
    histo->GetXaxis()->SetTitleSize(xTitleSize);
    histo->GetXaxis()->SetTitleOffset(xTitleOffset);
    histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    histo->GetYaxis()->SetDecimals();
    histo->GetYaxis()->SetLabelSize(yLableSize);
    histo->GetYaxis()->SetTitleSize(yTitleSize);
    histo->GetYaxis()->SetTitleOffset(yTitleOffset);
    histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);

    histo->GetZaxis()->SetLabelSize(zLableSize);
    histo->GetZaxis()->SetTitleSize(zTitleSize);
    histo->GetZaxis()->SetTitleOffset(zTitleOffset);
    histo->GetZaxis()->SetNdivisions(zNDivisions,kTRUE);

}

//__________________________________________________________________________________________________________
void SetStyleHistoTH1ForGraphs(
    TH1* histo,
    TString XTitle,
    TString YTitle,
    Size_t xLableSize,
    Size_t xTitleSize,
    Size_t yLableSize,
    Size_t yTitleSize,
    Float_t xTitleOffset    = 1,
    Float_t yTitleOffset    = 1,
    Int_t xNDivisions       = 510,
    Int_t yNDivisions       = 510,
    Font_t textFontLabel    = 42,
    Font_t textFontTitle    = 62    
){
    histo->SetXTitle(XTitle);
    histo->SetYTitle(YTitle);
    histo->SetTitle("");

    histo->GetXaxis()->SetLabelFont(textFontLabel);
    histo->GetYaxis()->SetLabelFont(textFontLabel);
    histo->GetXaxis()->SetTitleFont(textFontTitle);
    histo->GetYaxis()->SetTitleFont(textFontTitle);

    histo->GetXaxis()->SetLabelSize(xLableSize);
    histo->GetXaxis()->SetTitleSize(xTitleSize);
    histo->GetXaxis()->SetTitleOffset(xTitleOffset);
    histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    histo->GetYaxis()->SetDecimals();
    histo->GetYaxis()->SetLabelSize(yLableSize);
    histo->GetYaxis()->SetTitleSize(yTitleSize);
    histo->GetYaxis()->SetTitleOffset(yTitleOffset);
    histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);

}


//
//  draw a line to guide the eye (good for some comparisons)
//  min x, max x, min y, max y, line width, line color, line style)
//  DrawLines(minX, maxX, 1., 1., 1, kGray+2, 7);
//
//__________________________________________________________________________________________________________
TLine* DrawLines(
    Float_t startX, Float_t endX,
    Float_t startY, Float_t endY,
    Float_t linew, Float_t lineColor = 4, Style_t lineStyle = 1
){
    TLine * l1 = new TLine (startX,startY,endX,endY);
    l1->SetLineColor(lineColor);
    l1->SetLineWidth(linew);
    l1->SetLineStyle(lineStyle);
    l1->Draw("same");
    return l1;
}

//
//  draw a line to guide the eye (good for some comparisons)
//  min x, max x, min y, max y, line width, line color, line style)
//  DrawLines(minX, maxX, 1., 1., 1, kGray+2, 7);
//
//__________________________________________________________________________________________________________
TArrow* DrawArrow(
    Float_t startX, Float_t endX,
    Float_t startY, Float_t endY,
    Float_t angle = 0, Float_t lengthArrow = 0.05,
    Float_t linew = 1, Float_t lineColor = kBlack, Style_t lineStyle = 1, TString arrowHead = ">",Color_t fillColor = kBlack
){
    
    TArrow *ar4 = new TArrow(startX,startY,endX,endY,lengthArrow,arrowHead);
    ar4->SetLineColor(lineColor);
    ar4->SetLineWidth(linew);
    ar4->SetLineStyle(lineStyle);
    ar4->SetFillColor(fillColor);
    ar4->Draw("same");
    return ar4;
}

// function to set the plotting style (takes graph, markerstyle, markersize, and two times the color (for marker and lines))
// Example
//
//  Color_t  colorData          = kBlack;
//  Style_t  markerStyleData    = 20; // do not use asymmetric markers (like triangles [numbers 23, 26, 32 are forbidden!])
//  Size_t   markerSize         = 1.5; //1.5 or 2 make reasonable sizes on most paper plots
//   DrawSetMarker(  histoData,   30,               markerSize,     kOrange-8,     kOrange-8);
//
//__________________________________________________________________________________________________________
void DrawSetMarker(
    TH1* histo1,
    Style_t markerStyle,
    Size_t markerSize,
    Color_t markerColor,
    Color_t lineColor,
    Size_t lineWidth = 0,
    Style_t lineStyle = 1
) {
    histo1->SetMarkerStyle(markerStyle);
    histo1->SetMarkerSize(markerSize);
    histo1->SetMarkerColor(markerColor);
    histo1->SetLineColor(lineColor);
    if (lineStyle != 1) histo1->SetLineStyle(lineStyle);
    if (lineWidth > 0) histo1->SetLineWidth(lineWidth);
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);
}
//
// function to set the plotting style (takes graph, markerstyle, markersize, and two times the color (for marker and lines))
//  Color_t  colorData          = kBlack;
//  Style_t  markerStyleData    = 20; // do not use asymmetric markers (like triangles [numbers 23, 26, 32 are forbidden!])
//  Size_t   markerSize         = 1.5; //1.5 or 2 make reasonable sizes on most paper plots
//
//  DrawSetMarkerTGraphErr(graphData,        markerStyleData,  markerSize,     colorData,     colorData);
//
//__________________________________________________________________________________________________________
void DrawSetMarkerTGraphErr(
    TGraphErrors* graph,
    Style_t markerStyle,
    Size_t markerSize,
    Color_t markerColor,
    Color_t lineColor,
    Width_t lineWidth       = 1,
    Bool_t boxes            = kFALSE,
    Color_t fillColor       = 0,
    Bool_t isHollow         = kFALSE,
    Style_t lineStyle       = 1
) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
    graph->SetLineStyle(lineStyle);
    if (boxes){
        graph->SetFillColor(fillColor);
        if (fillColor!=0){
            if (!isHollow){
                graph->SetFillStyle(1001);
            } else {
                graph->SetFillStyle(0);
            }
        } else {
            graph->SetFillStyle(0);
        }
    }
}

//__________________________________________________________________________________________________________
void DrawSetMarkerTGraph(
    TGraph* graph,
    Style_t markerStyle,
    Size_t markerSize,
    Color_t markerColor,
    Color_t lineColor,
    Width_t lineWidth       = 1,
    Bool_t boxes            = kFALSE,
    Color_t fillColor       = 0,
    Bool_t isHollow         = kFALSE,
    Style_t lineStyle       = 1
) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
    graph->SetLineStyle(lineStyle);
    if (boxes){
        graph->SetFillColor(fillColor);
        if (fillColor!=0){
            if (!isHollow){
                graph->SetFillStyle(1001);
            } else {
                graph->SetFillStyle(0);
            }
        } else {
            graph->SetFillStyle(0);
        }
    }
}


//__________________________________________________________________________________________________________
void DrawSetMarkerTGraphAsym(
    TGraphAsymmErrors* graph,
    Style_t markerStyle,
    Size_t markerSize,
    Color_t markerColor,
    Color_t lineColor,
    Width_t lineWidth   =1,
    Bool_t boxes        = kFALSE,
    Color_t fillColor   = 0,
    Bool_t isHollow     = kFALSE,
    Style_t lineStyle   = 1
) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
    graph->SetLineStyle(lineStyle);
    if (boxes){
        graph->SetFillColor(fillColor);
        if (fillColor!=0){
            if (!isHollow){
                graph->SetFillStyle(1001);
            } else {
                graph->SetFillStyle(0);
            }
        } else {
            graph->SetFillStyle(0);
        }
    }
}


//
//     Example
//     Double_t textSizeLabelsRel    = 35./canvasheight;
//     drawLatexAdd("ALICE Performance",0.15,0.92,textSizeLabelsRel,kFALSE);
//     drawLatexAdd("Xe-Xe, #sqrt{#it{s}_{NN}} = 90 TeV",0.15,0.92-textSizeLabelsRel,textSizeLabelsRel,kFALSE);
//     drawLatexAdd("e^{#pm} rec. with EMCal",0.15,0.92-2*textSizeLabelsRel,textSizeLabelsRel,kFALSE);
//

//__________________________________________________________________________________________________________
void SetStyleTLatex(
    TLatex* text,
    Size_t textSize,
    Width_t lineWidth,
    Color_t textColor = 1,
    Font_t textFont = 42,
    Bool_t kNDC = kTRUE,
    Short_t align = 11
){
    if (kNDC) {text->SetNDC();}
    text->SetTextFont(textFont);
    text->SetTextColor(textColor);
    text->SetTextSize(textSize);
    text->SetLineWidth(lineWidth);
    text->SetTextAlign(align);
}

//__________________________________________________________________________________________________________

void drawLatexAdd(TString latextext, Double_t textcolumn, Double_t textrow, Double_t textSizePixel,Bool_t setFont = kFALSE, Bool_t setFont2 = kFALSE, Bool_t alignRight = kFALSE, Color_t textcolor = kBlack, Float_t angle = 0., Bool_t alignBottomRight = kFALSE){
    TLatex *latexDummy                  = new TLatex(textcolumn ,textrow,latextext);
    SetStyleTLatex( latexDummy, textSizePixel,4, textcolor);
    if(setFont)
        latexDummy->SetTextFont(62);
    if(setFont2)
        latexDummy->SetTextFont(43);
    if(alignRight)
        latexDummy->SetTextAlign(31);
    if (alignBottomRight)
        latexDummy->SetTextAlign(33);
    if (angle != 0.)
        latexDummy->SetTextAngle(angle);
    latexDummy->SetTextColor(textcolor);
    latexDummy->Draw();
}


//__________________________________________________________________________________________________________
void DrawGammaSetMarkerTF1(
    TF1* fit1,
    Style_t lineStyle,
    Size_t lineWidth,
    Color_t lineColor
) {
    fit1->SetLineColor(lineColor);
    fit1->SetLineStyle(lineStyle);
    fit1->SetLineWidth(lineWidth);
}

//__________________________________________________________________________________________________________
void DrawGammaSetMarkerTSpline(
    TSpline* fit1,
    Style_t lineStyle,
    Size_t lineWidth,
    Color_t lineColor
) {
    fit1->SetLineColor(lineColor);
    fit1->SetLineStyle(lineStyle);
    fit1->SetLineWidth(lineWidth);
}


//__________________________________________________________________________________________________________

TLegend *GetAndSetLegend2(
    Double_t positionX,
    Double_t positionY,
    Double_t positionXRight,
    Double_t positionYUp,
    Size_t textSize,
    Int_t columns               = 1,
    TString header              = "",
    Font_t textFont             = 43,
    Double_t margin             = 0
){
    TLegend *legend = new TLegend(positionX,positionY,positionXRight,positionYUp);
    legend->SetNColumns(columns);
    legend->SetLineColor(0);
    legend->SetLineWidth(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetLineStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(textFont);
    legend->SetTextSize(textSize);
    if (margin != 0) legend->SetMargin(margin);
    if (header.CompareTo("")!= 0) legend->SetHeader(header);
    return legend;
}


//__________________________________________________________________________________________________________
void ReturnCorrectValuesForCanvasScaling(   Int_t sizeX,
                                            Int_t sizeY,
                                            Int_t nCols,
                                            Int_t nRows,
                                            Double_t leftMargin,
                                            Double_t rightMargin,
                                            Double_t upperMargin,
                                            Double_t lowerMargin,
                                            Double_t* arrayBoundariesX,
                                            Double_t* arrayBoundariesY,
                                            Double_t* relativeMarginsX,
                                            Double_t* relativeMarginsY,
                                            Bool_t verbose = kTRUE){
    Int_t realsizeX             = sizeX- (Int_t)(sizeX*leftMargin)- (Int_t)(sizeX*rightMargin);
    Int_t realsizeY             = sizeY- (Int_t)(sizeY*upperMargin)- (Int_t)(sizeY*lowerMargin);

    Int_t nPixelsLeftColumn     = (Int_t)(sizeX*leftMargin);
    Int_t nPixelsRightColumn    = (Int_t)(sizeX*rightMargin);
    Int_t nPixelsUpperColumn    = (Int_t)(sizeY*upperMargin);
    Int_t nPixelsLowerColumn    = (Int_t)(sizeY*lowerMargin);

    Int_t nPixelsSinglePlotX    = (Int_t) (realsizeX/nCols);
    Int_t nPixelsSinglePlotY    = (Int_t) (realsizeY/nRows);
    if(verbose){
        cout << realsizeX << "\t" << nPixelsSinglePlotX << endl;
        cout << realsizeY << "\t" << nPixelsSinglePlotY << endl;
        cout << nPixelsLeftColumn << "\t" << nPixelsRightColumn  << "\t" << nPixelsLowerColumn << "\t" << nPixelsUpperColumn << endl;
    }
    Int_t pixel = 0;
    if(verbose)cout << "boundaries X" << endl;
    for (Int_t i = 0; i < nCols+1; i++){
        if (i == 0){
            arrayBoundariesX[i] = 0.;
            pixel = pixel+nPixelsLeftColumn+nPixelsSinglePlotX;
        } else if (i == nCols){
            arrayBoundariesX[i] = 1.;
            pixel = pixel+nPixelsRightColumn;
        } else {
            arrayBoundariesX[i] = (Double_t)pixel/sizeX;
            pixel = pixel+nPixelsSinglePlotX;
        }
        if(verbose)cout << i << "\t" << arrayBoundariesX[i] << "\t" << pixel<<endl;
    }

    if(verbose)cout << "boundaries Y" << endl;
    pixel = sizeY;
    for (Int_t i = 0; i < nRows+1; i++){
        if (i == 0){
            arrayBoundariesY[i] = 1.;
            pixel = pixel-nPixelsUpperColumn-nPixelsSinglePlotY;
        } else if (i == nRows){
            arrayBoundariesY[i] = 0.;
            pixel = pixel-nPixelsLowerColumn;
        } else {
            arrayBoundariesY[i] = (Double_t)pixel/sizeY;
            pixel = pixel-nPixelsSinglePlotY;
        }
        if(verbose)cout << i << "\t" << arrayBoundariesY[i] <<"\t" << pixel<<endl;
    }

    relativeMarginsX[0]         = (Double_t)nPixelsLeftColumn/(nPixelsLeftColumn+nPixelsSinglePlotX);
    relativeMarginsX[1]         = 0;
    relativeMarginsX[2]         = (Double_t)nPixelsRightColumn/(nPixelsRightColumn+nPixelsSinglePlotX);;

    relativeMarginsY[0]         = (Double_t)nPixelsUpperColumn/(nPixelsUpperColumn+nPixelsSinglePlotY);
    relativeMarginsY[1]         = 0;
    relativeMarginsY[2]         = (Double_t)nPixelsLowerColumn/(nPixelsLowerColumn+nPixelsSinglePlotY);;

    return;
}

//__________________________________________________________________________________________________________
void ReturnCorrectValuesTextSize(   TPad * pad,
                                    Double_t &textsizeLabels,
                                    Double_t &textsizeFac,
                                    Int_t textSizeLabelsPixel,
                                    Double_t margin){

    textsizeLabels = 0;
    textsizeFac = 0;
    if (pad->XtoPixel(pad->GetX2()) < pad->YtoPixel(pad->GetY1())){
        textsizeLabels = (Double_t)textSizeLabelsPixel/pad->XtoPixel(pad->GetX2()) ;
        textsizeFac = (Double_t)1./pad->XtoPixel(pad->GetX2()) ;
    } else {
        textsizeLabels = (Double_t)textSizeLabelsPixel/pad->YtoPixel(pad->GetY1());
        textsizeFac = (Double_t)1./pad->YtoPixel(pad->GetY1());
    }
    cout << textsizeLabels << endl;
    cout << textsizeFac << endl;

    return;

}

//__________________________________________________________________________________________________________
Double_t ReturnTextSize(   TPad * pad,
                                        Int_t textSizeLabelsPixel
                                    ){

    Double_t textsizeLabels = 0;
    if (pad->XtoPixel(pad->GetX2()) < pad->YtoPixel(pad->GetY1())){
        textsizeLabels = (Double_t)textSizeLabelsPixel/pad->XtoPixel(pad->GetX2()) ;
    } else {
        textsizeLabels = (Double_t)textSizeLabelsPixel/pad->YtoPixel(pad->GetY1());
    }
    cout << "textsize: " << textsizeLabels << endl;
    return textsizeLabels;
    
}


TBox* CreateBox( Double_t xStart, Double_t yStart, Double_t xEnd, Double_t yEnd, Color_t color = kGray+2 ) {
    TBox* box = new TBox(xStart ,yStart , xEnd, yEnd);
    box->SetLineColor(color);
    box->SetFillColor(color);
    box->SetLineWidth(1);
    box->SetFillStyle(1001);
    return box;
}


//================================================================================================================
//Scale TF1 with constant
//================================================================================================================
TF1* ScaleTF1(TF1* func, Double_t constant, TString name) {

    if (!func) return NULL;

    Double_t    xMin, xMax;
    TString     formula         = func->GetExpFormula();
    func->GetRange(xMin, xMax);
        #if !defined (__CINT__) || defined (__CLING__)
        for (Int_t i=0; i<func->GetNpar(); i++) {
            formula.ReplaceAll(Form("[p%d]", i), Form("[placeholder%d]",i+1));
        }
        for (Int_t i=1; i<func->GetNpar()+1; i++) {
            formula.ReplaceAll(Form("[placeholder%d]", i), Form("[p%d]",i));
        }
    #else
        for (Int_t i=0; i<func->GetNpar(); i++) {
            formula.ReplaceAll(Form("[%d]", i), Form("[placeholder%d]",i+1));
        }
        for (Int_t i=1; i<func->GetNpar()+1; i++) {
            formula.ReplaceAll(Form("[placeholder%d]", i), Form("[%d]",i));
        }
    #endif

    TF1* result                 = new TF1(name.Data(), Form("[0] * (%s)", formula.Data()), xMin, xMax);
    for (Int_t i=0; i<func->GetNpar()+1; i++) {
        if (i==0)   result->SetParameter(i, constant);
        else        result->SetParameter(i, func->GetParameter(i-1));
    }

    return result;
}
