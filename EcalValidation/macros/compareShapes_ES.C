#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TPaveStats.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

void compareShapes_ES()
{
  gROOT->ProcessLine(".x ~/public/style.C");

  //    TFile* fD = TFile::Open("../test/mcRERECO_ES_DIGI_DIGIRECO_Default.root");
  //  TFile* fD = TFile::Open("../test/mcRERECO_ES_DIGI_DIGIRECO_m10BX.root");

  
  TFile* f0 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_0BX.root");
  TFile* fm1 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_m1BX.root");
  TFile* fm2 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_m2BX.root");
  TFile* fm3 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_m3BX.root");
  TFile* fm4 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_m4BX.root");
  TFile* fm5 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_m5BX.root");
  TFile* fm6 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_m6BX.root");
  TFile* fm7 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_m7BX.root");
  TFile* fp1 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_p1BX.root");
  TFile* fp2 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_p2BX.root");
  TFile* fp3 = TFile::Open("../test/HG/mcRERECO_ES_DIGI_DIGIRECO_p3BX.root");


  std::cout << " OK 1 " << std::endl;

  std::vector<string> typeSet;
  typeSet.push_back("0BXPu");
  typeSet.push_back("m1BXPu");
  typeSet.push_back("m2BXPu");
  typeSet.push_back("m3BXPu");
  typeSet.push_back("m4BXPu");
  typeSet.push_back("m5BXPu");
  typeSet.push_back("m6BXPu");
  typeSet.push_back("m7BXPu");
  typeSet.push_back("p1BXPu");
  typeSet.push_back("p2BXPu");
  typeSet.push_back("p3BXPu");

  std::cout << " OK 2 " << std::endl;

  //h_PulseShape_ES
  //h_PulseShape_ES_kGood

  TProfile** hES = new TProfile*[11];
  hES[0] = (TProfile*)f0->Get("ecalvalidationSamplesMC/h_PulseShape_ES");
  hES[1] = (TProfile*)fm1->Get("ecalvalidationSamplesMC/h_PulseShape_ES");
  hES[2] = (TProfile*)fm2->Get("ecalvalidationSamplesMC/h_PulseShape_ES");
  hES[3] = (TProfile*)fm3->Get("ecalvalidationSamplesMC/h_PulseShape_ES");
  hES[4] = (TProfile*)fm4->Get("ecalvalidationSamplesMC/h_PulseShape_ES");
  hES[5] = (TProfile*)fm5->Get("ecalvalidationSamplesMC/h_PulseShape_ES");
  hES[6] = (TProfile*)fm6->Get("ecalvalidationSamplesMC/h_PulseShape_ES");
  hES[7] = (TProfile*)fm7->Get("ecalvalidationSamplesMC/h_PulseShape_ES");
  hES[8] = (TProfile*)fp1->Get("ecalvalidationSamplesMC/h_PulseShape_ES");
  hES[9] = (TProfile*)fp2->Get("ecalvalidationSamplesMC/h_PulseShape_ES");
  hES[10] = (TProfile*)fp3->Get("ecalvalidationSamplesMC/h_PulseShape_ES");

  TProfile** hESkG = new TProfile*[11];
  hESkG[0] = (TProfile*)f0->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");
  hESkG[1] = (TProfile*)fm1->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");
  hESkG[2] = (TProfile*)fm2->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");
  hESkG[3] = (TProfile*)fm3->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");
  hESkG[4] = (TProfile*)fm4->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");
  hESkG[5] = (TProfile*)fm5->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");
  hESkG[6] = (TProfile*)fm6->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");
  hESkG[7] = (TProfile*)fm7->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");
  hESkG[8] = (TProfile*)fp1->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");
  hESkG[9] = (TProfile*)fp2->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");
  hESkG[10] = (TProfile*)fp3->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGood");


  TProfile** hESpedSub = new TProfile*[11];
  hESpedSub[0] = (TProfile*)f0->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");
  hESpedSub[1] = (TProfile*)fm1->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");
  hESpedSub[2] = (TProfile*)fm2->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");
  hESpedSub[3] = (TProfile*)fm3->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");
  hESpedSub[4] = (TProfile*)fm4->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");
  hESpedSub[5] = (TProfile*)fm5->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");
  hESpedSub[6] = (TProfile*)fm6->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");
  hESpedSub[7] = (TProfile*)fm7->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");
  hESpedSub[8] = (TProfile*)fp1->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");
  hESpedSub[9] = (TProfile*)fp2->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");
  hESpedSub[10] = (TProfile*)fp3->Get("ecalvalidationSamplesMC/h_PulseShape_ESpedSub");

  TProfile** hESkGpedSub = new TProfile*[11];
  hESkGpedSub[0] = (TProfile*)f0->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");
  hESkGpedSub[1] = (TProfile*)fm1->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");
  hESkGpedSub[2] = (TProfile*)fm2->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");
  hESkGpedSub[3] = (TProfile*)fm3->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");
  hESkGpedSub[4] = (TProfile*)fm4->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");
  hESkGpedSub[5] = (TProfile*)fm5->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");
  hESkGpedSub[6] = (TProfile*)fm6->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");
  hESkGpedSub[7] = (TProfile*)fm7->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");
  hESkGpedSub[8] = (TProfile*)fp1->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");
  hESkGpedSub[9] = (TProfile*)fp2->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");
  hESkGpedSub[10] = (TProfile*)fp3->Get("ecalvalidationSamplesMC/h_PulseShape_ES_kGoodpedSub");


  TProfile** hEB = new TProfile*[11];
  hEB[0] = (TProfile*)f0->Get("ecalvalidationSamplesMC/h_PulseShape_EB");
  hEB[1] = (TProfile*)fm1->Get("ecalvalidationSamplesMC/h_PulseShape_EB");
  hEB[2] = (TProfile*)fm2->Get("ecalvalidationSamplesMC/h_PulseShape_EB");
  hEB[3] = (TProfile*)fm3->Get("ecalvalidationSamplesMC/h_PulseShape_EB");
  hEB[4] = (TProfile*)fm4->Get("ecalvalidationSamplesMC/h_PulseShape_EB");
  hEB[5] = (TProfile*)fm5->Get("ecalvalidationSamplesMC/h_PulseShape_EB");
  hEB[6] = (TProfile*)fm6->Get("ecalvalidationSamplesMC/h_PulseShape_EB");
  hEB[7] = (TProfile*)fm7->Get("ecalvalidationSamplesMC/h_PulseShape_EB");
  hEB[8] = (TProfile*)fp1->Get("ecalvalidationSamplesMC/h_PulseShape_EB");
  hEB[9] = (TProfile*)fp2->Get("ecalvalidationSamplesMC/h_PulseShape_EB");
  hEB[10] = (TProfile*)fp3->Get("ecalvalidationSamplesMC/h_PulseShape_EB");



  TProfile** hEE = new TProfile*[11];
  hEE[0] = (TProfile*)f0->Get("ecalvalidationSamplesMC/h_PulseShape_EE");
  hEE[1] = (TProfile*)fm1->Get("ecalvalidationSamplesMC/h_PulseShape_EE");
  hEE[2] = (TProfile*)fm2->Get("ecalvalidationSamplesMC/h_PulseShape_EE");
  hEE[3] = (TProfile*)fm3->Get("ecalvalidationSamplesMC/h_PulseShape_EE");
  hEE[4] = (TProfile*)fm4->Get("ecalvalidationSamplesMC/h_PulseShape_EE");
  hEE[5] = (TProfile*)fm5->Get("ecalvalidationSamplesMC/h_PulseShape_EE");
  hEE[6] = (TProfile*)fm6->Get("ecalvalidationSamplesMC/h_PulseShape_EE");
  hEE[7] = (TProfile*)fm7->Get("ecalvalidationSamplesMC/h_PulseShape_EE");
  hEE[8] = (TProfile*)fp1->Get("ecalvalidationSamplesMC/h_PulseShape_EE");
  hEE[9] = (TProfile*)fp2->Get("ecalvalidationSamplesMC/h_PulseShape_EE");
  hEE[10] = (TProfile*)fp3->Get("ecalvalidationSamplesMC/h_PulseShape_EE");



  std::cout << " OK 3 " << std::endl;

  int lineWidthExt[11] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };
  int markerStyleExt[11] = {27, 20, 20, 20, 20, 20, 20, 20, 4, 4, 4};
  int colorsExt[11] = {kGreen, kBlack, kBlue, kBlue-7, kAzure+1, kCyan+1, kCyan+3, kGreen+2, kRed, kRed+1, kRed+2 };
  std::vector<int> lineWidth;
  std::vector<int> markerStyle;
  std::vector<int> colors;
  for(int posVec = 0; posVec<11; ++posVec){
    lineWidth.push_back(lineWidthExt[posVec]);
    markerStyle.push_back(markerStyleExt[posVec]);
    colors.push_back(colorsExt[posVec]);
  }

  std::cout << " OK 4 " << std::endl;

  for(int i=0; i<11; ++i){
    //    std::cout << " name = " << typeSet.at(i) << std::endl;
    hEB[i]->SetName((typeSet.at(i)).c_str());
    hEB[i]->SetDirectory(0);
    hEB[i]->SetLineColor(colors.at(i));
    hEB[i]->SetMarkerColor(colors.at(i));
    hEB[i]->SetLineWidth(lineWidth.at(i));
    hEB[i]->SetMarkerStyle(markerStyle.at(i));

    hEE[i]->SetName((typeSet.at(i)).c_str());
    hEE[i]->SetDirectory(0);
    hEE[i]->SetLineColor(colors.at(i));
    hEE[i]->SetMarkerColor(colors.at(i));
    hEE[i]->SetLineWidth(lineWidth.at(i));
    hEE[i]->SetMarkerStyle(markerStyle.at(i));

    hESpedSub[i]->SetName((typeSet.at(i)).c_str());
    hESpedSub[i]->SetDirectory(0);
    hESpedSub[i]->SetLineColor(colors.at(i));
    hESpedSub[i]->SetMarkerColor(colors.at(i));
    hESpedSub[i]->SetLineWidth(lineWidth.at(i));
    hESpedSub[i]->SetMarkerStyle(markerStyle.at(i));

    hESkGpedSub[i]->SetName((typeSet.at(i)).c_str());
    hESkGpedSub[i]->SetDirectory(0);
    hESkGpedSub[i]->SetLineColor(colors.at(i));
    hESkGpedSub[i]->SetMarkerColor(colors.at(i));
    hESkGpedSub[i]->SetLineWidth(lineWidth.at(i));
    hESkGpedSub[i]->SetMarkerStyle(markerStyle.at(i));
    
    hES[i]->SetName((typeSet.at(i)).c_str());
    hES[i]->SetDirectory(0);
    hES[i]->SetLineColor(colors.at(i));
    hES[i]->SetLineWidth(lineWidth.at(i));
    hES[i]->SetMarkerColor(colors.at(i));
    hES[i]->SetMarkerStyle(markerStyle.at(i));
    
    hESkG[i]->SetName((typeSet.at(i)).c_str());
    hESkG[i]->SetDirectory(0);
    hESkG[i]->SetLineColor(colors.at(i));
    hESkG[i]->SetLineWidth(lineWidth.at(i));
    hESkG[i]->SetMarkerColor(colors.at(i));
    hESkG[i]->SetMarkerStyle(markerStyle.at(i));
    
  }

  std::cout << " OK 5 " << std::endl;
  f0->Close();
  fm1->Close();
  fm2->Close();
  fm3->Close();
  fm4->Close();
  fm5->Close();
  fm6->Close();
  fm7->Close();
  fp1->Close();
  fp2->Close();
  fp3->Close();

  std::cout << " OK 6 " << std::endl;

  TLegend *leg = new TLegend(0.70,0.80,0.99,1.,NULL,"brNDC");
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  for(int ItSet=0; ItSet<11; ++ItSet)
    leg->AddEntry(hEB[ItSet], (typeSet.at(ItSet)).c_str(), "l");

  std::cout << " OK 7 " << std::endl;

  TCanvas* cEB = new TCanvas();
  cEB->cd();
  hEB[0]->GetYaxis()->SetRangeUser(180, 220);
  hEB[0]->GetXaxis()->SetRangeUser(-1, 11);
  hEB[0]->GetXaxis()->SetTitle("EB pulse");
  hEB[0]->Draw("");
  for(int i=1; i<11; ++i){
    hEB[i]->Draw("same");
  }
  leg->Draw("same");

  TCanvas* cEE = new TCanvas();
  cEE->cd();
  hEE[0]->GetYaxis()->SetRangeUser(180, 220);
  hEE[0]->GetXaxis()->SetRangeUser(-1, 11);
  hEE[0]->GetXaxis()->SetTitle("EE pulse");
  hEE[0]->Draw("");
  for(int i=1; i<11; ++i){
    hEE[i]->Draw("same");
  }
  leg->Draw("same");

  TCanvas* cESpedSub = new TCanvas();
  cESpedSub->cd();
  hESpedSub[0]->GetYaxis()->SetRangeUser(0., 20.);
  hESpedSub[0]->GetXaxis()->SetRangeUser(-1, 11);
  hESpedSub[0]->GetXaxis()->SetTitle("ESpedSub pulse");
  hESpedSub[0]->Draw("");
  for(int i=1; i<11; ++i){
    hESpedSub[i]->Draw("same");
  }
  leg->Draw("same");

  TCanvas* cESkGpedSub = new TCanvas();
  cESkGpedSub->cd();
  hESkGpedSub[0]->GetYaxis()->SetRangeUser(0., 20.);
  hESkGpedSub[0]->GetXaxis()->SetRangeUser(-1, 11);
  hESkGpedSub[0]->GetXaxis()->SetTitle("ESkGpedSub pulse");
  hESkGpedSub[0]->Draw("");
  for(int i=1; i<11; ++i){
    hESkGpedSub[i]->Draw("same");
  }
  leg->Draw("same");

  TCanvas* cES = new TCanvas();
  cES->cd();
  hES[0]->GetYaxis()->SetRangeUser(0., 20.);
  hES[0]->GetXaxis()->SetRangeUser(-1, 11);
  hES[0]->GetXaxis()->SetTitle("ES pulse");
  hES[0]->Draw("");
  for(int i=1; i<11; ++i){
    hES[i]->Draw("same");
  }
  leg->Draw("same");

  TCanvas* cESkG = new TCanvas();
  cESkG->cd();
  hESkG[0]->GetYaxis()->SetRangeUser(0., 35.);
  hESkG[0]->GetXaxis()->SetRangeUser(-1, 11);
  hESkG[0]->GetXaxis()->SetTitle("ESkG pulse");
  hESkG[0]->Draw("");
  for(int i=1; i<11; ++i){
    hESkG[i]->Draw("same");
  }
  leg->Draw("same");

  std::cout << " OK fine " << std::endl;

}
