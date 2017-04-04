#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "utilities.h"

void plotCorrelation(){
  TH1D::SetDefaultSumw2();
  TFile * file = new TFile("test.root");
  TH1D * hist_MI[22];
  TH1D * hist_CF[22];
  for(int iChrom = 0; iChrom < 22; iChrom++){
    hist_MI[iChrom] = (TH1D*) file->Get(Form("histMI_chrom%d",iChrom+1));
    hist_CF[iChrom] = (TH1D*) file->Get(Form("histCF_chrom%d",iChrom+1));
  } 
 
  TH1D* hist_rat[22];
  TCanvas * canvas = new TCanvas("canvas","", 6*300, 4*300);
  makeMultiPanelCanvas(canvas,6,4,0.0,0.0,0.25,0.25,0.02);
  //canvas->Divide(6,4);
  for(int iChrom = 0; iChrom < 22; iChrom++){
    hist_rat[iChrom] = (TH1D*)hist_CF[iChrom]->Clone(Form("hist_rat_chrom%d",iChrom));
    hist_rat[iChrom]->Divide(hist_MI[iChrom]);
    canvas->cd(iChrom+3);
    hist_rat[iChrom]->GetXaxis()->SetTitle("#Delta(Position)/L_{Chrom}");
    hist_rat[iChrom]->GetYaxis()->SetTitle("CF / ME");
    hist_rat[iChrom]->GetYaxis()->CenterTitle();
    hist_rat[iChrom]->GetXaxis()->CenterTitle();
    makePretty(hist_rat[iChrom],3.4);
    hist_rat[iChrom]->SetMarkerSize(2);
    hist_rat[iChrom]->SetMarkerStyle(20);
    hist_rat[iChrom]->SetMarkerColor(kBlue-3);
    hist_rat[iChrom]->SetLineColor(kBlue-3);
    hist_rat[iChrom]->SetMaximum(0.999);
    hist_rat[iChrom]->SetMinimum(0.001);
    hist_rat[iChrom]->Draw();
    drawText(Form("CHROM %d",iChrom+1),0.25,0.28);
    canvas->cd(iChrom+1)->RedrawAxis();
  } 
  canvas->SaveAs("test.pdf");
}