#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "utilities.h"

void plotCorrelation(){
  TH1D::SetDefaultSumw2();
  TFile * file[22];
  for(int iChrom = 0; iChrom < 22; iChrom++){
    if(iChrom<4) file[iChrom] = new TFile(Form("testSameBaseInOnly%d-Mix20-40.root",iChrom+1));
    else file[iChrom] = new TFile(Form("testSameBaseInOnly%d.root",iChrom+1));
  }
  TH1D * hist_MI[22];
  TH1D * hist_CF[22];
  for(int iChrom = 0; iChrom < 22; iChrom++){
    cout << 1 << endl;
    hist_MI[iChrom] = (TH1D*) file[iChrom]->Get(Form("histMI_chrom%d",iChrom+1));
    hist_CF[iChrom] = (TH1D*) file[iChrom]->Get(Form("histCF_chrom%d",iChrom+1));
    cout << 2 << endl;
  } 
 
  TH1D* hist_rat[22];
  TCanvas * canvas = new TCanvas("canvas","", 6*300, 4*300);
  cout << 3 << endl;
  makeMultiPanelCanvas(canvas,6,4,0.0,0.0,0.25,0.25,0.02);
  cout << 4 << endl;
  //canvas->Divide(6,4);
  TLine line(0.,1.,1.,1.);
  line.SetLineWidth(2);
  for(int iChrom = 0; iChrom < 22; iChrom++){
    hist_rat[iChrom] = (TH1D*)hist_CF[iChrom]->Clone(Form("hist_rat_chrom%d",iChrom));
    hist_rat[iChrom]->Divide(hist_MI[iChrom]);
    cout << hist_rat[iChrom]->Integral() << endl;
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
    hist_rat[iChrom]->SetMaximum(1.1999);
    hist_rat[iChrom]->SetMinimum(0.850001);
    hist_rat[iChrom]->Draw();
    line.Draw("same");
    if(iChrom == 4 || iChrom == 10) drawText(Form("CHROM %d",iChrom+1),0.55,0.08);
    else if(iChrom == 16) drawText(Form("CHROM %d",iChrom+1),0.35,0.28);
    else if(iChrom > 16)  drawText(Form("CHROM %d",iChrom+1),0.15,0.28);
    else  drawText(Form("CHROM %d",iChrom+1),0.15,0.08);
    canvas->cd(iChrom+1)->RedrawAxis();
  } 
  cout << 5 << endl;
  canvas->SaveAs("test.pdf");
  cout << 6 << endl;

  TCanvas * canvasCI = new TCanvas("canvasCI","", 6*300, 4*300);
  cout << 3 << endl;
  makeMultiPanelCanvas(canvasCI,6,4,0.0,0.0,0.25,0.25,0.02);
  for(int iChrom = 0; iChrom < 22; iChrom++){
    canvasCI->cd(iChrom+1);
    hist_CF[iChrom]->Draw();
    hist_MI[iChrom]->Draw("same hist");
  }
  canvasCI->SaveAs("testDist.pdf");
 
}
