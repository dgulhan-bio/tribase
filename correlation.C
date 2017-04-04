#include "TH1D.h"
#include "TFile.h"
#include "parseTribaseFiles.h"

extern "C" const int chromSize[] = { 240, 240, 200, 190, 180,
				     170, 150, 140, 130, 130,
                                     130, 130, 110, 100, 100,
                                     90,  80,  70,  60,  60,
                                     50, 40, 150, 50};

#define nHistBin  40
#define totalMix  20

TH1D* calculateCF(vector<Individ> individs, int chromNum = 1){
  int nIndivid = individs.size();
  double length = ((double)chromSize[chromNum])*pow(10,6);
  cout << length << endl;
  TH1D * histCF = new TH1D(Form("histCF_chrom%d", chromNum), "", 40, 0.0001, 0.999999);
  for(int ind = 0; ind < nIndivid; ind++){
    int nMut = individs[ind].getNMuts();
    for(int iMut = 0; iMut < nMut; iMut++){
      Mutation * mut1 = individs[ind].getIthMut(iMut);
      if(mut1->getChrom() != chromNum) continue;
      for(int jMut = 0; jMut < nMut; jMut++){
        if(iMut == jMut) continue;
        Mutation * mut2 = individs[ind].getIthMut(jMut);
        if(mut2->getChrom() != chromNum) continue; 
        double deltaPos = fabs(((double)mut1->getPos()) - ((double) mut2->getPos()))/length;
        if(deltaPos==0) cout << deltaPos << endl;
        histCF->Fill(deltaPos);     
      }
    } 
  }
  histCF->Scale(1./2);
  return histCF;
}

TH1D* calculateMI(vector<Individ> individs, int chromNum = 1){
  int nIndivid = individs.size();
  double length = ((double)chromSize[chromNum])*pow(10,6);
  TH1D * histMI = new TH1D(Form("histMI_chrom%d", chromNum), "", 40, 0.00001, 0.99999);
  int nMI = 0; 
  for(int ind = 0; ind < nIndivid; ind++){
    int nMut1 = individs[ind].getNMuts();
    for(int iMut = 0; iMut < nMut1; iMut++){
      Mutation * mut1 = individs[ind].getIthMut(iMut);
      if(mut1->getChrom() != chromNum) continue;
      for(int indMI = 0; indMI < nIndivid; indMI++ ){
        if(indMI == ind) continue;
        int nMut2 = individs[indMI].getNMuts();
        for(int jMut = 0; jMut < nMut2; jMut++){
          Mutation * mut2 = individs[indMI].getIthMut(jMut);
          if(mut2->getChrom() != chromNum) continue; 
          //cout << "pos1 = " << mut1->getPos() << " pos2 =" << mut2->getPos() << endl;
	  double deltaPos = fabs(((double)mut1->getPos()) - ((double) mut2->getPos()))/length;
	  histMI->Fill(deltaPos);   
        }
      }
    }
  }
  histMI->Scale(1./((double)(nIndivid-1)));
  return histMI;
}

int main(int argc, char *argv[]){
  TH1D::SetDefaultSumw2();
  cout << argv[1] << endl; 
  string fileList = argv[1];
  Parser *parser = new Parser(fileList);
  vector<Individ> individs = parser->readIndivids();
  cout << "size of individs = " << individs.size() << endl; 
  TH1D* histsCF[22];
  TH1D* histsMI[22];
  for(int i = 1; i < 23; i++){
    histsCF[i-1] = calculateCF(individs, i);
    histsMI[i-1] = calculateMI(individs, i);
  }
  TFile * outf = new TFile("test.root", "recreate");
  for(int i = 1; i < 23; i++){
    histsCF[i-1]->Write();
    histsMI[i-1]->Write();
  }
  outf->Close();
  return 0;
}