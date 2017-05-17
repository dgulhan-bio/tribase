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
  cout << "nIndivid =" << nIndivid << endl;
  double nPair = 0;
  TH1D * histCF = new TH1D(Form("histCF_chrom%d", chromNum), "", 40, 0.00001, 0.999999);
  for(int ind = 0; ind < nIndivid; ind++){
    cout << "ind=" << ind << endl;
    int nMut = individs[ind].getNMuts();
    for(int iMut = 0; iMut < nMut; iMut++){
      Mutation * mut1 = individs[ind].getIthMut(iMut);
      if(mut1->getChrom() != chromNum) continue;
      for(int jMut = 0; jMut < nMut; jMut++){
        if(iMut == jMut) continue;
        Mutation * mut2 = individs[ind].getIthMut(jMut);
        if(mut2->getChrom() != chromNum) continue; 
        //if(!(mut1->getBaseIn() == mut2->getBaseIn() && mut1->getStrand() == mut2->getStrand())) continue;
        double deltaPos = fabs(((double)mut1->getPos()) - ((double) mut2->getPos()))/length;
        //if(deltaPos==0) cout << deltaPos << endl;
        histCF->Fill(deltaPos,0.01);
        nPair+=0.01;  
        //if((int)nPair%1000==0) cout << "npair = "<< nPair << endl;
      }
    } 
  }
  histCF->Scale(1./nPair);
  cout << "CF npair = " << ((double)nPair) << endl; 
  return histCF;
}

TH1D* calculateMI(vector<Individ> individs, int chromNum = 1){
  int nIndivid = individs.size();
  double length = ((double)chromSize[chromNum])*pow(10,6);
  TH1D * histMI = new TH1D(Form("histMI_chrom%d", chromNum), "", 40, 0.00001, 0.999999);
  int nMI = 0; 
  double nPair = 0;
  for(int ind = 0; ind < nIndivid; ind++){
    int nMut1 = individs[ind].getNMuts();
    if(ind%100==0) cout << "ind = "<< ind << "/" << nIndivid << "nMut 1= " << nMut1 << endl;
    for(int iMut = 0; iMut < nMut1; iMut++){
      Mutation * mut1 = individs[ind].getIthMut(iMut);
      if(mut1->getChrom() != chromNum) continue;
      for(int indMI = 20; indMI < 40; indMI++ ){
        if(indMI == ind) continue;
        int nMut2 = individs[indMI].getNMuts();
        for(int jMut = 0; jMut < nMut2; jMut++){
          Mutation * mut2 = individs[indMI].getIthMut(jMut);
          if(mut2->getChrom() != chromNum) continue; 
          //cout << "pos1 = " << mut1->getPos() << " pos2 =" << mut2->getPos() << endl;
	  double deltaPos = fabs(((double)mut1->getPos()) - ((double) mut2->getPos()))/length;
	  histMI->Fill(deltaPos,0.01);
          nPair+=0.01;   
        }
      }
    }
  }
  histMI->Scale(1./nPair);
  cout << "MI npair = " << (double)nPair << endl; 
  return histMI;
}

int main(int argc, char *argv[]){
  TH1D::SetDefaultSumw2();
  cout << argv[1] << endl; 
  string fileList = argv[1];
  int isCancer = atoi(argv[2]);
  int chrom = atoi(argv[3]);
  Parser *parser = new Parser(fileList);
  vector<Individ> individs;

  if(!isCancer) individs = parser->readIndividsFromTribase();
  else individs = parser->readIndividsSnvMnvVcf();
  cout << "size of individs = " << individs.size() << endl; 
  TH1D* histsCF[22];
  TH1D* histsMI[22];
  for(int i = chrom; i < chrom+1; i++){
    //for(int i = 4; i < 5; i++){
    cout << "i = " << i << endl; 
    histsCF[i-1] = calculateCF(individs, i);
    cout << "calculated CF" << endl;  
    histsMI[i-1] = calculateMI(individs, i);
    cout << "calculated MI" << endl;
  }
  TFile * outf = new TFile(Form("testSameBaseInOnly%d-Mix20-40.root", chrom), "recreate");
  //for(int i = 4; i < 5; i++){
  for(int i = chrom; i < chrom+1; i++){
    histsCF[i-1]->Write();
    histsMI[i-1]->Write();
  }
  individs.clear();
  outf->Close();
  return 0;
}
