#include "parseTribaseFiles.h"
#include "TTree.h"
#include "TFile.h"

struct IndividFlat{
public:
  int nMut;
  vector<int> pos;
  vector<int> baseIn;
  vector<int> baseOut;
  vector<int> tribase;
  vector<int> chrom;
  vector<int> strand;
  IndividFlat(){
  }
  void set(Individ *individ){
    nMut = individ->getNMuts();
    cout << nMut << endl;
    for(int i = 0; i < nMut; i++){
      Mutation * mut = individ->getIthMut(i);
      pos.push_back(mut->getPos());
      baseIn.push_back(mut->getBaseIn());
      baseOut.push_back(mut->getBaseOut());
      strand.push_back(mut->getStrand());
      tribase.push_back(mut->getTribase());
      chrom.push_back(mut->getChrom());
      cout << "pos = " << pos[i] << "baseIn=" << baseIn[i] << "baseOut=" << baseOut[i] << "strand=" << strand[i] << "tribase=" << tribase[i] << "chrom=" << chrom[i] << endl;
    }
  }
  void reset(){
    pos.clear();
    baseIn.clear();
    baseOut.clear();
    strand.clear();
    tribase.clear();
    chrom.clear();
    nMut = 0;
  }
};

int main(int argc, char *argv[]){
  cout << argv[1] << endl;
  string fileList = argv[1];
  Parser *parser = new Parser(fileList);
  vector<Individ> individs = parser->readIndivids();
  cout << "size of individs = " << individs.size() << endl;

  IndividFlat individ;

  TFile * outf = new TFile("outfIndividuals.root","recreate");
  TTree * indTree = new TTree("indTree","tree of individuals and their snvs");
  cout << 1 << endl;
  indTree->Branch("nMut", &(individ.nMut),"nMut/I");
  indTree->Branch("pos", &(individ.pos));
  indTree->Branch("baseIn", &(individ.baseIn));
  indTree->Branch("baseOut", &(individ.baseOut));
  indTree->Branch("chrom",&(individ.chrom));
  indTree->Branch("strand",&(individ.strand));
  indTree->Branch("tribase",&(individ.tribase));
  
  for(int i = 0; i < individs.size(); i++){
    individ.reset();
    cout << "nMut" << individs[i].getNMuts() << endl;
    individ.set(&individs[i]);
    indTree->Fill();
    cout << 4 << endl;    
  }
  cout << 5 << endl;
  indTree->Write();
  outf->Close();

  return 0;
}
