#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

using namespace std;

struct Mutation{
private:
  int baseIn_;
  int baseOut_;
  int strand_; 
  int pos_;
  int chrom_;
  int tribase_;  
public:  
  Mutation(int baseIn, int baseOut,  int pos, int chrom, int strand = 0, int tribase = 1000){
    baseIn_ = baseIn;
    baseOut_ = baseOut;
    strand_ = strand;
    pos_ = pos;
    chrom_ = chrom;
    tribase_ = tribase;
  }
  int getBaseIn(){ return baseIn_; }
  int getBaseOut(){ return baseOut_; }
  int getStrand(){ return strand_; }
  int getPos(){ return pos_; }
  int getChrom(){ return chrom_; }
  int getTribase(){ return tribase_; }
};

struct Individ{
private:
  vector<Mutation> muts_;
  int age_; 
  bool isFemale_; 
  bool isPatient_;
public:
  Individ(){
    age_ = 0; 
    isFemale_ = false;
    isPatient_ = false;
  }
  void addMut(Mutation mut){
    muts_.push_back(mut);
  }
  void setMuts(vector<Mutation> muts){ 
    for(int mutInd = 0; mutInd < muts.size(); mutInd++){
      muts_.push_back(muts[mutInd]);
    }
  }
  void setAge(int age){ age_ = age; }
  void setSex(bool isFemale){ isFemale_ = isFemale; }
  void setPatient(bool isPatient){ isPatient_ = isPatient; }
  int getAge(){ return age_; }
  int isFemale(){ return isFemale_; }
  int isPatient(){ return isPatient_; }
  int getNMuts(){ return muts_.size(); }
  Mutation * getIthMut(int index){ return &muts_[index]; }
};

