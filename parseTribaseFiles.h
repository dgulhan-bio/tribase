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
#include "mutation.h"

#ifdef __cplusplus
extern "C" {
#endif
  extern const int array[];
#ifdef __cplusplus
}
#endif

using namespace std;
// trim from start
static inline string &ltrim(string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
  return s;
}

// trim from end
static inline std::string &rtrim(string &s) {
  s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
  return s;
}

// trim from both ends
static inline string &trim(string &s) {
  return ltrim(rtrim(s));
}

int getBaseDigit(char baseChar){
  switch(baseChar){
  case 'A': 
    return 0;
  case 'C': 
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    cout << "invalid base" << endl;
    return -1;
  }
}

char convertDigitToBase(int digit){
  switch(digit){
  case 0: 
    return 'A';
  case 1: 
    return 'C';
  case 2:
    return 'G';
  case 3:
    return 'T';
  default:
    cout << "invalid base" << endl;
    return -1;
  }
}

int getBaseNum(const char* baseChar){ 
  int baseLength = strlen(baseChar);
  if(baseLength == 0) return -1;
  else if(baseLength == 1) return getBaseDigit(baseChar[0]);
  else{
    int num = pow(10,baseLength);
    for(int i = 0; i < baseLength; i++){
      int digit = getBaseDigit(baseChar[i]);
      num += digit*pow(10, baseLength - i - 1);
    }
    return num;
  }
}
  
int getStrandNum(const char* strand){
  switch(strand[0]){
  case '+':
    return 1;
  case '-':
    return -1;
  default:
    return 0;
  }
}

struct Parser{
private:
  vector<Individ> individs_;
  vector<string> inFiles;

public:
  vector<Individ> getIndivids(){
    return individs_;
  }

  Parser(string inputFiles){
    //Load the file list in the input txt file

    ifstream fileList(inputFiles);
    if(!fileList){
      cout<<"Error opening file list"<<endl;
      system("pause");
      return;
    }
    string line;
    while(getline(fileList, line)){
      inFiles.push_back(line);
      cout << "file in:" << line << endl;
    }
  }
  
  vector<Individ> readIndividsFromTribase(){
    for(int fileInd = 0; fileInd < inFiles.size(); fileInd++){
      ifstream inputFile(inFiles[fileInd]);
      int state = 0;
      Individ current;
      string lineTmp;
      while(getline(inputFile, lineTmp)){
        trim(lineTmp);
        if(lineTmp.length() == 0) continue;
        if(state == 0){
	  state = 1;
          getline(inputFile, lineTmp);
        }
        if(state == 1){
          string chromStr; 
          int pos;
          string baseInStr;
          string baseOutStr;
          string tribaseStr;
          string dummy; 
          string strandStr;
          stringstream lineStream;
          lineStream << lineTmp;
          lineStream >> chromStr >> pos >> baseInStr >> baseOutStr >> tribaseStr >> dummy >> strandStr;
          //cout << chromStr <<  " " << pos << " " << baseInStr << " " << baseOutStr << " "<< tribaseStr <<  " " << strandStr << endl; 
       
          int chrom;
          if(chromStr == "X") chrom = 23;
          else chrom = atoi(chromStr.c_str());
          
          int baseIn = getBaseNum(baseInStr.c_str());
          int baseOut = getBaseNum(baseOutStr.c_str());
          int tribase = getBaseNum(tribaseStr.c_str());
          int strand = getStrandNum(strandStr.c_str());
          
	  Mutation mut(baseIn, baseOut, pos, chrom, strand, tribase);          
          current.addMut(mut);
        }
        //cout << "nMuts" << current.getNMuts() << endl;
      }
      individs_.push_back(current);
    }
    individs_.size();
    return individs_;
  }
  
  vector<Individ> readIndividsSnvMnvVcf(){
    for(int fileInd = 0; fileInd < inFiles.size(); fileInd++){
      ifstream inputFile(inFiles[fileInd]);
      int state = 0;
      Individ current;
      string lineTmp;
      while(getline(inputFile, lineTmp)){
        //cout << lineTmp << endl;
        trim(lineTmp);
        if(lineTmp.length() == 0) continue;
        if(string::npos == lineTmp.find("#CHROM")){
	  state = 1;
          getline(inputFile, lineTmp);
        }
        if(state == 1){
          string chromStr; 
          int pos;
          string baseInStr;
          string baseOutStr;
          string strandStr;
          string dummy1; 
          string dummy2; 
          string dummy3; 
          string dummy4; 
          stringstream lineStream;
          lineStream << lineTmp;
          lineStream >> chromStr >> pos >> dummy1 >> baseInStr >> baseOutStr >> dummy2 >> dummy3 >> dummy4;
          //cout << chromStr << " " << pos << " " << dummy1 << " " << " " <<  baseInStr << " " <<  baseOutStr << " " << dummy2 << " " << dummy3 << " " << dummy4 << endl;
       
          int chrom;
          if(chromStr == "X") chrom = 23;
          else chrom = atoi(chromStr.c_str());
          
          int baseIn = getBaseNum(baseInStr.c_str());
          int baseOut = getBaseNum(baseOutStr.c_str());
          
          
	  Mutation mut(baseIn, baseOut, pos, chrom);          
          current.addMut(mut);
          //cout << "nMut for loop=" << current.getNMuts() << endl;
        }
      }
      //cout << "nMuts after loop=" << current.getNMuts() << endl;
      individs_.push_back(current);
      //cout << "size individs = " << individs_.size() << endl;
    }
    return individs_;
  }
};
