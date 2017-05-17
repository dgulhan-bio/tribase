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
#include "parseTribaseFiles.h"



void parse_reference_genome(string fileList){
  
  //parse snvs
  Parser *parser = new Parser(fileList);
  vector<Individ> individs;
  individs = parser->readIndividsFromTribase();
  
 
  int nIndivid = individs.size();  
  for(int ind = 0; ind < nIndivid; ind++){
    
    int nMut = individs[ind].getNMuts();
    int iMut = 0;
    Mutation *currentMut = individs[ind].getIthMut(iMut);
    int searchChrom = currentMut -> getChrom();
    int searchPos = currentMut -> getPos();
    int currentPos = 0; 
    int prevPos = 0;
    int lineAfter = 2;

    //parse reference
    ifstream inputFile("cleaned.fasta");
    string lineTmp;
    int chromNum = 0;
    int sizeChrom = 0;
    int state = 0;
    string stringMem1;
    string stringMem2;
    string composed;
    while (getline(inputFile, lineTmp)){
      trim(lineTmp);
      if (lineTmp.length() == 0) continue;
      if (string::npos == lineTmp.find("chromosome")){
        state = 1;
        string junk1, junk2, junk3, junk4, junk5, junk6, junk7;
        stringstream lineStream;
        lineStream << lineTmp;
        lineStream >> junk1 >> chromNum >> junk2 >> junk3 >> junk4 >> junk5 >> junk6 >> sizeChrom >> junk7;        
        currentPos = 0;
        prevPos = 0;
        getline(inputFile, lineTmp);
      }
      if (searchChrom < chromNum && state == 1){
        getline(inputFile, lineTmp);       
      }
      else if (searchChrom == chromNum && state == 1 ){
        state = 2;
      }
      if (state == 2){
        string dummyString;
        stringstream lineStream;
        lineStream << lineTmp;
        lineStream >> dummyString;
        currentPos += dummyString.size();
        
        while (currentPos >= searchPos){
          char base = dummyString.at(searchPos - prevPos);
          if (getBaseDigit(base) != currentMut->getBaseIn()){
            cout << "Error parser read doesn't match base in" << endl;
            return;
	  }
          string baseOutStr(1,  convertDigitToBase(currentMut->getBaseOut()));
	  dummyString.replace(searchPos - prevPos, 1, baseOutStr);
	  
	  if(lineAfter >= 2){
            composed = stringMem1 + stringMem2 + dummyString;
	  }else{
            composed = composed + dummyString;
          }

          iMut++;
          currentMut = individs[ind].getIthMut(iMut);
          searchChrom = currentMut -> getChrom();
          searchPos = currentMut -> getPos();          
          if(currentPos < searchPos) lineAfter = 0; 
        }
        
        if(lineAfter < 2){
	  composed = composed + dummyString;
        } 

        stringMem2 = stringMem1;
        stringMem1 = dummyString;
        
        prevPos = currentPos;
        lineAfter++;
      }
    }
  }
}

int main(int argc, char *argv[]){
  string fileList = argv[1];
  parse_reference_genome(fileList);
}
