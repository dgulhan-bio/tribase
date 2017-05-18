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

   individs = parser->readIndividsSnvMnvVcf();
   //individs = parser->readIndividsFromTribase();
 
  int nIndivid = individs.size();  
  cout << nIndivid << endl;
  
  for(int ind = 0; ind < nIndivid; ind++){
    
    int nMut = individs[ind].getNMuts();
    cout << nMut << endl;
    int iMut = 0;
    Mutation *currentMut = individs[ind].getIthMut(iMut);
     
    int searchChrom = currentMut -> getChrom();
    int searchPos = currentMut -> getPos();
    cout << "searchChrom = " << searchChrom << " searchPos = " <<  searchPos << endl;
    int currentPos = 0; 
    int prevPos = 0;
    int lineAfter = 2;

    //parse reference
    ifstream inputFileGenome("cleaned.fasta");
    if (!inputFileGenome.is_open()){
      cout << "file is not open" << endl;
      return;
    }

    string lineTmp;
    int chromNum = 0;
    int sizeChrom = 0;
    int state = 0;
    string stringMem1;
    string stringMem2;
    string composed;
    while (getline(inputFileGenome, lineTmp)){
      //cout << state << endl;
      if(iMut == nMut) break;
      trim(lineTmp);
      if (lineTmp.length() == 0) continue;
      if (string::npos != lineTmp.find("chromosome")){
        state = 1;
        string junk1, junk2, junk3, junk4, junk5, junk6, junk7, chromNumStr, sizeChromStr;
        stringstream lineStream;
        lineStream << lineTmp;
        lineStream >> chromNumStr >> junk1 >> junk2 >> junk3 >> junk4 >> junk5 >> junk6 >> sizeChromStr >> junk7;      
         
        if(chromNumStr == "X") chromNum = 23;
        else if(chromNumStr == "Y") chromNum = 24;
        else{
          chromNum = atoi(chromNumStr.c_str());
        }
        cout << "state = " << 1 << " searchChrom = " << searchChrom << " chromNum = " << chromNum << endl;
        sizeChrom = atoi(sizeChromStr.c_str());
        currentPos = 0;
        prevPos = 0;
        getline(inputFileGenome, lineTmp);
      }
      if (searchChrom < chromNum && state == 1){
        getline(inputFileGenome, lineTmp);       
      }
      else if (searchChrom == chromNum && state == 1 ){
        cout << "state = " << 2 << " searchChrom = " << searchChrom << " chromNum = " << chromNum << endl;   
        state = 2;
      }
      if (state == 2){
	float percentPos = (100*((float)currentPos)/((float)searchPos));
        //if(percentPos - (int)percentPos == 0) cout << "currentPos/searchPos % = " << percentPos << endl;
        string dummyString;
        stringstream lineStream;
        lineStream << lineTmp;
        lineStream >> dummyString;
        currentPos += dummyString.size();
        while (currentPos >= searchPos && searchChrom == chromNum){
	  char base = dummyString.at(searchPos - prevPos - 1);
          if (getBaseDigit(base) != currentMut->getBaseIn()){
            cout << "Error parser read doesn't match base in " << searchPos << endl;
            return;
	  }
  
          string baseOutStr(1,  convertDigitToBase(currentMut->getBaseOut()));
	  dummyString.replace(searchPos - prevPos - 1, 1, baseOutStr);
	  if(searchPos == 248551710) cout << "dummy = " << dummyString << endl;
          
	  if(lineAfter >= 2){
            composed = stringMem1 + stringMem2 + dummyString;
	  }else{
            composed = composed + dummyString;
          }
          
          iMut++;
          currentMut = individs[ind].getIthMut(iMut);
          
          searchChrom = currentMut -> getChrom();
          searchPos = currentMut -> getPos(); 
          
          if(currentPos < searchPos && searchChrom == chromNum) lineAfter = 0; 
          
        }
        
        if(lineAfter < 2){
	  composed = composed + dummyString;
        } 

        stringMem2 = stringMem1;
        stringMem1 = dummyString;
        
        prevPos = currentPos;
        if(lineAfter < 3) lineAfter++;
        if(lineAfter == 2){
          cout << "position: " << prevPos << endl;
          cout << composed << endl;
        }
        if(currentPos + 80 > sizeChrom ){
          state = 0;
          cout << "position: " << prevPos << endl;
          cout << composed << endl;
        }
      }
    }
  }
}

int main(int argc, char *argv[]){
  string fileList = argv[1];
  parse_reference_genome(fileList);
}
