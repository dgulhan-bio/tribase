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

const char replacement[][3] =
  {                 
    {'B', 'D', 'E'},
    {'F', 'H', 'I'},
    {'J', 'K', 'L'},
    {'M', 'N', 'O'}
  };
const char base_in_array[][3] =
  {                 
    {'A', 'A', 'A'},
    {'C', 'C', 'C'},
    {'G', 'G', 'G'},
    {'T', 'T', 'T'}
  };
const char base_out_array[][3] =
  {                 
    {'C', 'G', 'T'},
    {'A', 'G', 'T'},
    {'A', 'C', 'T'},
    {'A', 'C', 'G'}
  };

const char assign_alphabet(char base_in, char base_out){
  for(int row = 0; row < 4; row++){
    for(int column = 0; column < 3; column++){
      if(base_in == base_in_array[row][column] && base_out == base_out_array[row][column]){ 
        return replacement[row][column];
      }
    }
  }
  return 'Z';
}

void parse_reference_genome(string fileList){
 
  //parse snvs
  Parser *parser = new Parser(fileList);
  vector<Individ> individs;
  
   individs = parser->readIndividsSnvMnvVcf();
   //individs = parser->readIndividsFromTribase();
  
  int nIndivid = individs.size();  
  
  for(int ind = 0; ind < nIndivid; ind++){
    ofstream outputfile;
    std::ostringstream streamOutputFile;
    streamOutputFile << "output/" << parser->getIthFileName(ind) << ".txt";
    std::string outFileName = streamOutputFile.str();
    outputfile.open(outFileName);
    
    int nMut = individs[ind].getNMuts();
    //outputfile << "#Individ " << parser->getIthFileName(ind) << " totalMut" << nMut << "\n";
    int iMut = 0;
    Mutation *currentMut = individs[ind].getIthMut(iMut);
     
    int searchChrom = currentMut -> getChrom();
    int searchPos = currentMut -> getPos();
    //cout << "searchChrom = " << searchChrom << " searchPos = " <<  searchPos << endl;
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
        //cout << "state = " << 1 << " searchChrom = " << searchChrom << " chromNum = " << chromNum << endl;
        sizeChrom = atoi(sizeChromStr.c_str());
        currentPos = 0;
        prevPos = 0;
        getline(inputFileGenome, lineTmp);
      }
      if (searchChrom < chromNum && state == 1){
        getline(inputFileGenome, lineTmp);       
      }
      else if (searchChrom == chromNum && state == 1 ){
        outputfile << "#Chrom" << chromNum << "\n" ; 
        //cout << "state = " << 2 << " searchChrom = " << searchChrom << " chromNum = " << chromNum << endl;   
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
          
          char snv = assign_alphabet(base, convertDigitToBase(currentMut->getBaseOut()));
          string baseOutStr(1,  snv);
	  dummyString.replace(searchPos - prevPos - 1, 1, baseOutStr);
	  
	  if(lineAfter >= 2){
            composed = stringMem1 + stringMem2 + dummyString;
	  }
	  
          iMut++;
          currentMut = individs[ind].getIthMut(iMut);
          
          searchChrom = currentMut -> getChrom();
          searchPos = currentMut -> getPos(); 
          
          if(currentPos < searchPos && searchChrom == chromNum) lineAfter = 0; 
          
        }
        
        if(lineAfter <= 2 && lineAfter>0){
	  composed += dummyString;
        } 

        stringMem2 = stringMem1;
        stringMem1 = dummyString;
        
        prevPos = currentPos;
        if(lineAfter < 3) lineAfter++;
        if(lineAfter == 2){
          outputfile << "position: " << prevPos << "\n";
          outputfile << composed << "\n";
        }
        if(currentPos + 80 > sizeChrom ){
          state = 0;
          outputfile << "position: " << prevPos << "\n";
          outputfile << composed << "\n";
        }
      }
    }
  }
}

int main(int argc, char *argv[]){
  string fileList = argv[1];
  parse_reference_genome(fileList);
}

