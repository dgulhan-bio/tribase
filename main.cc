#include "parseTribaseFiles.h"

int main(int argc, char *argv[]){
  cout << argv[1] << endl; 
  string fileList = argv[1];
   Parser *parser = new Parser(fileList);
  return 0;
}
