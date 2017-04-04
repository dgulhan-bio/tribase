#include "parseTribaseFiles.h"

int main(int argc, char *argv[]){
  cout << argv[1] << endl; 
  string fileList = argv[1];
  Parser *parser = new Parser(fileList);
  vector<Individ> individs = parser->readIndivids();
  cout << "size of individs = " << individs.size();
  individs[0]
  return 0;
}
