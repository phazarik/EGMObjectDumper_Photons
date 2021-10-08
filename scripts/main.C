#include "photonTreeProducer.C"
#include "TROOT.h"
#include <algorithm>


int main(int argc, const char *argv[]){
  gROOT->ProcessLine("#include <vector>");
  //gROOT->ProcessLine("#include <map>");

  string outfile = argv[2];
  long int nTotEvt = atof(argv[3]);
  long int nPrintEvt = atof(argv[4]);
  float evtWt = atof(argv[5]);

  photonTreeProducer m(argv[1]);
  m.Loop(outfile,nTotEvt,nPrintEvt,evtWt);
}
