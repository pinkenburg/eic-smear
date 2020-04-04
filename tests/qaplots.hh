#ifndef QAPLOTS_H
#define QAPLOTS_H

#include <TString.h>

#include <string>
#include <vector>

struct qaparameters{
  // std::string txtfilename="./tests/beagle_eD.txt";
  std::string txtfilename="./tests/ep_noradcorr.20x250.small.txt";
  TString outfilebase="./qaplots";
  long nevents=-1;
  std::vector<int> pids = {11 , 211, 2112, 2212}; // sign will be ignored. leave empty for everything. 

  std::string detstring = "BeAST"; // Capitalization does not matter
};

qaparameters ParseArguments ( int argc, char* argv[] );

#endif // QAPLOTS_H
