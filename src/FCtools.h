#ifndef _FCTOOLS_
#define _FCTOOLS_

#include <iostream>
#include "Vec.h"
#include <vector>
#include <string>
#include "TGraph.h"
using namespace std;


class FCtools{
 public:
  static vector<string> StrReadFile (string);
  static vector<vector<double> > VecReadFile(string);
};

#endif
