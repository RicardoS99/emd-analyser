#ifndef __DataPoints__
#define __DataPoints__

#include "cFCgraphics.h"
#include <string>
#include <vector>
using namespace std;

class DataPoints{
public:
  DataPoints();
  DataPoints(int N, double*xpoints, double*ypoints);
  DataPoints(vector<double> xpoints, vector<double> ypoints);
  virtual ~DataPoints();

  virtual void Print();

protected:
  int N; // number of data points
  double * x,* y; //arrays
};

#endif
