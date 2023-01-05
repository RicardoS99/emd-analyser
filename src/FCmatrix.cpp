#include "FCmatrix.h"
#include "Vec.h"
#include "FCmatrixFull.h"
#include <vector>
#include <string>
#include <iostream>

using namespace std;
FCmatrix::FCmatrix(){;}

FCmatrix::FCmatrix(double **fM, int fm, int fn)
{
  for(int i=0; i<fm ;i++)
  {
    Vec temp(fn, fM[i]);
    M.push_back(temp);
  }
}

FCmatrix::FCmatrix (double * fM, int fm, int fn)
{
  for(int i=0; i<fm; i++)
  {
    double * line = new double[fn];
    for(int j=0; j<fn; j++)
    {
      line[j]=fM[j+i*fn];
    }
    Vec temp(fn, line);
    M.push_back(temp);
    delete [] line;
  }
}


FCmatrix::FCmatrix (vector<Vec> fM)
{
  M = fM;
}

string FCmatrix::GetClassName() const
{
  return classname;
}

FCmatrix& FCmatrix::operator =(const FCmatrix &matrix)
{
  if(this==&matrix)
  {
    return *this;
  }
  if(matrix.GetClassName()=="FCmatrixFull")
  {
    M.clear();
    int n_rows = matrix.n_rows();
    FCmatrixFull *mtemp = (FCmatrixFull*) &matrix;
    Vec vtemp;
    for(int i = 0; i<n_rows; i++)
    {
      vtemp = mtemp->GetRow(i);
      M.push_back(vtemp);
    }
  }
  return *this;
}
