#ifndef __FCmatrix__
#define __FCmatrix__

#include "Vec.h"
#include <vector>
#include <string>
#include <iostream>

using namespace std;

class FCmatrix {

 public:
  //constructors
  FCmatrix(); //default constructor
  FCmatrix(double** fM, int fm, int fn); //constructor, parameters explained e derived classes
  FCmatrix(double* fM, int fm, int fn); //constructor, parameters explained e derived classes
  FCmatrix(vector <Vec>); //constructor, directly from a vector<Vec>
  //destructor
  virtual ~FCmatrix() {M.clear();classname.clear();}

  //operators
  virtual Vec& operator[] (int) = 0;
  virtual FCmatrix& operator =(const FCmatrix &matrix);

  //methods
  string GetClassName() const;
  virtual void Print() = 0;
  virtual int n_rows()const = 0;
  virtual int n_cols()const = 0;
  virtual Vec GetRow(int) const = 0;
  virtual Vec GetCol(int) const = 0;
  virtual double Determinant() = 0;
  virtual int GetRowMax(int i=0) = 0;
  virtual int GetColMax(int i=0) = 0;
 protected:
  vector<Vec> M;
  string classname;
};


#endif
