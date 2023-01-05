#ifndef __FCmatrixBd__
#define __FCmatrixBd__

#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include <vector>
#include "Vec.h"
#include <iostream>

using namespace std;

class FCmatrixBanded : public FCmatrix{
 public:
  //constructorsa
  FCmatrixBanded(); //default constructor
  FCmatrixBanded(double **, int fm, int fn); //constructor from 3 vectors of dim (fm must be 3! fn must be the lenght of dim and all vectors must have size fn)
  FCmatrixBanded(double * fM, int fm, int fn); //constructor, fm = 3, fn = dimension of matrix, dimension of fM = (3*dimension of matrix) with 0 on postion fn-1 and 3*fn-1
  FCmatrixBanded(vector<Vec>); //constructor, directly from a vector<Vec>

  //copy constructor
  FCmatrixBanded(const FCmatrixBanded &); //create matrix from matrix

  //destructor
  ~FCmatrixBanded(){;}

  //converters
  FCmatrixFull toFull();

  //operators
  Vec& operator[] (int i){return M[i];} //retrieve row i
  FCmatrixBanded& operator=(const FCmatrix&){return *this;} //equivalent to copy constructor

  void Print(); //print this matrix
  int n_rows()const;
  int n_cols()const;
  Vec GetRow(int) const{return M[0];}
  Vec GetCol(int) const{return M[0];}
  double Determinant();
  int GetRowMax(int i=0){return 0;}
  int GetColMax(int i=0){return 0;}
  Vec GetDiagonal(int);


};
#endif
