#ifndef __FCmatrixF__
#define __FCmatrixF__

#include "FCmatrix.h"
#include <vector>
#include "Vec.h"
#include <iostream>

using namespace std;

class FCmatrixFull : public FCmatrix{
 public:
  //constructors
  FCmatrixFull(); //default constructor
  FCmatrixFull(double **, int , int ); //constructor, mx fm x fn
  FCmatrixFull(double * fM, int fm, int fn); //constructor, mx fm x fn (warning: LENGTH OF fM MUST EQUAL fm *fn !!!)
  FCmatrixFull(vector<Vec>); //constructor, directly from a vector<Vec>

  //copy constructor
  FCmatrixFull(const FCmatrixFull &); //create matrix from matrix

  //destructor
  ~FCmatrixFull();

  //operators
  Vec& operator[] (int); //retrieve row i
  FCmatrixFull& operator=(const FCmatrixFull&); //equivalent to copy constructor
  FCmatrixFull operator+(const FCmatrix&);// add 2 matrices of any kind
  FCmatrixFull operator-(const FCmatrix&);// sub 2 matrices of any kind
  FCmatrixFull operator*(const FCmatrix&); // mul 2 matrices of any kind
  FCmatrixFull operator*(double lambda); // mul matrix of any kind by scalar
  FCmatrixFull operator*(const Vec&); // mul matrix by Vec


  // virtual inherited
  int n_rows()const; // return number of rows
  int n_cols()const; // return number of columns
  Vec GetRow(int i) const; // retrieve row i (equivalent to [i])
  Vec GetCol(int i) const; // retrieve column i
  int GetRowMax(int i); //return index of maximum of row i
  int GetColMax(int i); //return index of relative maximum of column i
  int GetColMaxAbs(int j); //return the index of absolute maximum of column j
  double Determinant(); //return determinant of this matrix
  void Print(); //print this matrix
  void swapRows(int i,int j); //swap rows of indices i and j

private:
  int *rowindices; //row indices (0,1,2,3...)
};
#endif
