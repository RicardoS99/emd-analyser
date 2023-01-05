#ifndef _A_
#define _A_

#include "Vec.h"
#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include <vector>
using namespace std;

class EqSolver {
public:
  EqSolver();
  EqSolver(const FCmatrix&, const Vec&); // matrix M and vector of constants B, functional
  ~EqSolver() {;}
  //set
  void SetConstants(const Vec&); //changing the vector of constants of the system, functional
  void SetMatrix(const FCmatrix&); //changing the matrix of the system, functional
  //methods
  Vec GaussEliminationSolver();
  Vec LUdecompositionSolver();
  Vec ThomasAlgorithm();
  Vec JacobiIterator(Vec x_init, int max_it = 1000, double tol=1e-9);
  Vec GaussSeidelIterator(Vec x_init, int max_it = 1000, double tol=1e-9);
  void Print();

private:
  //Decomposition LU with |L|=1
  void LUdecomposition(FCmatrix&, vector<int>& index);
  //return triangular matrix and changed vector of constants
  void GaussElimination(FCmatrix&, Vec&);
  FCmatrix *M; //matrix of coefficients
  Vec b; //vector of constants
};

#endif
