#include "Vec.h"
#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include "EqSolver.h"
#include "FCmatrixBanded.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

EqSolver::EqSolver(){}

EqSolver::EqSolver(const FCmatrix& A, const Vec& B)
{
  b = B;

  if(A.GetClassName()=="FCmatrixFull")
  {
    M=(FCmatrixFull *) &A;
  }
  if(A.GetClassName()=="FCmatrixBanded")
  {
    M =(FCmatrixBanded*) &A;
  }
}

void EqSolver::SetConstants(const Vec& B)
{
  b = B;
}

void EqSolver::SetMatrix(const FCmatrix& A)
{
  if(A.GetClassName()=="FCmatrixFull")
  {
    M = (FCmatrixFull *) &A;
  }
  if(A.GetClassName()=="FCmatrixBanded")
  {
    M= (FCmatrixBanded*) &A;
  }
}

void EqSolver::Print()
{
  cout<<"\nMatrix :"<<endl;
  M->Print();
  cout<<"Vec : ";
  b.print();
  //cout<<endl;
}

void EqSolver::GaussElimination(FCmatrix& matrix, Vec& b) //acho que devíamos ter uma função triângulação que fizesse a primeira parte por nós, porque isto também é usado no determinante, é so chato...
{
  int n_rows = matrix.n_rows();
  int n_cols = matrix.n_cols();

  FCmatrixFull A = * (FCmatrixFull*) &matrix;

  for(int i=0; i<n_cols;i++) //loop running through columns
  {
    //Get column relative maximum starting at row i
    int maxRow = A.GetColMax(i);

    //swap maximum row with current row (row switches are taken into acount in the constant Vec)
    A.swapRows(i, maxRow);
    b.swap(i, maxRow);

    //make all rows below this one 0 (considering possible changes made to constant Vec)
    for(int k=i+1; k<n_rows;k++)
    {
      double c = -A[k][i]/A[i][i];
      b[k] += c*b[i];
      for (int j=i; j<n_rows; j++)
      {
        if (i==j)
        A[k][j] = 0;
        else
        {
          A[k][j] += c * A[i][j];
        }
      }
    }
  }

  b[n_rows-1] = b[n_rows-1]/(A[n_rows-1][n_rows-1]);
  for(int i=n_rows-2;i>=0;i--)
  {
    double s=0;
    for(int c=i+1; c<n_rows;c++) //ciclo sobre as colunas da matriz
    s+=A[i][c]*b[c];
    double fac = (1/A[i][i]);
    b[i]=fac*(b[i]-s);
  }
  matrix = A;
}

Vec EqSolver::GaussEliminationSolver()
{
  Vec output(b);
  if(M->GetClassName() == "FCmatrixFull")
  {
    FCmatrixFull matrix = *(FCmatrixFull*) M ;
    int n_rows = matrix.n_rows();

    if(matrix.Determinant()==0)
    {
      cout<<"ERROR: Not possible to solve for solution vector(matrix is singular)!"<<endl;
      Vec Error;
      return Error;
    }
    else if(n_rows!=b.size())
    {
      cout<<"ERROR: Not possible to solve for solution vector (different number of equations and variables)!"<<endl;
      Vec Error;
      return Error;
    }
    else
    {
      GaussElimination(matrix, output);
    }
  }
  else if(M->GetClassName() == "FCmatrixBanded")
  {
    FCmatrixBanded temp = *(FCmatrixBanded*) M ;
    FCmatrixFull matrix = temp.toFull();
    int n_rows = matrix.n_rows();

    if(matrix.Determinant()==0)
    {
      cout<<"ERROR: Not possible to solve for solution vector(matrix is singular)!"<<endl;
      Vec Error;
      return Error;
    }
    else if(n_rows!=b.size())
    {
      cout<<"ERROR: Not possible to solve for solution vector (different number of equations and variables)!"<<endl;
      Vec Error;
      return Error;
    }
    else
    {
      GaussElimination(matrix, output);
    }
  }
  return output;
}

void EqSolver::LUdecomposition(FCmatrix& matrix, vector<int>& index)
{
  FCmatrixFull A = * (FCmatrixFull*) &matrix;
  FCmatrixFull temp(A);

  int n = A.n_rows();
  vector<Vec> v;
  for(int i=0;i<n;i++)
  {
    v.push_back(Vec(n));
  }
  FCmatrixFull L(v);

  for(int i=0;i<n;i++)
  index.push_back(i);
  for(int i=0; i<n-1;i++) //loop running through columns
  {
    if(i!=0)
    for(int j=0; j<n;j++)
    temp[i-1][j]=0;

    //Get column relative maximum starting at row i
    int maxRow = temp.GetColMax(i);

    //swap maximum row with current row
    temp.swapRows(i, maxRow);
    int t;
    t=index[i];
    index[i] = index[maxRow];
    index[maxRow] = t;


    //make all rows below this one 0
    for(int k=i+1; k<n;k++)
    {
      double c = -temp[k][i]/temp[i][i];
      for (int j=i; j<n; j++)
      {
        if (i==j) temp[k][j] = 0;
        else
        {
          temp[k][j] += c * temp[i][j];
        }
      }
    }
  }
  //trocas de linhas em A, de modo a ficar pivoted
  FCmatrixFull B(A);
  for(int i=0; i<n;i++)
  {
    if((index[i])!=i)
    {
      B[i] = A[index[i]];
    }
  }
  A=B;

  //diagonalização de A
  for(int i=0; i<n-1;i++) //loop running through columns
  {
    for(int k=i+1; k<n;k++)
    {
      double c = -A[k][i]/A[i][i];
      L[k][i] = -c;
      for (int j=i; j<n; j++)
      {
        if (i==j) A[k][j] = 0;
        else
        {
          A[k][j] += c * A[i][j];
        }
      }
    }
  }

  FCmatrixFull output (A+L);
  matrix = output;
}



Vec EqSolver::LUdecompositionSolver()
{
  if(M->GetClassName() == "FCmatrixFull" || M->GetClassName() == "FCmatrixBanded")
  {

    if(M->Determinant()==0)
    {
      cout<<"ERROR: Not possible to solve for solution vector(matrix is singular)!"<<endl;
      Vec Error;
      return Error;
    }
    else if(M->n_rows()!=b.size())
    {
      cout<<"ERROR: Not possible to solve for solution vector (different number of equations and variables)!"<<endl;
      Vec Error;
      return Error;
    }
    else
    {
      FCmatrixFull matrix;
      if(M->GetClassName() == "FCmatrixFull")
      {
        matrix = *(FCmatrixFull*) M ;
      }
      else if(M->GetClassName() == "FCmatrixBanded")
      {
        FCmatrixBanded temp = *(FCmatrixBanded*) M ;
        matrix = temp.toFull();
      }
      Vec y(b.size());
      Vec x(b.size());
      vector<int> index;
      int n = M->n_rows();
      for(int i=0; i<matrix.n_rows();i++)
      index.push_back(i);

      LUdecomposition(matrix,index);
      for(int i=0; i<n;i++)
      {
        double Sum = 0.;
        for(int j=0; j<i;j++)
        Sum+= matrix[i][j] * y[j];
        y[i] = b[index[i]] - Sum;

      }
      for(int i=n-1; i>=0 ; i--)
      {
        double Sum = 0.;
        for(int j=n-1; j>i;j--)
        {
          Sum += matrix[i][j]*x[j];
        }
        x[i] = (y[i]-Sum)/matrix[i][i];
      }
      return x;
    }
  }
  cout << "\nSomething went wrong!\n";
  Vec Error;
  return Error;
}


Vec EqSolver::ThomasAlgorithm()
{
  if(M->GetClassName()=="FCmatrixBanded")
  {
    FCmatrixBanded Mx = *(FCmatrixBanded*) M;
    Vec output(Mx.n_rows());
    if(Mx.Determinant()==0)
    {
      cout<<"determinant : "<<Mx.Determinant()<<endl;
      cout<<"The banded matrix which you are trying to use Thomas Algorithm on is singular."<<endl;
      return output;
    }
    else if(Mx.n_rows()!=b.size())
    {
      cout<<"ERROR: Not possible to solve for solution vector (different number of equations and variables)!"<<endl;
      Vec Error;
      return Error;
    }
    else
    {
    	//cout<<"Det!=0"<<endl;
      int n = Mx.n_rows();
      Vec d (n);
      Vec c (n);

      //first sweep(left to right) -> elimination of the bottom diagonal
      for (int i=0; i<n;i++)
      {
        if (i==0)
        {
          if(Mx[1][0]!=0)
          {
            c[0]=Mx[0][0]/Mx[1][0];
            d[0]=b[0]/Mx[1][0];
          }
          else
          {
            cout<<"\nPivot with 0. Executing Gauss Elimination instead.\n WARNING: This method may take more time to solve!\n";
            output = this->GaussEliminationSolver();
            return output;
          }
        }
        else if (i==n-1)
        {
          c[i] = 0;
          d[i] = (b[i]-Mx[2][i-1]*d[i-1])/(Mx[1][i]-Mx[2][i-1]*c[i-1]);
        }
        else
        {
          if(Mx[1][i] - Mx[2][i-1]*c[i-1]==0)
          {
            cout<<"\nPivot with 0. Executing Gauss Elimination instead.\n WARNING: This method may take more time to solve!\n";
            output = this->GaussEliminationSolver();
            return output;
          }
          c[i] = Mx[0][i]/(Mx[1][i] - Mx[2][i-1]*c[i-1]);
          d[i] = (b[i]-Mx[2][i-1]*d[i-1])/(Mx[1][i]-Mx[2][i-1]*c[i-1]);
        }
      }

      //right to left sweep->back substitution
      for (int i = n-1; i >= 0 ; i--)
      {
        if(i == n-1)
        {
          output[n-1] = d[n-1];
        }
        else
        output[i] = d[i] - c[i]*output[i+1];
      }
    }
    return output;
  }

  else
  {
    Vec output;
    cout<<"You are trying to call the Thomas Algorithm with a non-Banded matrix"<<endl;
    return output;
  }
}

Vec EqSolver::JacobiIterator(Vec x_init, int max_it, double tol)
{
  FCmatrixFull matrix;
  if(M->GetClassName() == "FCmatrixFull")
  {
    matrix = *(FCmatrixFull*) M ;
  }
  else if(M->GetClassName() == "FCmatrixBanded")
  {
    FCmatrixBanded temp = *(FCmatrixBanded*) M ;
    matrix = temp.toFull();
  }
  int m = x_init.size();

  // linear system of m unknowns
  Vec x(x_init);
  Vec x_aux(m); //zero’s
  bool btol = false;
  int it = 0.;
  double eps = tol;
  while (!btol && (it++ < max_it)) {
    x_aux = x;
    for (int i=0; i<m; i++) {
      x[i] = 0.;
      for (int j=0; j<m; j++)
      if (i != j) x[i] += -matrix[i][j]*x_aux[j];
      x[i] += b[i];
      x[i] /= matrix[i][i];
      //guarantee that all vector entries are converging equally
      if (fabs(x[i]-x_aux[i]) < eps) btol = true;
      else btol = false;
    }
  }

  return x;
}

Vec EqSolver::GaussSeidelIterator(Vec x_init, int max_it, double tol)
{
  FCmatrixFull matrix;
  if(M->GetClassName() == "FCmatrixFull")
  {
    matrix = *(FCmatrixFull*) M ;
  }
  else if(M->GetClassName() == "FCmatrixBanded")
  {
    FCmatrixBanded temp = *(FCmatrixBanded*) M ;
    matrix = temp.toFull();
  }
  int m = x_init.size();

  // linear system of m unknowns
  Vec x(x_init);
  Vec x_aux(m); //zero’s
  bool btol = false;
  int it = 0.;
  double eps = tol;
  while (!btol && (it++ < max_it)) {
    x_aux = x;
    for (int i=0; i<m; i++) {
      x[i] = 0.;
      for (int j=0; j<m; j++)
      if (i != j) x[i] += -matrix[i][j]*x[j];
      x[i] += b[i];
      x[i] /= matrix[i][i];
      //guarantee that all vector entries are converging equally
      if (fabs(x[i]-x_aux[i]) < eps) btol = true;
      else btol = false;
    }
  }

  return x;
}
