#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include "Vec.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//construtores

FCmatrixFull::FCmatrixFull() : FCmatrix()
{
  classname ="FCmatrixFull";
  rowindices = NULL;
}
FCmatrixFull::FCmatrixFull(double **fM, int fm, int fn) : FCmatrix(fM, fm, fn)
{
  classname = "FCmatrixFull";
  rowindices = new int[fm];
  for(int i=0; i<fm; i++)
  rowindices[i]=i;
}
FCmatrixFull::FCmatrixFull (double * fM, int fm, int fn) : FCmatrix(fM, fm, fn)
{
  classname = "FCmatrixFull";
  rowindices = new int[fm];
  for(int i=0; i<fm; i++)
  rowindices[i]=i;
}

FCmatrixFull::FCmatrixFull (vector<Vec> fM) : FCmatrix(fM)
{
  classname = "FCmatrixFull";
  rowindices = new int[(int)fM.size()];
  for(int i=0; i<(int)fM.size(); i++)
    rowindices[i]=i;
}

FCmatrixFull::FCmatrixFull(const FCmatrixFull & matrix)
{
  classname = "FCmatrixFull";
  int n = matrix.n_rows();
  rowindices = new int[n];
  for(int i=0; i<n; i++)
  {
    M.push_back(matrix.GetRow(i));
    rowindices[i]=i;
  }
}

FCmatrixFull::~FCmatrixFull()
{
  M.clear();
  if(rowindices !=NULL)
  {
    delete [] rowindices;
  }
}


//OPERADORES

//IGUALDADE ENTRE MATRIZES

FCmatrixFull& FCmatrixFull:: operator=(const FCmatrixFull &matrix)
{
  if(this==&matrix)
  {
    return *this;
  }
  if(matrix.GetClassName()=="FCmatrixFull")
  {
    if(rowindices != NULL)
    {
      M.clear();
      delete [] rowindices;
    }
    int n = matrix.n_rows();
    rowindices = new int[n];
    for(int i=0; i<n; i++)
    {
      rowindices[i]=i;
      M.push_back(matrix.GetRow(i));
    }
    return *this;
  }
  else
  {
    cout << "\nDifferent kind of matrices! Operation impossible!\n";
    return *this;
  }
}

//SOMA MATRIZES
FCmatrixFull FCmatrixFull:: operator +(const FCmatrix& matrix)
{

  if(((int)M.size()!=n_rows()) || ((int)M[0].size()!=matrix.n_cols()))
  {
    cout<<"\nERROR: Sum of matrices of different dimensions!"<<endl;
    FCmatrixFull error;
    return error;
  }

  else
  {
    if(matrix.GetClassName()=="FCmatrixFull")
    {
      FCmatrixFull *temp = (FCmatrixFull *) &matrix;
      vector<Vec> result;

      for(int i=0; i<(int)M.size(); i++)
      result.push_back(this->GetRow(i)+temp->GetRow(i));

      FCmatrixFull output(result);
      return output;
    }
  }
  cout<<"\nERROR: Something has gone wrong!"<<endl;
  FCmatrixFull error;
  return error;
}

//SUBTRAÇÃO MATRIZES
FCmatrixFull FCmatrixFull:: operator -(const FCmatrix& matrix)
{

  if(((int)M.size()!=n_rows()) || ((int)M[0].size()!=matrix.n_cols()))
  {
    cout<<"\nERROR: Subtraction of matrices of different dimensions!"<<endl;
    FCmatrixFull error;
    return error;
  }
  else
  {
    if(matrix.GetClassName()=="FCmatrixFull")
    {
      FCmatrixFull *temp = (FCmatrixFull *) &matrix;
      vector<Vec> result;

      for(int i=0; i<(int)M.size(); i++)
      result.push_back(this->GetRow(i)-temp->GetRow(i));

      FCmatrixFull output(result);
      return output;
    }
  }
  cout<<"\nERROR: Something has gone wrong!"<<endl;
  FCmatrixFull error;
  return error;
}

//PRODUTO POR ESCALAR
FCmatrixFull FCmatrixFull :: operator *(double lambda)
{
  FCmatrixFull temp(*this);
  for(int i=0; i<(int)M.size(); i++)
  temp[i] = temp[i] * lambda;
  return temp;
}


//PRODUTO POR VETOR
FCmatrixFull FCmatrixFull :: operator *(const Vec& vetor)
{
  if(this->n_cols()!=vetor.size())
  {
    cout<<"\nERROR: Product of matrix by incompatible vector!"<<endl;
    FCmatrixFull error;
    return error;
  }
  else
  {
    Vec entrada_output(1);
    vector<Vec> output;
    for(int i=0; i<(int)M.size(); i++)
    output.push_back(entrada_output);
    FCmatrixFull temp(output);

    for(int i=0; i<(int)M.size(); i++)
    {
      Vec valor_temp(1);
      valor_temp[0]= this->M[i].dot(vetor);
      temp[i] = valor_temp;
    }
    return temp;
  }
}

//PRODUTO POR MATRIX...........
FCmatrixFull FCmatrixFull :: operator *(const FCmatrix& matrix)
{
  if(this->n_cols() != matrix.n_rows())
  {
    cout<<"\nERROR: Product of incompatible matrices!"<<endl;
    FCmatrixFull error;
    return error;
  }
  else
  {
    if(matrix.GetClassName()=="FCmatrixFull")
    {
      FCmatrixFull *entry = (FCmatrixFull *) &matrix;
      vector<Vec> pre_output;
      for(int i=0; i<this->n_rows();i++)
      {
        double entrada;
        Vec temp(entry->n_cols(), 0.);
        for(int j=0; j<entry->n_cols();j++)
        {
          entrada = M[i].dot(entry->GetCol(j));
          temp[j] = entrada;
        }
        pre_output.push_back(temp);
      }
      FCmatrixFull output(pre_output);
      return output;
    }
  }
  cout<<"\nERROR: Something has gone wrong!"<<endl;
  FCmatrixFull error;
  return error;
}



//Methods

double FCmatrixFull::Determinant()
{
  //GaussElimination
  FCmatrixFull A(*this);
  int n_rows = A.n_rows();
  int n_cols = A.n_cols();

  if(n_rows!=n_cols)
  {
    return 0;
  }
  else
  {
    int swapfactor = 1;
    for(int i=0; i<n_cols-1;i++) //loop running through columns
    {
      //Get column relative maximum starting at row i
      int maxRow = A.GetColMax(i);

      //swap maximum row with current row
      A.swapRows(i, maxRow);

      //Updates swapfactor if rows were swapped
      if(i!=maxRow)
      {
        swapfactor = swapfactor * (-1);
      }

      //make all rows below this one 0
      for(int k=i+1; k<n_rows;k++)
      {
        double c = -A[k][i]/A[i][i];
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

    //Calculate Determinant by Multiplying pivots
    double det = 1;
    for(int i=0; i<n_rows; i++)
    {
      det= det * A[i][i];
    }
    det = det * swapfactor;
    return det;
  }
}


int FCmatrixFull::n_rows() const {return M.size();}

int FCmatrixFull::n_cols() const {return M[0].size();}

Vec& FCmatrixFull::operator[] (int i){return M[rowindices[i]];} //returns row n

Vec FCmatrixFull::GetRow(int i)const
{
  return M[rowindices[i]];
}

Vec FCmatrixFull::GetCol(int i)const
{
  Vec output(M.size());
  for(int j=0; j<(int)M.size(); j++)
  {
    output[j] = M[rowindices[j]].At(i);
  }
  return output;
}

void FCmatrixFull::swapRows(int i, int j)
{
  int temp;
  temp = rowindices[i];
  rowindices[i] = rowindices[j];
  rowindices[j] = temp;
}

int FCmatrixFull::GetRowMax(int j)
{
  Vec temp = GetRow(j);
  double valmax = fabs(temp[0]);
  int indexmax = 0;
  for(int i=1; i<temp.size(); i++)
  {
    if(valmax<fabs(temp[i]))
    {
      valmax = fabs(temp[i]);
      indexmax = i;
    }
  }
  return indexmax;
}

int FCmatrixFull::GetColMax(int j)
{
  Vec col = GetCol(j);
  double valmax = fabs(col[j]/(GetRow(j).norm(1)));
  int indexmax = j;

  for(int i=j; i<col.size(); i++)
    {
	     if(valmax<fabs(col[i]/(GetRow(i).norm(1))))
		{
		  valmax = fabs(col[i]/(GetRow(i).norm(1)));
		  indexmax = i;
		}
    }
  return indexmax;
}

int FCmatrixFull::GetColMaxAbs(int j)
{
  Vec col = GetCol(j);
  double valmax = fabs(col[0]);
  int indexmax = 0;

  for(int i=j; i<col.size(); i++)
    {
	     if(valmax<fabs(col[i]))
		{
		  valmax = fabs(col[i]);
		  indexmax = i;
		}
    }
  return indexmax;
}

void FCmatrixFull:: Print(){
  for(int i=0; i<this->n_rows(); i++)
    {
      this->GetRow(i).print();
    }
  //cout<<endl;
}
