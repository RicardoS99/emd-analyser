#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include "FCmatrixBanded.h"
#include "Vec.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//construtores

FCmatrixBanded::FCmatrixBanded():FCmatrix()
{
	classname = "FCmatrixBanded";
}

FCmatrixBanded::FCmatrixBanded(double **fM, int fm, int fn):FCmatrix(fM, fm, fn)
{
	classname = "FCmatrixBanded";
	M[0].Del(fn-1);
	M[2].Del(fn-1);
}

FCmatrixBanded::FCmatrixBanded (double * fM, int fm, int fn):FCmatrix(fM, fm, fn)
{
	classname = "FCmatrixBanded";
	M[0].Del(fn-1);
	M[2].Del(fn-1);
}
FCmatrixBanded::FCmatrixBanded (vector<Vec> fM):FCmatrix(fM)
{
	classname = "FCmatrixBanded";
}

FCmatrixBanded::FCmatrixBanded(const FCmatrixBanded & matrix) : FCmatrix(matrix)
{
	classname = "FCmatrixBanded";
}

FCmatrixFull FCmatrixBanded::toFull()
{
	int dim = M[1].size();
	int dim2 = dim*dim;
	double * temp = new double [dim2];
	for(int i=0; i<dim2; i++)
	{
		temp[i]=0;
	}
	temp[0] = M[1][0];
	temp[1] = M[0][0];
	for(int i=0; i<dim-2; i++)
	{
		temp[i+dim*(i+1)] = M[2][i];
		temp[i+1+dim*(i+1)] = M[1][i+1];
		temp[i+2+dim*(i+1)] = M[0][i+1];
	}
	temp[dim2-2] = M[2][dim-2];
	temp[dim2-1] = M[1][dim-1];

	FCmatrixFull output(temp,dim,dim);
	return output;
}
/*
//OPERADORES

//IGUALDADE ENTRE MATRIZES

void FCmatrixBanded :: operator=(const FCmatrix &matrix)
{
  M.clear();
  int n_rows = matrix.n_rows();
  delete [] rowindices;

  if(matrix.GetClassName()=="FCmatrixBanded")
  {
    rowindices = new int[n_rows];
    for(int i=0; i<n_rows; i++)
    {
      rowindices[i] = i;
      M.push_back(matrix.GetRow(i));
    }
  }
}

//SOMA MATRIZES
FCmatrixBanded FCmatrixBanded:: operator +(const FCmatrix& matrix)
{

  if(((int)M.size()!=n_rows()) || ((int)M[0].size()!=matrix.n_cols()))
  {
    cout<<"\nERROR: Sum of matrices of different dimensions!"<<endl;
    FCmatrixBanded error;
    return error;
  }

  else
  {
    if(matrix.GetClassName()=="FCmatrixBanded")
    {
      FCmatrixBanded *temp = (FCmatrixBanded *) &matrix;
      vector<Vec> result;

      for(int i=0; i<(int)M.size(); i++)
      result.push_back(this->GetRow(i)+temp->GetRow(i));

      FCmatrixBanded output(result);
      return output;
    }
  }
  cout<<"\nERROR: Something has gone wrong!"<<endl;
  FCmatrixBanded error;
  return error;
}

//SUBTRAÇÃO MATRIZES
FCmatrixBanded FCmatrixBanded:: operator -(const FCmatrix& matrix)
{

  if(((int)M.size()!=n_rows()) || ((int)M[0].size()!=matrix.n_cols()))
  {
    cout<<"\nERROR: Subtraction of matrices of different dimensions!"<<endl;
    FCmatrixBanded error;
    return error;
  }
  else
  {
    if(matrix.GetClassName()=="FCmatrixBanded")
    {
      FCmatrixBanded *temp = (FCmatrixBanded *) &matrix;
      vector<Vec> result;

      for(int i=0; i<(int)M.size(); i++)
      result.push_back(this->GetRow(i)-temp->GetRow(i));

      FCmatrixBanded output(result);
      return output;
    }
  }
  cout<<"\nERROR: Something has gone wrong!"<<endl;
  FCmatrixBanded error;
  return error;
}

//PRODUTO POR ESCALAR
FCmatrixBanded FCmatrixBanded :: operator *(double lambda)
{
  FCmatrixBanded temp(*this);
  for(int i=0; i<(int)M.size(); i++)
  temp[i] = temp[i] * lambda;
  return temp;
}


//PRODUTO POR VETOR
FCmatrixBanded FCmatrixBanded :: operator *(const Vec& vetor)
{
  if(this->n_cols()!=vetor.size())
  {
    cout<<"\nERROR: Product of matrix by incompatible vector!"<<endl;
    FCmatrixBanded error;
    return error;
  }
  else
  {
    Vec entrada_output(1);
    vector<Vec> output;
    for(int i=0; i<(int)M.size(); i++)
    output.push_back(entrada_output);
    FCmatrixBanded temp(output);

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
FCmatrixBanded FCmatrixBanded :: operator *(const FCmatrix& matrix)
{
  if(this->n_cols() != matrix.n_rows())
  {
    cout<<"\nERROR: Product of incompatible matrices!"<<endl;
    FCmatrixBanded error;
    return error;
  }
  else
  {
    if(matrix.GetClassName()=="FCmatrixBanded")
    {
      FCmatrixBanded *entry = (FCmatrixBanded *) &matrix;
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
      FCmatrixBanded output(pre_output);
      return output;
    }
  }
  cout<<"\nERROR: Something has gone wrong!"<<endl;
  FCmatrixBanded error;
  return error;
}


*/
//Methods

double FCmatrixBanded::Determinant()
{
	//contninuant
	FCmatrixBanded A(*this);
	int n_rows = A.n_rows();
	int n_cols = A.n_cols();

	if(n_rows!=n_cols)
	{
		return 0;
	}
	else
	{
  	//the determinant of a banded matrix is a fibbonacci sequence
		vector <double> K(n_rows);
		K[0] = 1;
		K[1] = M[1][0];
		for(int i=2; i<n_rows;i++)
		{
			K[i] = M[1][i]*K[i-1] - M[0][i-1]*M[2][i-1]*K[i-2];
      //cout<<" i :"<<K[i]<<endl;
		}
    //cout<<endl;
		return K[n_rows-1];
	}
}


int FCmatrixBanded::n_rows() const {return M[1].size();}

int FCmatrixBanded::n_cols() const {return M[1].size();}
/*

Vec FCmatrixBanded::GetRow(int i)const
{
  return M[rowindices[i]];
}

Vec FCmatrixBanded::GetCol(int i)const
{
  Vec output(M.size());
  for(int j=0; j<(int)M.size(); j++)
  {
    output[j] = M[rowindices[j]].At(i);
  }
  return output;
}

void FCmatrixBanded::swapRows(int i, int j)
{
  int temp;
  temp = rowindices[i];
  rowindices[i] = rowindices[j];
  rowindices[j] = temp;
}

int FCmatrixBanded::GetRowMax(int i)
{
  Vec temp = GetRow(i);
  return temp.norm(0);
}

int FCmatrixBanded::GetColMax(int i)
{
  Vec temp = GetCol(i);
  return temp.norm(0);
}
*/
Vec FCmatrixBanded::GetDiagonal(int i)
{
	if(i==1 || i==2 || i==3)
	{
		return M[i-1];
	}
	else
	{
		cout << "There is no diagonal " << i << ".! i must be 1,2 or 3.\n";
		Vec error;
		return error;
	}
}
void FCmatrixBanded:: Print()
{
	FCmatrixFull temp = this->toFull();
	temp.Print();
}
