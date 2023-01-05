#include <iostream>
#include "Vec.h"
#include <vector>
#include <string>
#include <cmath>

using namespace std;


Vec::Vec()
{
  entries = NULL;
  N = 0;
}

Vec::Vec(int n) :N(n)
{
  entries = new double[N];
  for(int i=0;i<N;i++)
  {
    entries[i]=0;
  }
}

Vec::Vec(int n, double value)
{
  N = n;
  entries = new double[N];
  for(int i=0;i<N;i++)
  {
    entries[i]=value;
  }
}

Vec::Vec(int n, double * values)
{
  N = n;
  entries = new double[N];
  for(int i=0;i<N;i++)
  {
    entries[i]=values[i];
  }
}

Vec::Vec(const Vec & v)
{
  N = v.size();
  entries = new double[N];
  for(int i=0;i<N;i++)
  {
    entries[i]=v.At(i);
  }
}

Vec::~Vec()
{
  delete [] entries;
}

int Vec::size() const{
  return N;
}

double & Vec::operator[](int i)
{
  return entries[i];
}

double Vec::At(int i) const{
  return entries[i];
}

void Vec::print(){
  for(int i=0; i<N; i++)
    cout << entries[i] << "   ";
  cout << endl;
}

void Vec::SetEntries (int n, double* a){
  delete [] entries;//pos ac
  entries = new double[n];//pos ac
  N=n;//pos ac
  //if (n<=N)
  //{
    for(int i=0; i<n; i++)
    {
      entries[i]=a[i];
    }
  //}
  /*else
  {
    for(int i=0; i<N; i++)
    {
      entries[i]=a[i];
    }
    cout << "W\ARNING: Array size is bigger than Vector size"<<endl;
  }*///pos ac 
}

Vec& Vec::Add(double value)
{
  N++;
  double * temp = new double[N];
  for(int i=0; i<N-1; i++)
  {
    temp[i] = entries[i];
  }
  temp[N-1] = value;
  delete [] entries;
  entries = new double[N];
  for(int i=0; i<N; i++)
  {
    entries[i] = temp[i];
  }
  delete [] temp;
  return *this;
}

Vec& Vec::Del(int n)
{
  if(n<0 || n>=N)
  {
    cout << "\nIndex out of range! Couldn't delete!\n";
    return *this;
  }
  double * temp = new double[N-1];
  for(int i=0; i<n; i++)
  {
    temp[i] = entries[i];
  }
  for(int i=n+1; i<N; i++)
  {
    temp[i-1] = entries[i];
  }
  N = N-1;
  delete [] entries;
  entries = new double[N];
  for(int i=0; i<N; i++)
  {
    entries[i] = temp[i];
  }
  delete [] temp;
  return *this;
}

Vec& Vec::Del(int n1, int n2)
{
  if(n1<0 || n1>=N || n2<0 || n2>=N)
  {
    cout << "\nIndexes out of range! Couldn't delete!\n";
    return *this;
  }
  int tN = n2-n1+1;
  double * temp = new double[N-tN];
  for(int i=0; i<n1; i++)
  {
    temp[i] = entries[i];
  }
  for(int i=n2+1; i<N; i++)
  {
    temp[i-tN] = entries[i];
  }
  N = N-tN;
  delete [] entries;
  entries = new double[N];
  for(int i=0; i<N; i++)
  {
    entries[i] = temp[i];
  }
  delete [] temp;
  return *this;
}

Vec& Vec::operator = (const Vec & v)
{
  if(this==&v)
  {
    return *this;
  }
  N = (int) v.size();
  if(entries!=NULL)
  {
    delete [] entries;
  }
  entries = new double[N];
  for(int i=0; i<N; i++)
  {
    entries[i] = v.At(i);
  }
  return *this;
}

Vec& Vec::operator += (const Vec & v)
{
  if(N != (int) v.size())
    cout << "Vectors do not have the same dimension. Impossible operation!\n";
  else
  {
    for(int i=0; i<N; i++)
    {
      entries[i] += v.At(i);
    }
  }
  return *this;
}

Vec& Vec::operator -= (const Vec & v)
{
  if(N != (int) v.size())
    cout << "Vectors do not have the same dimension. Impossible operation!\n";
  else
  {
    for(int i=0; i<N; i++)
    {
      entries[i] -= v.At(i);
    }
  }
  return *this;
}

Vec Vec::operator + (const Vec & v)
{
  Vec temp(N);
  if(N != (int) v.size())
    cout << "Vectors do not have the same dimension. Impossible operation!\n";
  else
  {
    for(int i=0; i<N; i++)
    {
      temp[i] = entries[i] + v.At(i);
    }
  }
  return temp;
}

Vec Vec::operator - (const Vec & v)
{
  Vec temp(N);
  if(N != (int) v.size())
    cout << "Vectors do not have the same dimension. Impossible operation!\n";
  else
  {
    for(int i=0; i<N; i++)
    {
      temp[i] = entries[i] - v.At(i);
    }
  }
  return temp;
}

Vec Vec::operator - ()
{
  Vec temp(N);
  for(int i=0; i<N; i++)
  {
    temp[i] = - entries[i];
  }
  return temp;
}

Vec Vec::operator * (const Vec & v)
{
  Vec temp(N);
  if(N != (int) v.size())
    cout << "Vectors do not have the same dimension. Impossible operation!\n";
  else
  {
    for(int i=0; i<N; i++)
    {
      temp[i] = entries[i] * v.At(i);
    }
  }
  return temp;
}

Vec Vec::operator * (double value)
{
  Vec temp(N);
  for(int i=0; i<N; i++)
  {
      temp[i] = entries[i] * value;
  }
  return temp;
}

double Vec::dot(const Vec & v)
{
  double Sum = 0;
  for(int i=0; i<N; i++)
  {
      Sum += entries[i] * v.At(i);
  }
  return Sum;
}

void Vec::swap(int i1, int i2)
{
  if(0<=i1 && i1<N && 0<=i2 && i2<N)
  {
    double temp = entries [i1];
    entries[i1] = entries[i2];
    entries[i2] = temp;
  }
  else
  {
    cout << "Impossible operation! Indexes out of bounds.\n";
  }
}

double Vec::norm()
{
  double sum=0;
  for(int i=0; i <N; i++)
  {
    sum+=pow(entries[i],2);
  }
  sum = pow(sum,0.5);
  return sum;
}

double Vec::norm(int n)
{
  double sum=0;
  if(n>0)
  {
    for(int i=0; i <N; i++)
    {
      sum+=pow(entries[i],n);
    }
    sum = pow(sum,(double) 1/n);
  }
  else
  {
    sum = entries[0];
    for(int i=1; i <N; i++)
    {
      if(entries[i] > entries[i-1])
      {
        sum = entries[i];
      }
    }
  }

  return sum;
}
