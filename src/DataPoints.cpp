#include "cFCgraphics.h"
#include "DataPoints.h"
#include <string>
#include "TGraph.h"
#include "TPad.h"
#include <iostream>
#include <vector>

using namespace std;

DataPoints::DataPoints() {
  N=0;
  x=NULL;
  y=NULL;
}

DataPoints::DataPoints(int n, double*xpoints, double*ypoints) : N(n) //automatic assignment of N
{
  x = new double[N];
  y = new double[N];
  for(int i=0; i<N; i++)
  {
    x[i]=xpoints[i];
    y[i]=ypoints[i];
  }
}

DataPoints::DataPoints(vector<double> xpoints, vector<double> ypoints) : N(xpoints.size()) //automatic assignment of N
{
  x = new double[N];
  y = new double[N];
  for(int i=0; i<N; i++)
  {
    x[i]=xpoints[i];
    y[i]=ypoints[i];
  }
}

DataPoints::~DataPoints()
{
  delete [] x;
  delete [] y;
}

void DataPoints::Print()
{
  cout<<"There are "<< N <<" points stored:\nFormat(x,y)"<<endl;
  for(int i=0; i<N; i++)
    cout<<"( "<<x[i]<<" , "<<y[i]<<" )"<<endl;
}
