#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "Spline3Interpolator.h"
#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include "FCmatrixBanded.h"
#include "EqSolver.h"
#include <TF1.h>
#include <TGraph.h>
#include "cFCgraphics.h"
#include <vector>
using namespace std;

using namespace std;

Spline3Interpolator::Spline3Interpolator()
{;}

Spline3Interpolator::Spline3Interpolator(int fN, double *fx, double *fy, TF1* fF0) : DataPoints(fN,fx,fy)
{
	SetCurvatureLines();
	F0 = fF0;
	FInterpolator_f = new TF1("FInterpolator", this, &Spline3Interpolator::FInterpolator, x[0] ,x[N-1], 0);
	int flag = 0;
	for(int i=1;i<fN;i++)
	{
		if(fy[i]!=fy[i-1]) flag=1;
	}
	if(flag==0)
	{
		FInterpolator_f = new TF1("FInterpolator", "[0]", x[0] ,x[N-1]);
		FInterpolator_f->SetParameter(0,fy[0]);
	}
}

Spline3Interpolator::Spline3Interpolator(vector<double> fx, vector<double> fy, TF1* fF0) : DataPoints(fx,fy)
{
	SetCurvatureLines();
	F0 = fF0;
	FInterpolator_f = new TF1("FInterpolator", this, &Spline3Interpolator::FInterpolator, x[0] ,x[N-1], 0);
	int flag = 0;
	for(int i=1;i<(int)fx.size();i++)
	{
		if(fy[i]!=fy[i-1]) flag=1;
	}
	if(flag==0)
	{
		FInterpolator_f = new TF1("FInterpolator", "[0]", x[0] ,x[N-1]);
		FInterpolator_f->SetParameter(0,fy[0]);
	}
}

void Spline3Interpolator::SetCurvatureLines()
{
	K = new double[N];
	K[0]=0.;
	K[N-1]=0.;
	Vec b(N-2); //constants

	vector <Vec> banded; //trigonal matrix problem

	Vec upper(N-3);
	Vec lower(N-3);
	Vec diag(N-2);
	for(int i=0;i<N-2;i++)
	{
		b[i] = 6*((y[i]-y[i+1])/(x[i]-x[i+1]) -(y[i+1]-y[i+2])/(x[i+1]-x[i+2]));
		if(i==0)
		{
			diag[i]=2*(x[i]-x[i+2]);
			upper[i] = x[i+1]-x[i+2];
		}
		else if(i==N-3)
		{
			diag[i]=2*(x[i]-x[i+2]);
			lower[i-1]=x[i]-x[i+1];
		}
		else
		{
			diag[i]=2*(x[i]-x[i+2]);
			lower[i-1]=x[i]-x[i+1];
			upper[i] = x[i+1]-x[i+2];
		}
	}
	banded.push_back(upper);
	banded.push_back(diag);
	banded.push_back(lower);

	FCmatrixBanded mx_banded (banded);

	EqSolver system_banded(mx_banded, b);
	//system_banded.Print();

	//cout<<"The vector containing the desired curvatures is the solution to the above system:";
	Vec curv = system_banded.ThomasAlgorithm();
	//curv.print();
	//getchar();
//	double temp[N];
	//K.print();
	//temp[0]=0;
	for(int i=1; i<N-1; i++)
	{
		K[i] = curv[i-1];
	}
	//temp[N-1]=0;
	//K.SetEntries(N,temp);
	//K.print();
	//getchar();
}

double Spline3Interpolator::Interpolate(double xval)
{
	double f,A,B,C;

	// detect in which segment is x
	for (int i=0; i<N; i++)
	{
		if ((xval-x[N-1])>0.) // out of range
		{
			return 0.;
		}
		if ((xval-x[i])<0.)
		{
			//upper bound returned
			if (i==0) // out of range
			{
				return 0.;
			}
			else
			{
				A = K[i-1]/6*(pow(xval-x[i],3)/(x[i-1]-x[i])-(xval-x[i])*(x[i-1]-x[i]));
				B = K[i]/6*(pow(xval-x[i-1],3)/(x[i-1]-x[i])-(xval-x[i-1])*(x[i-1]-x[i]));
				C = (y[i-1]*(xval-x[i])-y[i]*(xval-x[i-1]))/(x[i-1]-x[i]);
				f = A - B + C;
				//printf("A:%lf	B:%lf	C:%lf	f:%lf\n",A,B,C,f);
				break;
			}
		}
	}

	//retrieve segment interpolator and return function value


	return f;

}

void Spline3Interpolator::Draw()
{
	cFCgraphics G;
	TGraph *g = new TGraph(N,x,y);
	g->SetMarkerStyle(20);
	g->SetMarkerColor(kRed-3);
	g->SetMarkerSize(1.5);

	TPad *pad1 = G.CreatePad("pad1");

	FInterpolator_f->SetLineColor(38);
	FInterpolator_f->SetLineWidth(4);
	FInterpolator_f->SetTitle("Spline3");
	G.AddObject(FInterpolator_f,"pad1");

	if (F0)
	{
		F0->SetNpx(FInterpolator_f->GetNpx());
		F0->SetLineColor(kGray+3);
		F0->SetLineWidth(4);
		F0->SetLineStyle(5);
		G.AddObject(F0,"pad1", "same");
	}

	G.AddObject(g,"pad1","P");
	G.DumpPad("pad1");
	G.AddObject(pad1);
	G.Draw();


	delete g;
}

void Spline3Interpolator::SetFunction(TF1* f)
{
	F0 = f;
}

void Spline3Interpolator::Print(string FILE)
{
	cFCgraphics G;
	TGraph *g = new TGraph(N,x,y);
	g->SetMarkerStyle(20);
	g->SetMarkerColor(kRed-3);
	g->SetMarkerSize(1.5);

	TPad *pad1 = G.CreatePad("pad1");

	FInterpolator_f->SetLineColor(38);
	FInterpolator_f->SetLineWidth(4);
	FInterpolator_f->SetTitle("Spline3");
	G.AddObject(FInterpolator_f,"pad1");

	if (F0)
	{
		F0->SetNpx(FInterpolator_f->GetNpx());
		F0->SetLineColor(kGray+3);
		F0->SetLineWidth(4);
		F0->SetLineStyle(5);
		G.AddObject(F0,"pad1", "same");
	}

	G.AddObject(g,"pad1","P");
	G.DumpPad("pad1");
	G.AddObject(pad1);
	G.Print(FILE);
}
