#ifndef __Spline3Interpolator__
#define __Spline3Interpolator__

#include "TF1.h"
#include "DataPoints.h"
#include "Vec.h"
#include <vector>
using namespace std;

class Spline3Interpolator : public DataPoints {
public:
	Spline3Interpolator();
	Spline3Interpolator(int N=0, double* x=NULL, double* y=NULL, TF1* fF0=NULL);
	Spline3Interpolator(vector<double> x, vector<double> y, TF1* fF0=NULL);
	~Spline3Interpolator() {delete FInterpolator_f; delete K;}
	double Interpolate(double x);			//assumed to be return interpolated function value at x
	TF1* GetFInterpolator() {return FInterpolator_f;}
	void Draw(); //draw points and interpolating function
	void SetFunction (TF1*);
	void Print(string FILE = ""); // print results
	void SetResolution(int n) {FInterpolator_f->SetNpx(n);}
private:
	void SetCurvatureLines();
	double FInterpolator(double*fx, double*par) {return Interpolate(fx[0]);}
	TF1* FInterpolator_f;
	TF1* F0; //eventual underlying function
	double * K; //2nd derivatives
};

#endif
