#include "IMF.h"
#include <vector>
#include <iostream>
#include "Spline3Interpolator.h"
#include <TGraph.h>
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TAxis.h"
#include "cFCgraphics.h"
#include <cmath>
#include <string>
#include <ctime>

using namespace std;



IMF::IMF(double * x, double * y, int n, int R) : N(n), EN(0), Resolution(R){
	for(int i=0; i<N; ++i){
		X.push_back(x[i]);
		Y.push_back(y[i]);
	}
	if(N>0){
		SetMaximum();
		SetMinimum();
		SetMedium();
		SetAmplitude();
	}
}



IMF::IMF(vector<double> x, vector<double> y, int R) : N(x.size()), EN(0), Resolution(R){
	if(x.size()==y.size())
	{
		X=x;
		Y=y;
		SetMaximum();
		SetMinimum();
		SetMedium();
		SetAmplitude();
	}
	else cout << "VECTORS DO NOT HAVE THE SAME SIZE!" << endl;
}



IMF::IMF(TF1 * F, double xmin, double xmax, int n, int R) : N(n), EN(0), Resolution(R){
	if(xmin>xmax)
	{
		double temp = xmin;
		xmin = xmax;
		xmax = temp;
	}
	if(N<0.) N = -N;
	double step = (xmax-xmin)/(double)N;
	for(int i=0; i<N; ++i){
		X.push_back(xmin+i*step);
		Y.push_back(F->Eval(xmin+i*step));
	}
	SetMaximum();
	SetMinimum();
	SetMedium();
	SetAmplitude();
}



IMF::IMF(const IMF & original) : N(original.GetN()), Resolution(original.Resolution){
	EN = 0;
	X=original.vectorX();
	Y=original.vectorY();
	SetMaximum();
	SetMinimum();
	SetMedium();
	SetAmplitude();
}

IMF::~IMF()
{
	EN=0;
	delete MED;
	delete AMP;
}

IMF IMF::operator+(const IMF & original){
	EN = 0;

	if(N==0)
	{
		IMF output = original;
		return output;
	}

	else if(N!=(int)original.vectorX().size())
	{
		cout<<"It was not possible to sum the two IMFs! (different size)"<<endl;
		return *this;
	}
	vector<double> temp=original.vectorY();
	vector<double> y;
	for(int i=0; i<N;i++)
	{
		y.push_back(temp[i]+Y[i]);
	}
	IMF output(X,y);

	return output;
}

IMF IMF::operator-(const IMF & original){
	EN = 0;

	if(N!=(int)original.vectorX().size())
	{
		cout<<"It was not possible to sum the two IMFs! (different size)"<<endl;
		return *this;
	}
	vector<double> temp=original.vectorY();
	vector<double> y;
	for(int i=0; i<N;i++)
	{
		y.push_back(Y[i]-temp[i]);
	}
	IMF output(X,y);

	return output;
}




IMF& IMF::operator = (const IMF & original){
	if(this!=&original)
	{
		N = original.GetN();
		EN = 0;
		X.clear();
		Y.clear();
		Resolution = original.Resolution;
		/*if(N>0){
		//delete MAX;
		//delete MIN;
		delete MED;
		delete AMP;
	}*/
	X = original.vectorX();
	Y = original.vectorY();
	SetMaximum();
	SetMinimum();
	SetMedium();
	SetAmplitude();
}
return *this;
}

void IMF::SetPoints(double * x, double * y, int N){
	EN=0;
	X.clear();
	Y.clear();
	if(N>0){
		delete MAX;
		delete MIN;
		delete MED;
		delete AMP;
	}
	for(int i=0; i<N; ++i){
		X.push_back(x[i]);
		Y.push_back(y[i]);
	}
	SetMaximum();
	SetMinimum();
	SetMedium();
	SetAmplitude();
}



void IMF::SetPoints(vector<double> x, vector<double> y){
	if(x.size()==y.size())
	{
		EN=0;
		X.clear();
		Y.clear();
		if(N>0){
			delete MAX;
			delete MIN;
			delete MED;
			delete AMP;
		}
		X=x;
		Y=y;
		SetMaximum();
		SetMinimum();
		SetMedium();
		SetAmplitude();
	}
	else cout << "VECTORS DO NOT HAVE THE SAME SIZE!" << endl;
}



double * IMF::GetX(){
	double * output = new double [N];
	for(int i=0; i<N; ++i){
		output[i] = X[i];
	}
	return output;
}



double * IMF::GetY(){
	double * output = new double [N];
	for(int i=0; i<N; ++i){
		output[i] = Y[i];
	}
	return output;
}

vector<double> IMF::vectorX() const {return X;}

vector<double> IMF::vectorY() const {return Y;}

TGraph * IMF::GetGraph(){
	TGraph * graph = new TGraph(N,this->GetX(),this->GetY());
	return graph;
}


void IMF::SetMedium(){
	MED = new TF1("Medium",this,&IMF::Medium,X.front(),X.back(),0);
	MED->SetNpx(Resolution);
}

void IMF::SetAmplitude(){
	AMP = new TF1("Amplitude",this,&IMF::Amplitude,X.front(),X.back(),0);
	AMP->SetNpx(Resolution);
}


void IMF::SetMaximum(){
	vector <double> y_max;
	vector <double> x_max;
	int k=0;
	for(int i=1; i<N-1;i++)
	{
		if(Y[i]>=Y[i-1] && Y[i]>=Y[i+1])
		{
			if((Y[i]!=Y[i-1] || Y[i]!=Y[i+1]))
			{
				k++;//
				y_max.push_back(Y[i]);
				x_max.push_back(X[i]);
			EN++;
		}

		}
	}

	if (k==0){
		MAX = new TF1("output","[0]",X[0],X[N-1]);
		if(Y.front()>Y.back()) MAX->SetParameter(0,Y.front());
		else MAX->SetParameter(0,Y.back());
		cout<<"NENHUM MAXIMO encontrado"<<endl; 		//			//
	}
	else{
		//EN += x_max.size();

		vector<vector<double> > mirror = Mirror(x_max, y_max,1);

		//generating the interpolator
		Spline3Interpolator max_spline (mirror[0], mirror[1]);
		max_spline.SetResolution(Resolution);
		MAX = (TF1*)((max_spline.GetFInterpolator())->Clone("max"));
		MAX->SetNpx(Resolution);
	}
}


void IMF::SetMinimum(){
	vector <double> y_min;
	vector <double> x_min;
	int k=0; 			//deteta se ha ou nao mínimos //

	for(int i=1; i<N-1;i++)
	{
		if(Y[i]<=Y[i-1] && Y[i]<=Y[i+1])
		{
			if((Y[i]!=Y[i-1] || Y[i]!=Y[i+1]))
			{
				k++;
				y_min.push_back(Y[i]);
				x_min.push_back(X[i]);
				EN++;
			}
		}
	}
	if (k==0){
		MIN = new TF1("output","[0]",X[0],X[N-1]);
		if(Y.front()<Y.back()) MIN->SetParameter(0,Y.front());
		else MIN->SetParameter(0,Y.back());
		cout<<"NENHUM MINIMO encontrado"<<endl; 		//
	}
	else{
		//EN += x_min.size();

		vector<vector<double> > mirror = Mirror(x_min, y_min,0);

		//generating the interpolator
		Spline3Interpolator min_spline (mirror[0], mirror[1]);
		min_spline.SetResolution(Resolution);
		MIN = (TF1*)((min_spline.GetFInterpolator())->Clone("min"));
		MIN->SetNpx(Resolution);
	}
}

vector<double> IMF::GetRoots()         //Returns vector with roots of function
{
	vector<double> output;
	for(int i=0; i<N-1; i++)
	{
		if(Y[i]*Y[i+1]<0.)
		{
			output.push_back((X[i]+X[i+1])/2);
		}
		if(Y[i]==0 && Y[i]!=0 && Y[i-1]!=0 && i!=0)
		{
			output.push_back((X[i]+X[i+1])/2);
		}
	}

	//else if(Y[i]==0 && Y[i]==0 && Y[i-1]==0)output.push_back((X[i]+X[i+1])/2);
	//}
	//if(Y[N-1]==0)output.push_back((X[N-2]+X[N-1])/2);

/*
		else if(Y[i+1]==0. && Y[i]!=0)
		{
			int j=2;
			while(1)
			{
				if(Y[i+j]==0.)
				{
					j++;
				}
				else
				{
					output.push_back((X[i+1]+X[i+j-1])/2.);
					break;
				}
			}
		}
	}*/
	return output;
}

vector<vector<double> > IMF::Mirror(vector<double>x, vector<double>y, int flag){
	vector<double> MX;
	vector<double> MY;
	int N = x.size();
	MX.push_back(2*X.front()-x[0]);
	MY.push_back(y[0]);

	if((Y.front()>=y[0]) && flag == 1)
	{
		MX.push_back(X.front());
		MY.push_back(Y.front());
	}
	if((Y.front()<=y[0]) && flag == 0)
	{
		MX.push_back(X.front());
		MY.push_back(Y.front());
	}
	for(int i=0;i<N;i++){
		MX.push_back(x[i]);
		MY.push_back(y[i]);
	}

	if((Y.back()>=y[N-1]) && flag == 1)
	{
		MX.push_back(X.back());
		MY.push_back(Y.back());
	}
	if((Y.back()<=y[N-1]) && flag == 0)
	{
		MX.push_back(X.back());
		MY.push_back(Y.back());
	}
	MX.push_back(2*X.back()-x[N-1]);
	MY.push_back(y[N-1]);
	vector<vector<double> > output;
	output.push_back(MX);
	output.push_back(MY);
	return output;
}

TGraph * IMF::GetFrequency(){
	//will first calculate the frequencies and then
	//proceed to splining over them
	vector<double> roots = GetRoots();
	vector<double> xf;
	vector<double> yf;

	for(int i=1;i<(int)roots.size();i++){
		xf.push_back((roots[i]+roots[i-1])/2.);
		yf.push_back(0.5/(roots[i]-roots[i-1]));
	}
	//cout << "Average Frequency: " << avg << endl;
	if((int)roots.size()>0 && (int)roots.size()<N/2){

		vector<vector<double> > mirror = Mirror(xf, yf,2);
		Spline3Interpolator spline(mirror[0],mirror[1]);
		if(N>2000) spline.SetResolution(N);
		else spline.SetResolution(2000);
		TF1* freqfunc = spline.GetFInterpolator();
		TGraph * output = new TGraph(freqfunc);
		return output;
	}
	else if((int)roots.size()>=N/2)
	{
		double * xd = new double[xf.size()];
		double * yd = new double[xf.size()];
		for(int i=0; i<(int)xf.size(); i++){
			xd[i]=xf[i];
			yd[i]=yf[i];
		}
		TGraph * output = new TGraph(xf.size(),xd,yd);
		delete xd;
		delete yd;
		return output;
	}
	TF1 * zero= new TF1("zero","0",X.front(),X.back());
	TGraph * output = new TGraph(zero);
	delete zero;
	return output;
}

double IMF::Medium(double*fx, double*par){
	double output = (MAX->Eval(fx[0])+MIN->Eval(fx[0]))/2.;
	return output;
}

double IMF::Amplitude(double*fx, double*par){
	double output = (MAX->Eval(fx[0])-MIN->Eval(fx[0]))/2.;
	return output;
}

int IMF::isIMF(){
	double tol=0.05;
	double weight=0.95;
	//Extrema and Zero test
	int ZN = (int)GetRoots().size();
	if(abs(ZN-EN)>1){
		return 0;
	}

	//Rilling Criteria

	int counter = 0;
	vector<double> var;
	for(int i=0; i<N; ++i){
		var.push_back(fabs(MED->Eval(X[i])/AMP->Eval(X[i])));
		if(var[i]>tol*10.){
			return 0;
		}
		if(var[i]<tol) counter ++;
	}
	if(((double)counter/(double)N)<weight) return 0;

	return 1;
}

int IMF::isLastIMF(){
	int counter = 0;
	double tol = 0.1;
	double avg = 0;
	double weight = 1;

	//3 Extrema Test
	if(EN<=3) return 1;

	//Monotone Test

	for(int i=1;i<N;i++){
		if(Y[i]>=Y[i-1]) counter ++;
	}
	if(counter==N) return 1;
	counter = 0;
	for(int i=1;i<N;i++){
		if(Y[i]<=Y[i-1]) counter ++;
	}
	if(counter==N) return 1;

	//Constant Test
	counter = 0;
	for(int i=0;i<N;i++){
		avg += Y[i];
	}
	avg = avg/N;
	for(int i=0; i<N; ++i){
		if(fabs(Y[i]-avg)<tol) counter ++;
	}
	if((double)counter/(double)N>=weight) return 1;

	return 0;
}

void IMF::GetDetail(){
	for(int i=0; i<N;i++)
	{
		Y[i]-=MED->Eval(X[i]);
	}
	delete MED;
	delete AMP;
	EN=0;
	SetMaximum();
	SetMinimum();
	SetMedium();
	SetAmplitude();
}
