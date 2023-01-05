#include <TMath.h>
#include <stdio.h>
#include "IMF.h"
#include "EMD.h"
#include "cFCgraphics.h"
#include "TF1.h"
#include "TH1F.h"
#include "Spline3Interpolator.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "FCtools.h"
using namespace std;


int main()
{
	//Introduction
	cout<<"\n\nThis program is an implementation of the Empiric Mode Decomposition (EMD) method, used to separate a signal into Intrinsic Mode Functions-IMF (press any key to continue)."<<endl;
	getchar();

	system("clear");
	//																		//1. SINES
	cout<<"\n\nEXERCISE 1:\nFirst, we will test the method on a set of data resulting of calculating an explicitly given function (linear combination of sines) in 400 equally spaced points, ranging from 0 to 4 (press any key to see the function). \n"<<endl;
	getchar();


	//constructing function
	TF1* s = new TF1 ("f", "[0]*sin(2*TMath::Pi()*[1]*x)+[2]*sin(2*TMath::Pi()*[3]*x)+[4]*sin(2*TMath::Pi()*[5]*x)", 0, 4);
	s->SetParameter(0,0.5);
	s->SetParameter(1,1);
	s->SetParameter(2,0.8);
	s->SetParameter(3,3);
	s->SetParameter(4,1.5);
	s->SetParameter(5,5);
	s->SetNpx(400);

  	//Constructing Signal
	IMF s_sines(s,0,4,400,10000);
	//constructing EMD
	EMD EMD_sines(s_sines);



	//Drawing signal
	EMD_sines.DrawIMF(0,"sine_signal.pdf");

	cout<<"\n The method yielded these results (press any key to continue):";
	getchar();
	EMD_sines.GetAnalysis(-1,"sine_analysis_default.pdf");
	EMD_sines.GetAnalysis(1,"sine_analysis_1.pdf");
	cout<<"The correlation values and frequencies for each IMF are (press any key to continue): \n"<<endl;

	EMD_sines.DrawCF("sine_CF.pdf");

	cout<<"\n\n(Press anybutton to continue)"<<endl;

	getchar();
	system("clear");





	//																		//2. SUNSPOTS

	cout<<"\nEXERCISE 2:\nNow the method will be used on a set of data representing the number of sunspots observed since 1749. This is its plot (might take a while):\n"<<endl;

	vector<vector<double> > data = FCtools::VecReadFile("sunspots.dat");
	vector<double> data_time;
	vector<double> data_sunspots;
	for(int i=0; i<(int)data.size();i++)
	{
		data_time.push_back(data[i][0]);
		data_sunspots.push_back(data[i][1]);
	}
	IMF s_sunspots (data_time,data_sunspots,30000);
	EMD EMD_sunspots(s_sunspots);
	EMD_sunspots.DrawIMF(0,"sunspots_signal.pdf");

	cout<<"\nSince there are alot more IMFs in this problem then there were in the previous one, we will start by presenting all of them, as well as the final residue and IMF Correlation factors (press any key to continue):";
	getchar();

	EMD_sunspots.GetAnalysis(0,"sunspots_analysis_0.pdf");

	getchar();
	cout<<"Now, we present all of the IMFs and their corresponding Hilbert Spectrum. After that, we show the Correlation factor for each IMF and their average frequencies."<<endl;
	EMD_sunspots.GetAnalysis(1,"sunspots_analysis_1.pdf");
	cout<<"The correlation values are: \n"<<endl;
	vector <double> corr = EMD_sunspots.GetCorrelation();

	EMD_sunspots.DrawCF("sunspots_CF.pdf");

	cout<<"\n Now that we have the Correlation Factors, these are the IMFs and corresponding Hilbert Spectra of the most significant IMFs (namely IMFs 5,6 and 9):"<<endl;
	getchar();

	EMD_sunspots.DrawIMF(5,"sunspots_IMF_5.pdf");
	EMD_sunspots.DrawHS(5,"sunspots_IMF_5_HS.pdf");
	EMD_sunspots.DrawIMF(6,"sunspots_IMF_6.pdf");
	EMD_sunspots.DrawHS(6,"sunspots_IMF_6_HS.pdf");
	EMD_sunspots.DrawIMF(9,"sunspots_IMF_9.pdf");
	EMD_sunspots.DrawHS(9,"sunspots_IMF_9_HS.pdf");

	cout << "At last, we present all the Hilbert Spectra in one plot."<<endl;

	EMD_sunspots.DrawHS(-1,"sunspots_HS.pdf");


}
