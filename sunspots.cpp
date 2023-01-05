#include "EMD.h"
#include "IMF.h"
#include <vector>
#include <iostream>
#include "FCtools.h"
#include "Vec.h"
#include "cFCgraphics.h"
#include "TAxis.h"
#include <ctime>

using namespace std;

int main()
{
	//1st step: importing the data.

	vector<vector<double> > data = FCtools::VecReadFile("sunspots.dat");
	vector<double> data_time;
	vector<double> data_sunspots;
	for(int i=0; i<(int)data.size();i++)
	{
		data_time.push_back(data[i][0]);
		data_sunspots.push_back(data[i][1]);
		//cout << data_time[i] << "	" << data_sunspots[i] << endl;
	}



	//2nd step: plotting the data:
	IMF s (data_time,data_sunspots,100000);
	/*
	TGraph* signal_graph = s.GetGraph();
	signal_graph->SetLineColor(kRed+2);
	signal_graph->SetMarkerColor(kRed+2);
	signal_graph->GetXaxis()->SetLimits(1749,2018);
	signal_graph->SetTitle("Manchas Solares");


	TF1* max_graph = s.GetMAX();
	max_graph->SetLineColor(kRed-2);
	max_graph->SetMarkerColor(kRed-2);
	TF1* min_graph = s.GetMIN();
	min_graph->SetLineColor(kBlue-2);
	min_graph->SetMarkerColor(kBlue-2);

	cFCgraphics C;
	TPad *pad = C.CreatePad("pad");
	C.AddObject(signal_graph,"pad");
	//C.AddObject(s.GetMAX(),"pad","same");
	//C.AddObject(s.GetMIN(),"pad","same");
	C.AddObject(pad);
	C.Draw();*/

	//3rd step: EMD analysis
	EMD EMDdata(s);

	vector<IMF> IMFs = EMDdata.GetIMFs(); 	//vector containing all the IMFs
	EMDdata.DrawCF();
	EMDdata.DrawCorrelations();
	EMDdata.GetAnalysis();
	EMDdata.GetAnalysis(0);
	EMDdata.GetAnalysis(1);
	EMDdata.DrawHS(5);
	EMDdata.DrawHS(6);
	EMDdata.DrawHS(7);
}
