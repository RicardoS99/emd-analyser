#include "EMD.h"
#include <vector>
#include <iostream>
#include "Spline3Interpolator.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TBox.h"
#include "cFCgraphics.h"
#include <cmath>
#include <ctime>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
using namespace std;

//constructors

EMD::EMD(IMF fsignal) : residue(fsignal){
	IMFs.push_back(fsignal);
	IMFanalyser();
	for(int i=1; i<(int)IMFs.size();i++)
	{
		correlation.push_back(CalculateCorrelation(i));
		frequencies.push_back(CalculateFrequency(i));
	}
}

IMF EMD::GetIMF(int n){
	return IMFs[n];
}


//methods

void EMD::IMFanalyser(){
	//cout << "\n-----STARTING EMD -----"<<endl;
	int i=0;
	int j;
	//clock_t begin;
	//clock_t end;
	while(!residue.isLastIMF()) //loop that generates the vector IMFs
	{
		IMF detail = residue;
		j=0;
		//cout << endl;
		while(!detail.isIMF())	//loop that search for next IMF
		{
			//begin = clock();
			detail.GetDetail();
			//end = clock();
			j++;
			//cout << j << "	Time: " << double(end-begin)/CLOCKS_PER_SEC <<endl;
		}

		IMFs.push_back(detail); //store IMF in vector
		residue = residue - detail; //calculate residue
		i++;
		//		cout << "\nIMF CALCULATED	i: "<< i << "	Iterations: " << j << endl;
	}
	//	cout << "\n-----EMD FINISHED-----" << endl;
}

double EMD::CalculateCorrelation(int n){
	vector<double> IMFn = IMFs[n].vectorY();
	vector<double> IMF0 = IMFs[0].vectorY();
	double sum1=0;
	double sum2=0;
	double sum3=0;
	for(int i=0; i<(int)IMFn.size();i++)
	{
		sum1+= IMF0[i]*IMFn[i];
		sum2+= IMF0[i]*IMF0[i];
		sum3+= IMFn[i]*IMFn[i];
	}
	return sum1/sqrt(sum2*sum3);
}

double EMD::CalculateFrequency(int n){
	TGraph * F = IMFs[n].GetFrequency();
	int N = F->GetN();
	double * y = F->GetY();
	double avg=0;
	for(int i=0; i<N; i++){
		avg+=y[i];
	}
	avg = avg/(double)N;
	return avg;
}

void EMD::GetAnalysis(int flag, string filename){
	if(flag==0)
	{
		//Getting TGraphs and Setting Colors
		int N = IMFs.size();
		vector<TGraph *> gimf;
		for(int i=0; i<N;i++){
			gimf.push_back(IMFs[i].GetGraph());
			gimf[i]->SetLineColor(i+40);
			gimf[i]->SetMarkerColor(i+40);
			gimf[i]->GetXaxis()->SetTitle("Time");
			gimf[i]->GetYaxis()->SetTitle("Amplitude");
		}
		TGraph * gres = residue.GetGraph();
		gres->GetXaxis()->SetTitle("Time");
		gres->GetYaxis()->SetTitle("Amplitude");

		TLegend * lc = new TLegend(0.7,0.7,0.95,0.95);
		vector <double> correlation = GetCorrelation();
		double width = 0.25;
		TH1D * Frame = new TH1D("frame","Correlations",1,0.5,N-0.5);
		Frame->GetXaxis()->SetTitle("IMF");
		Frame->GetYaxis()->SetTitle("Correlation");
		TH1D ** h = new TH1D*[N-1];
		for(int i=1;i<N;i++){
			string cname = "cor" + to_string(i);
			h[i-1] = new TH1D(cname.c_str(),"Correlations",1,i-width,i+width);
			h[i-1]->Fill(i,correlation[i-1]);
			h[i-1]->SetFillColor(i+40);
			h[i-1]->SetMaximum(1);
			lc->AddEntry(h[i-1],to_string(correlation[i-1]).c_str(),"f");
		}

		double Xmin = IMFs[0].vectorX().front();
		double Xmax = IMFs[0].vectorX().back();
		double Ymin = gimf[0]->GetYaxis()->GetXmin();
		double Ymax = gimf[0]->GetYaxis()->GetXmax();
		//cFCgraphics
		cFCgraphics G;

		//Create Pads
		TPad *pad1 = G.CreatePad("pad1");
		TPad *pad2 = G.CreatePad("pad2");

		//Edit First Pad (Signal, Residue and Correlations)
		pad1->Divide(1,3);
		pad1->cd(1);
		gimf[0]->SetNameTitle("Signal","Signal");
		gimf[0]->GetXaxis()->SetLimits(Xmin,Xmax);
		gimf[0]->Draw();

		pad1->cd(2);
		gres->SetNameTitle("Residue","Residue");
		gres->GetXaxis()->SetLimits(Xmin,Xmax);
		gres->GetXaxis()->SetTitle("Time");
		gres->GetYaxis()->SetTitle("Amplitude");
		gres->GetHistogram()->SetMinimum(Ymin);
		gres->GetHistogram()->SetMaximum(Ymax);
		gres->Draw();

		pad1->cd(3);
		Frame->Draw("HIST");
		lc->Draw();
		for(int i=0;i<N-1;i++){
			h[i]->Draw("HIST SAME");
		}

		//Edit Second Pad (IMFs)
		pad2->Divide(1,N-1);
		for(int i=1;i<N;i++){
			pad2->cd(i);
			string name = "IMF " + to_string(i);
			gimf[i]->SetNameTitle("IMF",name.c_str());
			gimf[i]->GetXaxis()->SetLimits(Xmin,Xmax);
			gimf[i]->Draw();
		}

		//Add Pads
		G.AddObject(pad1);
		G.AddObject(pad2);

		//DrawCanvas
		G.Update();
		G.Draw();
		if(filename.length()>0) G.Print(filename);

		for(int i=0;i<N-1;i++){
			delete h[i];
		}
		delete [] h;
		delete Frame;
	}

	else if(flag==1)
	{
		int N = IMFs.size();
		int factor = 3000;
		string title = "Hilbert Spectrum (All IMFs)";
		vector<double> T = IMFs[0].vectorX();
		double Xmin = T.front();
		double Xmax = T.back();

		vector<TGraph *> gfreqs;
		vector<TF1 *> gamps;
		for(int i=0; i<N;i++){
			gfreqs.push_back(IMFs[i].GetFrequency());
			gamps.push_back(IMFs[i].GetAMP());
		}
		gStyle->SetPalette(107);
		double step = (Xmax-Xmin)/((double)factor);
		TH2D ** spec = new TH2D*[N-1];
		for(int i=1;i<N;i++){
			string hsname = "hsname" + to_string(i);
			string hstitle = "Hilbert Spectrum - IMF " + to_string(i);
			spec[i-1] = new TH2D(hsname.c_str(),hstitle.c_str(),factor,Xmin,Xmax,factor/4,gfreqs[i]->GetYaxis()->GetXmin(),gfreqs[i]->GetYaxis()->GetXmax());
			spec[i-1]->SetFillColor(i+40);
			spec[i-1]->GetXaxis()->SetTitle("Time");
			spec[i-1]->GetYaxis()->SetTitle("Frequency");
			spec[i-1]->SetContour(255);
			for(int j=0;j<factor;j++)
			{
				spec[i-1]->Fill(Xmin+j*step,gfreqs[i]->Eval(Xmin+j*step),fabs(gamps[i]->Eval(Xmin+j*step)));
			}
		}
		//Getting TGraphs and Setting Colors

		vector<TGraph *> gimf;
		vector<TGraph *> gper;
		for(int i=0; i<N;i++){
			gimf.push_back(IMFs[i].GetGraph());
			gimf[i]->SetLineColor(i+40);
			gimf[i]->SetMarkerColor(i+40);
		}
		//cFCgraphics
		cFCgraphics G;

		//Create Pads
		TPad *pad1 = G.CreatePad("pad1");
		TPad *pad2 = G.CreatePad("pad2");

		//Edit First Pad (IMFs)
		pad1->Divide(1,N-1);
		for(int i=1;i<N;i++){
			pad1->cd(i);
			string name = "IMF " + to_string(i);
			gimf[i]->SetNameTitle("IMF",name.c_str());
			gimf[i]->GetXaxis()->SetLimits(Xmin,Xmax);
			gimf[i]->GetXaxis()->SetTitle("Time");
			gimf[i]->GetYaxis()->SetTitle("Amplitude");
			gimf[i]->Draw();
		}

		//Edit Second Pad (Periods)
		pad2->Divide(1,N-1);
		for(int i=1;i<N;i++){
			pad2->cd(i);
			spec[i-1]->Draw("COLZ0");
		}


		//Add Pads
		G.AddObject(pad1);
		G.AddObject(pad2);

		//DrawCanvas
		G.Update();
		G.Draw();
		if(filename.length()>0) G.Print(filename);
		for(int i=0;i<N-1;i++){
			delete spec[i];
		}
		delete [] spec;
	}

	else
	{
		//Getting TGraphs and Setting Colors
		int N = IMFs.size();
		vector<TGraph *> gimf;
		vector<TGraph *> gimfnotitle;
		vector<TGraph *> gfreqs;
		vector<TF1 *> gamps;
		for(int i=0; i<N;i++){
			gimf.push_back(IMFs[i].GetGraph());
			gimfnotitle.push_back(IMFs[i].GetGraph());
			gfreqs.push_back(IMFs[i].GetFrequency());
			gamps.push_back(IMFs[i].GetAMP());
			gimf[i]->SetLineColor(i+40);
			gimf[i]->SetMarkerColor(i+40);
			gimf[i]->GetXaxis()->SetTitle("Time");
			gimf[i]->GetYaxis()->SetTitle("Amplitude");
			gimfnotitle[i]->SetLineColor(i+40);
			gimfnotitle[i]->SetMarkerColor(i+40);
			gimfnotitle[i]->GetXaxis()->SetTitle("Time");
			gimfnotitle[i]->GetYaxis()->SetTitle("Amplitude");
		}
		TGraph * gres = residue.GetGraph();
		gres->GetXaxis()->SetTitle("Time");
		gres->GetYaxis()->SetTitle("Amplitude");
		vector<double> T = IMFs[0].vectorX();
		double Xmin = T.front();
		double Xmax = T.back();
		double Ymin = gimf[0]->GetYaxis()->GetXmin();
		double Ymax = gimf[0]->GetYaxis()->GetXmax();

		gStyle->SetPalette(107);
		int factor = 3000;
		double step = (Xmax-Xmin)/((double)factor);
		TH2D ** spec = new TH2D*[N-1];
		for(int i=1;i<N;i++){
			string hsname = "hsname" + to_string(i);
			spec[i-1] = new TH2D(hsname.c_str(),"Hilbert Spectrum",factor,Xmin,Xmax,factor/4,0,gfreqs[1]->GetYaxis()->GetXmax());
			spec[i-1]->SetContour(255);
			spec[i-1]->SetFillColor(i+40);
			spec[i-1]->GetXaxis()->SetTitle("Time");
			spec[i-1]->GetYaxis()->SetTitle("Frequency");
			for(int j=0;j<factor;j++)
			{
				spec[i-1]->Fill(Xmin+j*step,gfreqs[i]->Eval(Xmin+j*step),fabs(gamps[i]->Eval(Xmin+j*step)));
			}
		}

		vector <double> correlation = GetCorrelation();
		TLegend * lc = new TLegend(0.7,0.7,0.95,0.95);
		double width = 0.25;
		TH1D * Frame = new TH1D("frame","Correlations",1,0.5,N-0.5);
		Frame->GetXaxis()->SetTitle("IMF");
		Frame->GetYaxis()->SetTitle("Correlation");
		TH1D ** h = new TH1D*[N-1];
		for(int i=1;i<N;i++){
			string cname = "cname" + to_string(i);
			h[i-1] = new TH1D(cname.c_str(),"Correlations",1,i-width,i+width);
			h[i-1]->Fill(i,correlation[i-1]);
			h[i-1]->SetFillColor(i+40);
			h[i-1]->SetMaximum(1);
			h[i-1]->GetXaxis()->SetTitle("IMF");
			h[i-1]->GetYaxis()->SetTitle("Correlation");
			lc->AddEntry(h[i-1],to_string(correlation[i-1]).c_str(),"f");
		}

		//cFCgraphics
		cFCgraphics G;

		//Create Pads
		TPad *pad1 = G.CreatePad("pad1");
		TPad *pad2 = G.CreatePad("pad2");
		TPad *pad3 = G.CreatePad("pad3");

		//Edit First Pad (Signal and Residue)
		pad1->Divide(1,2);
		pad1->cd(1);
		gimf[0]->SetNameTitle("Signal","Signal");
		gimf[0]->GetXaxis()->SetLimits(Xmin,Xmax);
		gimf[0]->Draw();
		pad1->cd(2);
		gres->SetNameTitle("Residue","Residue");
		gres->GetXaxis()->SetLimits(Xmin,Xmax);
		gres->GetHistogram()->SetMinimum(Ymin);
		gres->GetHistogram()->SetMaximum(Ymax);
		gres->Draw();

		//Edit Second Pad (IMFs)
		if(N>10){
			for(int i=1;i<N;i++){
				gimf[i]->SetNameTitle("IMF","All IMFs");
				gimf[i]->GetXaxis()->SetLimits(Xmin,Xmax);
				G.AddObject(gimf[i],"pad2");
			}
		}
		else{
			pad2->Divide(1,N-1);
			for(int i=1;i<N;i++){
				pad2->cd(i);
				string name = "IMF " + to_string(i);
				gimf[i]->SetNameTitle("IMF",name.c_str());
				gimf[i]->GetXaxis()->SetLimits(Xmin,Xmax);
				gimf[i]->Draw();
			}
		}

		//Edit Third Pad (Correlation, Frequencies and SuperImposed IMFs)
		pad3->Divide(1,3);
		pad3->cd(1);
		spec[0]->Draw("COLZ0");
		//SF->Draw("COLZ");
		for(int i=1;i<N-1;i++){
			spec[i]->Draw("COL0 Same");
		}

		pad3->cd(2);
		Frame->Draw("HIST");
		lc->Draw();
		for(int i=0;i<N-1;i++){
			h[i]->Draw("HIST SAME");
		}

		pad3->cd(3);
		gimfnotitle[1]->SetNameTitle("","All IMFs");
		gimfnotitle[1]->Draw();
		for(int i=2; i<N; i++){
			gimfnotitle[i]->SetNameTitle();
			gimfnotitle[i]->Draw("SAME");
		}

		//Add Pads
		G.AddObject(pad1);
		G.AddObject(pad2);
		G.AddObject(pad3);

		//DrawCanvas
		G.Update();
		G.Draw();
		if(filename.length()>0) G.Print(filename);


		for(int i=0;i<N-1;i++){
			delete h[i];
			delete spec[i];
		}
		delete [] spec;
		delete [] h;
		delete Frame;

	}
}

void EMD::DrawIMF(int n, string filename){
	int N = IMFs.size();
	if(n>=0 && n<N)
	{
		TGraph * g = IMFs[n].GetGraph();
		string title;
		if(n==0) title = "Signal";
		else title = "IMF " + to_string(n);
		vector<double> T = IMFs[0].vectorX();
		double Xmin = T.front();
		double Xmax = T.back();
		g->SetLineColor(n+40);
		g->SetLineWidth(2);
		g->SetMarkerColor(n+40);
		g->GetXaxis()->SetLimits(Xmin,Xmax);
		g->GetXaxis()->SetTitle("Time");
		g->GetYaxis()->SetTitle("Amplitude");
		g->SetNameTitle("IMF",title.c_str());

		cFCgraphics G;
		//auto legend = new TLegend(0.1,0.7,0.48,0.9);
		//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
		TPad * pad = G.CreatePad("pad");
		//legend->AddEntry(pad,"Histogram filled with random numbers","f");
		G.AddObject(g,"pad");
		G.AddObject(pad);
		G.Draw();
		if(filename.length()>0) G.Print(filename);
	}
	else if (n==N)
	{
		TGraph * g = residue.GetGraph();
		string title = "Residue";
		vector<double> T = IMFs[0].vectorX();
		double Xmin = T.front();
		double Xmax = T.back();
		g->SetLineColor(n+40);
		g->SetLineWidth(2);
		g->SetMarkerColor(n+40);
		g->GetXaxis()->SetLimits(Xmin,Xmax);
		g->GetXaxis()->SetTitle("Time");
		g->GetYaxis()->SetTitle("Amplitude");
		g->SetNameTitle("IMF",title.c_str());

		cFCgraphics G;
		TPad * pad = G.CreatePad("pad");
		G.AddObject(g,"pad");
		G.AddObject(pad);
		G.Draw();
		if(filename.length()>0) G.Print(filename);
	}
	else{
		string title = "All IMFS";
		vector<double> T = IMFs[0].vectorX();
		double Xmin = T.front();
		double Xmax = T.back();

		vector<TGraph *> gimf;
		for(int i=1; i<N;i++){
			gimf.push_back(IMFs[i].GetGraph());
		}

		cFCgraphics G;
		TPad * pad = G.CreatePad("pad");
		for(int i=0;i<N-1;i++){
			gimf[i]->SetLineColor(i+41);
			gimf[i]->SetMarkerColor(i+41);
			gimf[i]->GetXaxis()->SetLimits(Xmin,Xmax);
			gimf[i]->GetXaxis()->SetTitle("Time");
			gimf[i]->GetYaxis()->SetTitle("Amplitude");
			gimf[i]->SetNameTitle("IMF",title.c_str());
			G.AddObject(gimf[i],"pad");
		}
		G.AddObject(pad);
		G.Draw();
		if(filename.length()>0) G.Print(filename);
	}
}

void EMD::DrawHS(int n, string filename){
	int N = IMFs.size();
	int factor=3000;

	gStyle->SetPalette(107);
	if(n>0 && n<N)
	{
		string title = "Hilbert Spectrum - IMF " + to_string(n);
		vector<double> T = IMFs[0].vectorX();
		double Xmin = T.front();
		double Xmax = T.back();

		TGraph * gfreqs=IMFs[n].GetFrequency();
		TF1 * gamps = IMFs[n].GetAMP();

		double step = (Xmax-Xmin)/((double)factor);
		string hstitle = "Hilbert Spectrum - IMF " + to_string(n);
		string hsname = "hsname" + to_string(n);
		TH2D * spec = new TH2D(hsname.c_str(),hstitle.c_str(),factor,Xmin,Xmax,factor/4,gfreqs->GetYaxis()->GetXmin(),gfreqs->GetYaxis()->GetXmax());
		spec->GetXaxis()->SetTitle("Time");
		spec->GetYaxis()->SetTitle("Frequency");
		spec->SetContour(255);
		for(int j=0;j<factor;j++)
		{
			spec->Fill(Xmin+j*step,gfreqs->Eval(Xmin+j*step),fabs(gamps->Eval(Xmin+j*step)));
		}


		cFCgraphics G;
		TPad * pad = G.CreatePad("pad");
		G.AddObject(spec,"pad","COLZ0");
		G.AddObject(pad);
		G.Draw();
		if(filename.length()>0) G.Print(filename);

		delete spec;
	}

	else{
		string title = "Hilbert Spectrum (All IMFs)";
		vector<double> T = IMFs[0].vectorX();
		double Xmin = T.front();
		double Xmax = T.back();

		vector<TGraph *> gfreqs;
		vector<TF1 *> gamps;
		for(int i=0; i<N;i++){
			gfreqs.push_back(IMFs[i].GetFrequency());
			gamps.push_back(IMFs[i].GetAMP());
		}

		double step = (Xmax-Xmin)/((double)factor);
		TH2D ** spec = new TH2D*[N-1];
		for(int i=1;i<N;i++){
			string hsname = "hsname" + to_string(i);
			spec[i-1] = new TH2D(hsname.c_str(),"Hilbert Spectrum",factor,Xmin,Xmax,factor/4,0,gfreqs[1]->GetYaxis()->GetXmax());
			spec[i-1]->SetFillColor(i+40);
			spec[i-1]->GetXaxis()->SetTitle("Time");
			spec[i-1]->GetYaxis()->SetTitle("Frequency");
			spec[i-1]->SetContour(255);
			for(int j=0;j<factor;j++)
			{
				spec[i-1]->Fill(Xmin+j*step,gfreqs[i]->Eval(Xmin+j*step),fabs(gamps[i]->Eval(Xmin+j*step)));
			}
		}

		cFCgraphics G;
		TPad * pad = G.CreatePad("pad");
		pad->Divide(1,1);
		pad->cd(1);
		spec[0]->Draw("COLZ0");
		for(int i=1;i<N-1;i++){
			spec[i]->Draw("COL0 Same");
		}
		G.AddObject(pad);
		G.Draw();
		if(filename.length()>0) G.Print(filename);

		for(int i=0;i<N-1;i++){
			delete spec[i];
		}
		delete [] spec;
	}
}

void EMD::DrawCorrelations(string filename){
	TLegend * lc = new TLegend(0.7,0.7,0.95,0.95);
	vector <double> correlation = GetCorrelation();
	int N = IMFs.size();
	double width = 0.25;
	TH1D * Frame = new TH1D("frame","Correlations",1,0.5,N-0.5);
	Frame->GetXaxis()->SetTitle("IMF");
	Frame->GetYaxis()->SetTitle("Correlation");
	TH1D ** h = new TH1D*[N-1];
	for(int i=1;i<N;i++){
		string cname = "cname" + to_string(i);
		h[i-1] = new TH1D(cname.c_str(),"Correlations",1,i-width,i+width);
		h[i-1]->Fill(i,correlation[i-1]);
		h[i-1]->SetFillColor(i+40);
		h[i-1]->SetMaximum(1);
		lc->AddEntry(h[i-1],to_string(correlation[i-1]).c_str(),"f");
	}

	cFCgraphics G;
	TPad * pad = G.CreatePad("pad");
	pad->Divide(1,1);
	pad->cd(1);
	Frame->Draw("HIST");
	lc->Draw();
	for(int i=0;i<N-1;i++){
		h[i]->Draw("HIST SAME");
	}
	G.AddObject(pad);

	//auto legend = new TLegend(0.1,0.7,0.48,0.9);
	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	//legend->AddEntry(pad,"Histogram filled with random numbers","f");

	//G.AddObject(legend);

	G.Draw();
	if(filename.length()>0) G.Print(filename);


	for(int i=0;i<N-1;i++){
		delete h[i];
	}
	delete [] h;
	delete Frame;
}

void EMD::DrawCF(string filename){
	int N = IMFs.size();
	double width = 0.25;
	TLegend * lc = new TLegend(0.7,0.7,0.95,0.95);
	TLegend * lf = new TLegend(0.7,0.7,0.95,0.95);


	vector <double> correlation = GetCorrelation();
	TH1D * Frame = new TH1D("frame","Correlations",1,0.5,N-0.5);
	Frame->GetXaxis()->SetTitle("IMF");
	Frame->GetYaxis()->SetTitle("Correlation");
	TH1D ** h = new TH1D*[N-1];
	for(int i=1;i<N;i++){
		string cname = "cname" + to_string(i);
		h[i-1] = new TH1D(cname.c_str(),"Correlations",1,i-width,i+width);
		h[i-1]->Fill(i,correlation[i-1]);
		h[i-1]->SetFillColor(i+40);
		h[i-1]->SetMaximum(1);
		lc->AddEntry(h[i-1],to_string(correlation[i-1]).c_str(),"f");
	}

	vector <double> avf = GetFrequencies();
	TH1D * Framef = new TH1D("framef","Average Frequencies",1,0.5,N-0.5);
	Framef->GetXaxis()->SetTitle("IMF");
	Framef->GetYaxis()->SetTitle("Frequency");
	Framef->SetMaximum(avf[0]*1.1);
	TH1D ** fr = new TH1D*[N-1];
	for(int i=1;i<N;i++){
		string cname = "fname" + to_string(i);
		fr[i-1] = new TH1D(cname.c_str(),"Frequencies",1,i-width,i+width);
		fr[i-1]->Fill(i,avf[i-1]);
		fr[i-1]->SetFillColor(i+40);
		lf->AddEntry(fr[i-1],to_string(avf[i-1]).c_str(),"f");
	}

	cFCgraphics G;
	TPad * pad = G.CreatePad("pad");
	pad->Divide(2,1);
	pad->cd(1);
	Frame->Draw("HIST");
	for(int i=0;i<N-1;i++){
		h[i]->Draw("HIST SAME");
	}
	lc->Draw();

	pad->cd(2);
	Framef->Draw("HIST");
	for(int i=0;i<N-1;i++){
		fr[i]->Draw("HIST SAME");
	}
	lf->Draw();

	G.AddObject(pad);





	G.Draw();
	if(filename.length()>0) G.Print(filename);


	for(int i=0;i<N-1;i++){
		delete h[i];
	}
	delete [] h;
	delete Frame;

	for(int i=0;i<N-1;i++){
		delete fr[i];
	}
	delete [] fr;
	delete Framef;
}
