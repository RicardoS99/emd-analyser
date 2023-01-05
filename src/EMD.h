#include <vector>
#include <string>
#include "IMF.h"
#include "TF1.h" //for the cubic splines
using namespace std;

class EMD{
public:
	//constructors and destructors
	EMD (){;} //default constructor (does nothing)
	EMD (IMF fsignal); //construct and decompose
	~EMD(){;} // destructor

	//methods
	IMF GetIMF(int n = 0); //return nth IMF (IMF[0] is the signal)
	IMF GetResidual(){return residue;} // return the residue
	vector<IMF> GetIMFs() {return IMFs;}
	vector <double> GetCorrelation() {return correlation;}
	vector <double> GetFrequencies() {return frequencies;}
	void GetAnalysis(int flag=-1,string filename = ""); //return the analysis on graphics (flag determines tha analysis presented)
	void DrawIMF(int n=-1, string filename = "");
	void DrawHS(int n=-1, string filename = "");
	void DrawCorrelations(string filename = "");
	void DrawCF(string filename = "");

private:
	//data members
	vector <IMF> IMFs; //vector with signal + IMFs
	IMF residue; // residual of decomposition
	vector<double> correlation;
	vector<double> frequencies;
	//methods
	void IMFanalyser ();	//constructs the IMF vector
	double CalculateCorrelation(int n); //returns the correlation of the nth IMF
	double CalculateFrequency(int n); //returns the correlation of the nth IMF
};
