#ifndef __IMF__
#define __IMF__

#include <TGraph.h>
#include <TF1.h>
#include <vector>
using namespace std;

class IMF {
public:
  IMF(double * x = NULL, double * y = NULL, int n=0, int R=10000); //Default constructor and constructor from double *
  IMF(vector<double> x, vector<double> y, int R=10000); // Constructor from vectors
  IMF(TF1 * F, double xmin, double xmax, int n, int R=10000); // Constructor from TF1
  IMF(const IMF &);// Copy Constructor
  ~IMF(); //Destructor

  //operators
  IMF& operator = (const IMF &); // copy assignment
  IMF operator + (const IMF &); // summing two IMF's
  IMF operator - (const IMF & ); //taking the difference between two IMF's


  //methods
  void SetPoints(double * x, double * y, int n); //Set points from double *
  void SetPoints(vector<double> x, vector<double> y); // Set points from vectors
  double * GetX(); //return x points in double *
  double * GetY(); //return y points in double *
  int GetN() const{return N;} //return number of points
  int GetResolution() const{return Resolution;} //return resolution
  int GetEN() const {return EN;} //return extrema number
  vector<double> vectorX() const; //return x points in vector
  vector<double> vectorY() const; //return y points in vector
  vector<double> GetRoots(); //return zeros in vector
  vector<vector<double> > GetMaxs(); //return maxs in vector
  vector<vector<double> > GetMins(); //return mins in vector
  TGraph * GetGraph(); //retrieve graph of IMF
  TGraph * GetFrequency(); //returns the instantaneous frequency of the nth IMF
  TF1* GetMAX() const {return MAX;} //return MAX envelope function
  TF1* GetMIN() const {return MIN;} //return MIN envelope function
  TF1* GetMED() const {return MED;} //return MED function
  TF1* GetAMP() const {return AMP;} //return AMP function
  int isIMF(); //is IMF test
  int isLastIMF(); //Test if residue has IMFs
  void GetDetail(); // GetDetail
private:
  //methods
  vector<vector<double> > Mirror(vector<double>x, vector<double>y, int flag); //Mirroring method returns pointer to vectors ([0]->x,[1]->y)
  void SetMinimum(); //return vector with interpolation on minimums
  void SetMaximum(); //return vector with interpolation on maximums
  void SetMedium(); //returns the medium vector
  void SetAmplitude(); //returns the amplitude vector
  double Medium(double*fx, double*par); //MED EVAL
  double Amplitude(double*fx, double*par); //AMP EVAL

  //data members
  int N; //number of points
  int EN; //number of extrema
  int Resolution; //resolution of splines
  vector<double> X; //x points
  vector<double> Y; //y points
  TF1* MAX; //MAX envelope function
  TF1* MIN; //MIN envelope function
  TF1* MED; //MED function
  TF1* AMP; //AMP function
};
#endif
