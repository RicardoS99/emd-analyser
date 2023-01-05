#ifndef __VEC__
#define __VEC__

class Vec {
public:
  Vec();//create a null vector
  Vec(int); //create null vector
  Vec(int, double); //create vector with a value
  Vec(int, double*); //create vector from array
  Vec(const Vec &); //create vector from vector
  ~Vec();
  void SetEntries (int, double*);
  int size() const; //return number of elements
  double & operator[](int); //return
  void print(); // print vector
  double At(int) const; // access entrie[i]
  Vec& Add(double); // similar to push_back method of vector class
  Vec& Del(int); // delete one element
  Vec& Del(int, int); // delete elements between two indexes
  Vec& operator = (const Vec &); // copy constructor
  Vec& operator += (const Vec &); // sum Vec
  Vec& operator -= (const Vec &); //subtract Vec
  Vec operator + (const Vec &); // return Vec of sum of two Vecs
  Vec operator - (const Vec &); // return Vec of difference of two Vecs
  Vec operator - (); // swap sign of all elements of vector
  Vec operator * (const Vec &); // return multiplication of two Vecs
  Vec operator * (double); // return multiplication of a Vec by a scalar
  double dot(const Vec &); // return dot product
  void swap(int, int); // swap elements of vector
  double norm(); //return Euclidian Norm
  double norm(int); // return Norm n (if n<=0, return Infinity Norm)

private:
  int N; //number of elements
  double *entries; // pointer to array of doubles
};

#endif
