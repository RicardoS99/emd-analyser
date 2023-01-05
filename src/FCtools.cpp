#include "FCtools.h"
#include "Vec.h"
#include "TGraph.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

vector<string> FCtools::StrReadFile (string file_name)
{
  fstream file(file_name.c_str());
  vector<string> output;

  //Check if file exists
  if(file.good())
  {
    string line;
    while (getline(file,line))
    {
      output.push_back(line);
      line.clear();
    }
    if(file.eof())
    {
      cout << "\nFile Reading Successful!\n";
    }
    else
    {
      cout << "\nError Reading File!\n";
    }
    file.close();
    return output;
  }
  else
  {
    cout << "\nFile not found!\n";
    file.close();
    return output;
  }
  file.close();
  cout << "\nSomething went wrong!\n";
  return output;
}



vector<vector<double> > FCtools::VecReadFile(string file_name)
{
  fstream file(file_name.c_str());
  vector<vector<double> > output;
  //Check if file exists
  if(file.good())
  {
    vector<string> lines = FCtools::StrReadFile(file_name.c_str()); //Initialize a vector with lines of file as strings

    int n = lines.size(); // return number of lines

    for(int i=0; i<n;i++)
    {
      if(lines[i].find("//",0)==string::npos)
      {
        stringstream ss;
        ss << lines[i]; //line to stringstream
        double element;
        vector<double> temp;
        while(ss>>element) //retrieve elements from stringstream
        {
          temp.push_back(element); //store elements in Vec
        }
        output.push_back(temp); //Add Vec to ouput vector<Vec>
      }
    } 
    file.close();
    return output;
  }
  else
  {
    cout << "\nFile not found!\n";
    file.close();
    return output;
  }
  file.close();
  cout << "\nSomething went wrong!\n";
  return output;
}