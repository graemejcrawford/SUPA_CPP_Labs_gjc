// Author: Graeme Crawford
// Date: 26.11.23

// This is the header file which contains the prototypes for the custom functions



#pragma once // Replaces the #ifndef and #define lines 

#include <vector>
#include <string>

using namespace std;

void Read(string filename, vector<double>& xValues, vector<double>& yValues);

void Mag( double xData, double yData, vector<double>& vecMag); 

void Lsq( vector<double>& xData, vector<double>& yData, vector<double>& errMag, double yErrors, vector<double>& yExp, double& redChisq, double& p, double& q);

void Power(double xData, double yData, vector<double>& power);

void print(vector<double> xData, vector<double>yData);

void print(vector<double> magnitudes);

void print(double p, double q, double chisq);

void printxy(vector<double> XY);

void printBestFit(vector<double> xData, vector<double> yExp);



