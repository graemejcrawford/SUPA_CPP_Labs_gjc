// Author: Graeme Crawford
// Date: 9/12/23
// Header File for the plotting of functions and Metropolis sampling

#include <string>
#include <vector>
#include "gnuplot-iostream.h"

#pragma once //Replacement for IFNDEF

void Read(std::string file, std::vector<double>& Values);  // Function to read in data

/*
###################
//Base class (default function)
###################
*/ 

class FiniteFunction{

public:
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor
  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  double integral(int Ndiv = 1000); 
  std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setOutfile(std::string outfile);
  void plotFunction(); //Plot the function using scanFunction
    
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)
  virtual void Sampling(); // GC: Performs Metropolis algorithm and plots the sampled data (Overridable)
  

//Protected members can be accessed by child classes but not users
protected:
  double m_RMin;
  double m_RMax;
  double m_Integral;
  int m_IntDiv = 0; //Number of division for performing integral
  std::string m_FunctionName;
  std::string m_OutData; //Output filename for data
  std::string m_OutPng; //Output filename for plot
  std::vector< std::pair<double,double> > m_data; //input data points to plot
  std::vector< std::pair<double,double> > m_samples; //Holder for randomly sampled data 
  std::vector< std::pair<double,double> > m_function_scan; //holder for data from scanFunction (slight hack needed to plot function in gnuplot)
  bool m_plotfunction = false; //Flag to determine whether to plot function
  bool m_plotdatapoints = false; //Flag to determine whether to plot input data
  bool m_plotsamplepoints = false; //Flag to determine whether to plot sampled data 
  double integrate(int Ndiv); // GC: Uses the left Riemann sum to approximate intergal of distribution
  std::vector< std::pair<double, double> > makeHist(std::vector<double> &points, int Nbins); //Helper function to turn data points into histogram with Nbins
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 
  std::vector<double> m_sampleData; // GC: Vector for storing the accepted sample data
  double m_randomX[10000]; // GC: Run up to 10000 samples
  double m_randomY[10000];
  double m_randomT[10000];
  
private:
  double invxsquared(double x); //The default functional form
};


/*
###################
//Inherited Classes
###################
*/ 


/*
###################
//Normal Distribution
###################
*/ 

class normalDistribution: public FiniteFunction{

public:
  normalDistribution(); //Empty constructor
  normalDistribution(double range_min, double range_max, std::string outfile, double mean_val, double stdev_val); // GC: Updated variable constructor
  ~normalDistribution(); //Destructor

  double mean_val(); // GC: Mean value of the normal distribution
  double stdev_val(); // GC: Standard deviation of the normal distribution
  void setMean(double mean);
  void setStdev(double stdev);
  void plotData(std::vector<double> &points, int NBins, bool isdata=true);
  double integrate(int Ndiv);
  double integral(int Ndiv); 
  void plotFunction();
  void checkPath(std::string outstring); 
  void generatePlot(Gnuplot &gp); 
  virtual void printInfo(); // GC: Print important info about the normal distribution 
  virtual double callFunction(double x);  // GC: Call normal distribution function at a specifed x value
  virtual void Sampling(); // GC: Using same sampling method as default function, requiring no repetition of code
  
  

protected: // new protected members for normal distribution
  double m_mean;
  double m_stdev;


private:
 double normDist(double x); // GC: Keep the functional form of the normal distribution private

};


/*
###################
//Cauchy Lorentz
###################
*/ 

class cauchyLorentz: public FiniteFunction{



public:
  cauchyLorentz(); //Empty constructor
  cauchyLorentz(double range_min, double range_max, std::string outfile, double peakLoc_val, double gamma_val); // GC: Updated variable constructor
  ~cauchyLorentz(); //Destructor

  double peakLoc_val(); // GC: location of the peak of the Cauchy Lorentz distribution
  double gamma_val(); // GC: value for the scale parameter gamma
  void setPeakLoc(double peakLoc);
  void setGamma(double gamma); 
  void plotData(std::vector<double> &points, int NBins, bool isdata=true);
  double integrate(int Ndiv);
  double integral(int Ndiv); 
  void plotFunction();
  void checkPath(std::string outstring); 
  void generatePlot(Gnuplot &gp); 
  virtual void printInfo(); 
  virtual double callFunction(double x); 
  virtual void Sampling();
 
  

protected: // new protected members for Cauchy Lorentz distribution

  double m_peakLoc;
  double m_gamma;


private:
// GC: Keep the functional form of the Cauchy Lorentz distribution private
 double cauchyLorentzDistribution(double x); 

};


/*
###################
//Negative Crystal Ball
###################
*/ 

class negativeCrystalBall: public FiniteFunction{

public:
  negativeCrystalBall(); //Empty constructor
  negativeCrystalBall(double range_min, double range_max, std::string outfile,  double mean_val, double stdev_val, double n_val, double alpha_val); // GC: Updated variable constructor
  ~negativeCrystalBall(); //Destructor

  double mean_val(); // GC: Mean value for the negative crystal ball distribution
  double stdev_val(); // GC: Standard deviation value for the negative crystal ball distribution
  double n_val(); // GC: n is a parameter which determines the shape
  double alpha_val(); // GC: alpha is another parameter which affects the shape of the distribution
  void setMean(double mean);
  void setStdev(double stdev);
  void setN(double n);
  void setAlpha(double alpha);
  double integrate(int Ndiv);
  double integral(int Ndiv); 
  void plotData(std::vector<double> &points, int NBins, bool isdata=true);
  void plotFunction();
  void checkPath(std::string outstring); 
  void generatePlot(Gnuplot &gp); 
  virtual void printInfo();  
  virtual double callFunction(double x); 
  virtual void Sampling(); 

// GC: Functions required to build up the negative Crystal Ball distribution
  double A();
  double B();
  double N(); // GC: Constant required for normalisation
  double C();
  double D();
  

protected: // new protected members for negatove Crystal Ball distribution

  double m_mean;
  double m_stdev;
  double m_n;
  double m_alpha;


private:
// GC: Keep the functional form of the negative Crystal Ball distribution private
 double negativeCrystalBallDistribution(double x); 

};
