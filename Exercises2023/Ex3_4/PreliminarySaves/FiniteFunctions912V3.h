// Author: Graeme Crawford
// Date: 9/12/23
// Header File for the plotting of functions and sampling

#include <string>
#include <vector>
#include "gnuplot-iostream.h"

#pragma once //Replacement for IFNDEF

class FiniteFunction{

public:
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor
  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  double integral(int Ndiv = 1000); 
  std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  void setOutfile(std::string outfile);
  void plotFunction(); //Plot the function using scanFunction
    
 
  
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)
  virtual void Sampling(); // (Overridable)

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
  double integrate(int Ndiv);
  std::vector< std::pair<double, double> > makeHist(std::vector<double> &points, int Nbins); //Helper function to turn data points into histogram with Nbins
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 
  std::vector<double> m_sampleData;
  double m_randomX[1000];
  double m_randomY[1000];
  double m_randomT[1000];
  
private:
  double invxsquared(double x); //The default functional form
};


// Making an inherited class for a normal distribution

class normalDistribution: public FiniteFunction{

public:
  normalDistribution(); //Empty constructor
  normalDistribution(double range_min, double range_max, std::string outfile, double mean_val, double stdev_val); //Variable constructor
  ~normalDistribution(); //Destructor

  double mean_val(); 
  double stdev_val(); 
  void setMean(double mean);
  void setStdev(double stdev);
  virtual void printInfo(); //Dump parameter info about the current function (Overridable) 
  virtual double callFunction(double x); //Call the function with value x (Overridable) 
  void plotData(std::vector<double> &points, int NBins, bool isdata=true);
  double integrate(int Ndiv);
  void plotFunction();
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 
  double integral(int Ndiv); 
  virtual void Sampling();
  

protected:

  double m_mean;
  double m_stdev;


private:

 double normDist(double x);

};


// Making an inherited class for the Cauchy-Lorentz function

class cauchyLorentz: public FiniteFunction{

public:
  cauchyLorentz(); //Empty constructor
  cauchyLorentz(double range_min, double range_max, std::string outfile, double peakLoc_val, double gamma_val); //Variable constructor
  ~cauchyLorentz(); //Destructor

  double peakLoc_val(); 
  double gamma_val(); 
  void setPeakLoc(double peakLoc);
  void setGamma(double gamma);
  //void setRandomX(int randomX);
  //void setRandomY(int randomY);
  virtual void printInfo(); //Dump parameter info about the current function (Overridable) 
  virtual double callFunction(double x); //Call the function with value x (Overridable) 
  void plotData(std::vector<double> &points, int NBins, bool isdata=true);
  virtual void Sampling();
  double integrate(int Ndiv);
  void plotFunction();
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 
  double integral(int Ndiv); 
 
  

protected:

  double m_peakLoc;
  double m_gamma;



   



 
  //double m_mean_sampling;
  //double m_width;
  // double m_randomY;

    // This should hopefully allow the sampling to be done for Cauchy Lorentz, should be accessible to other classes

  // double m_mean_sampling;
  // double m_width;
  // double m_randomY;
  // double m_fy;
  // double m_fxi;
  // double m_Tmin;
  // double m_Tmax;
  // double m_randomT;


private:

 double cauchyLorentzDistribution(double x);

};


// Making an inherited class for the Cauchy-Lorentz function

class negativeCrystalBall: public FiniteFunction{

public:
  negativeCrystalBall(); //Empty constructor
  negativeCrystalBall(double range_min, double range_max, std::string outfile,  double mean_val, double stdev_val, double n_val, double alpha_val); //Variable constructor
  ~negativeCrystalBall(); //Destructor

  double mean_val(); 
  double stdev_val(); 
  double n_val();
  double alpha_val();
  void setMean(double mean);
  void setStdev(double stdev);
  void setN(double n);
  void setAlpha(double alpha);
  virtual void printInfo(); //Dump parameter info about the current function (Overridable) 
  virtual double callFunction(double x); //Call the function with value x (Overridable) 
  void plotData(std::vector<double> &points, int NBins, bool isdata=true);
  double A();
  double B();
  double N();
  double C();
  double D();
  double integrate(int Ndiv);
  void plotFunction();
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 
  double integral(int Ndiv); 
  virtual void Sampling();
  

protected:

  double m_mean;
  double m_stdev;
  double m_n;
  double m_alpha;


private:

 double negativeCrystalBallDistribution(double x);

};
