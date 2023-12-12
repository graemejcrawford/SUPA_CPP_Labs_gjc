// Author: Graeme Crawford
// Date: 9/12/23
// Functions for the plotting of test distributions as well as executing Metropolis sampling

// Comments added to the initial script will be denoted by "GC: "


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <numbers>
#include "FiniteFunctionsFinal.h"
#include <filesystem> //  To check extensions in a nice way
#include <random> // GC: To generate random numbers
#include <vector>

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 

using std::filesystem::path;


/*
###################
//Reading data
###################
*/ 


// GC: Function to read in data from a file

void Read(std::string file, std::vector<double>& Values){

    std::ifstream mysteryData(file); // GC: Using ifstream to read in the mystery data
   
    
    if (mysteryData.fail()) { // GC: Checking if file has opened successfully
        std::cout << "File did not open" << std::endl;
        exit(1);

    }

    double Data; // GC: Variable to hold the mystery data
   
    while(!mysteryData.eof() ) { 
        mysteryData>> Data;
        Values.push_back(Data); // GC: Storing data in a vector
       
        }

   mysteryData.close(); // Closed after finishing

    }


/*
###################
// Default Function
###################
*/ 

//Empty constructor
FiniteFunction::FiniteFunction(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("Default");
  m_Integral = NULL;
}

//initialised constructor
FiniteFunction::FiniteFunction(double range_min, double range_max, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  this->checkPath(outfile); // Use provided string to name output files
}

//Plots are called in the destructor
//SUPACPP note: They syntax of the plotting code is not part of the course
FiniteFunction::~FiniteFunction(){
  Gnuplot gp; //Set up gnuplot object
  this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}

/*
###################
//Setters
###################
*/ 
void FiniteFunction::setRangeMin(double RMin) {m_RMin = RMin;};
void FiniteFunction::setRangeMax(double RMax) {m_RMax = RMax;};
void FiniteFunction::setOutfile(std::string Outfile) {this->checkPath(Outfile);};

/*
###################
//Getters
###################
*/ 
double FiniteFunction::rangeMin() {return m_RMin;};
double FiniteFunction::rangeMax() {return m_RMax;};

/*
###################
//Function eval
###################
*/ 
double FiniteFunction::invxsquared(double x) {return 1/(1+x*x);};
double FiniteFunction::callFunction(double x) {return this->invxsquared(x);}; //(overridable)

/*
###################
Integration by hand (output needed to normalise function when plotting)
###################
*/ 

// GC: Integration will be performed via left Riemann sum i.e finding area under curve by summing over rectangles using the left x point in each interval as the rectangle height

double FiniteFunction::integrate(int Ndiv){ //private
double sum = 0; // GC: Define a varible to calculate the sum 
double widths =  (this->rangeMax() - this->rangeMin()) / Ndiv ; // GC: Width of each rectangle determined by number of intervals
  for (int i = 0; i<Ndiv; i++) // GC: The last rectangle starts at Ndiv - 1, 
  {
    double xEvals = this->rangeMin() + (i*widths); // GC: Starting from the minimum the interval step size is the width of a rectangle
    double heights = this->callFunction(xEvals); // GC: Now call the function to get the height of the rectangle at each step
    double rectArea = heights * widths; // GC: Finding the area of each rectangle 

    sum += rectArea; // GC: Integral approximated to the sum of all rectangle areas
    
  }

  std::cout << "Width of each rectangle is:  " << widths << std::endl;
  return sum;

}
double FiniteFunction::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}

/*
###################
//Helper functions 
###################
*/
// Generate paths from user defined stem
void FiniteFunction::checkPath(std::string outfile){
 path fp = outfile;
 m_FunctionName = fp.stem(); 
 m_OutData = m_FunctionName+".data";
 m_OutPng = m_FunctionName+".png";
}


/*
###################
//Sampling 
###################
*/

// GC: Sampling function via Metropolis algortihm

void FiniteFunction::Sampling(){

    // GC: Generating initial random x value
   
    int random_num = 1; // GC: Want to generate one initial random number for x
    std::random_device rd; // GC: Used for initial seed
    std::mt19937 mtEngine{rd()}; // GC: Using the Mersenne twister random number generator
    std::uniform_real_distribution<float> rndNumber{m_RMin,m_RMax}; // GC: Random number sampled from uniform distribution between set min and max
    for (int j=0; j < random_num; j++){
    double rndX = rndNumber(mtEngine);
    m_randomX[j] = rndX;  

    }
     
  
    for (int i=0; i < 10000; i++){ // GC: Now want to sample over many different values

    // GC: Generating random Y values from normal distribution

    float m_centre_norm = m_randomX[i]; // GC: Random y value taken from normal distribution centred on x
    float m_width = 1.2; // GC: Arbitrarily chosen standard deviation for normal distribution, tweaked to find most suitable value
    std::normal_distribution<float> normalPDF{m_centre_norm,m_width}; 
    double rndY = normalPDF(mtEngine);
    m_randomY[i] = rndY;

    // GC: Generating a random value T between 0 and 1

    int m_Tmin = 0; // GC: Lower bound for T
    int m_Tmax = 1; // GC: Upper bound for T
    std::uniform_real_distribution<float> rndTNumber{m_Tmin,m_Tmax}; // GC: Random number generated from a uniform distribution
    double rndT = rndTNumber(mtEngine);
    m_randomT[i] = rndT;

    // GC: Now defining the acceptance ratio A

    double A;

   if (callFunction(m_randomY[i])/callFunction(m_randomX[i]) < 1){ // GC: Subbing in sample values for x and y into the test function

        A = callFunction(m_randomY[i])/callFunction(m_randomX[i]);
    }

    else { 

        A = 1; // GC: Acceptance ratio of 100%
    }

    if (m_randomT[i] < A){ // GC: Random value T lower than accepatnce ratio

      m_randomX[i+1] = m_randomY[i]; // GC: Accepts y and sets x(i+1) to y(i) in the next iteration of the loop      
    }

    else{

      m_randomX[i+1] = m_randomX[i]; // GC: Rejetcs y and sets x(i+1) to x(i) in the next iteration of the loop  
    }
        m_sampleData.push_back(m_randomX[i+1]); // GC: Storing all accepted values
        this->plotData(m_sampleData,105,false); // GC: Using the pre-defined plotting functionality to plot sampled data using 105 bins, avoiding troublesome bin numbers(50,70,90,110)
    }    
 }


//Print (overridable)
void FiniteFunction::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
  std::cout << std::endl; 
}



/*
###################
//Plotting
###################
*/

//Hack because gnuplot-io can't read in custom functions, just scan over function and connect points with a line... 
void FiniteFunction::plotFunction(){
  m_function_scan = this->scanFunction(10000);
  m_plotfunction = true;
}

//Transform data points into a format gnuplot can use (histogram) and set flag to enable drawing of data to output plot
//set isdata to true (default) to plot data points in black, set to false to plot sample points in blue
void FiniteFunction::plotData(std::vector<double> &points, int Nbins, bool isdata){
  if (isdata){
    m_data = this->makeHist(points,Nbins);
    m_plotdatapoints = true;
  }
  else{
    m_samples = this->makeHist(points,Nbins);
    m_plotsamplepoints = true;
  }
}


/*
  #######################################################################################################
  ## SUPACPP Note:
  ## The three helper functions below are needed to get the correct format for plotting with gnuplot
  ## In theory you shouldn't have to touch them
  ## However it might be helpful to read through them and understand what they are doing
  #######################################################################################################
 */

//Scan over range of function using range/Nscan steps (just a hack so we can plot the function)
std::vector< std::pair<double,double> > FiniteFunction::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

//Function to make histogram out of sampled x-values - use for input data and sampling
std::vector< std::pair<double,double> > FiniteFunction::makeHist(std::vector<double> &points, int Nbins){
  std::vector< std::pair<double,double> > histdata; //Plottable output shape: (midpoint,frequency)
  std::vector<int> bins(Nbins,0); //vector of Nbins ints with default value 0 
  int norm = 0;
  for (double point : points){
    //Get bin index (starting from 0) the point falls into using point value, range, and Nbins
    int bindex = static_cast<int>(floor((point-m_RMin)/((m_RMax-m_RMin)/(double)Nbins)));
    if (bindex<0 || bindex>Nbins){ // GC: Updated to fix initial bug
      continue;
    }
    bins[bindex]++; //weight of 1 for each data point
    norm++; //Total number of data points
  }
  double binwidth = (m_RMax-m_RMin)/(double)Nbins;
  for (int i=0; i<Nbins; i++){
    double midpoint = m_RMin + i*binwidth + binwidth/2; //Just put markers at the midpoint rather than drawing bars
    double normdata = bins[i]/((double)norm*binwidth); //Normalise with N = 1/(Ndata*binwidth)
    histdata.push_back(std::make_pair(midpoint,normdata));
  }
  return histdata;
}





//Function which handles generating the gnuplot output, called in destructor
//If an m_plot... flag is set, the we must have filled the related data vector
//SUPACPP note: They syntax of the plotting code is not part of the course
void FiniteFunction::generatePlot(Gnuplot &gp){

  if (m_plotfunction==true && m_plotdatapoints==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotdatapoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
  }
  else if (m_plotfunction==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title 'function'\n";
    gp.send1d(m_function_scan);
  }

  else if (m_plotdatapoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_data);
  }

  else if (m_plotsamplepoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_samples);
  }
}

/*
###################
// Inherited Classes
###################
*/ 




/*
###################
// Normal distribution
###################
*/ 

//Empty constructor
normalDistribution::normalDistribution(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  m_mean = 0; // GC: Setting the mean of the normal distribution to best match the mystery data
  m_stdev = 0.9; // GC: Setting the standard deviation of the normal distribution to best match the mystery data
  this->checkPath("NormalFunction");
  m_Integral = NULL;
}

  normalDistribution::normalDistribution(double range_min, double range_max, std::string outfile, double mean_val, double stdev_val){ // GC: Overloaded constructor
    FiniteFunction(range_min, range_max, outfile); // GC: Calls the parent class' constructor
    m_mean = mean_val; // GC: Mean value unique to the normal distribution class
    m_stdev = stdev_val; // GC: Standard deviation value unique to the normal distribution class
}

  normalDistribution::~normalDistribution(){
    Gnuplot gp; //Set up gnuplot object
    this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}


/*
###################
//Setters
###################
*/ 

// GC: Defining new setters for the members unique to the normal distribution
void normalDistribution::setMean(double mean) {m_mean = mean;};
void normalDistribution::setStdev(double stdev) {m_stdev = stdev;};




/*
###################
//Getters
###################
*/ 

// GC: Defining new getters for the members unique to the normal distribution
double normalDistribution::mean_val() {return m_mean;};
double normalDistribution::stdev_val() {return m_stdev;};



/*
###################
//Integrating
###################
*/ 

double normalDistribution::integrate(int Ndiv){
    FiniteFunction::integrate(Ndiv);
}

double normalDistribution::integral(int Ndiv) { 
  FiniteFunction::integral(Ndiv); 

}


// Check path

void normalDistribution::checkPath(std::string outfile){
  FiniteFunction::checkPath(outfile);
}


/*
###################
//Plotting
###################
*/ 

// GC: Plotting function

void normalDistribution::plotFunction(){
  FiniteFunction::plotFunction();
}


// Generate Plot

void normalDistribution::plotData(std::vector<double> &points, int Nbins, bool isdata){
  FiniteFunction::plotData(points, Nbins, isdata);
}

void normalDistribution::generatePlot(Gnuplot &gp){
  FiniteFunction::generatePlot(gp);
}


/*
###################
//Function Evals
###################
*/ 

// GC: Defining the normal distribution using new members
const double pi = 2 * acos(0); // GC: defining pi
double normalDistribution::normDist(double x) {return (1/(m_stdev * sqrt(2 * pi))) *  exp(pow((x-m_mean)/m_stdev,2)*(-0.5));}; 
double normalDistribution::callFunction(double x) {return this->normDist(x);}; //(overridable)


/*
###################
//Sampling
###################
*/ 

// GC: Using new sampling method without repeating any code

void normalDistribution::Sampling(){
  FiniteFunction::Sampling();
}

// GC: Printing useful information to the terminal

void normalDistribution::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "mean: " << m_mean << std::endl;
  std::cout << "stdev: " << m_stdev << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
  std::cout << std::endl; 
}



/*
###################
// Cauchy Lorentz
###################
*/ 


  //Empty constructor
cauchyLorentz::cauchyLorentz(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  m_peakLoc = 0; // GC: Setting the value of x0 (location of peak) which best describes the mystery data
  m_gamma = 0.76; // GC: Setting the value of gamma (related to full width half maximum) which best describes the mystery data, must be greater than zero
  this->checkPath("CauchyLorentz");
  m_Integral = NULL;
}



cauchyLorentz::cauchyLorentz(double range_min, double range_max, std::string outfile, double peakLoc_val, double gamma_val){ // GC: Overloaded constructor
    FiniteFunction(range_min, range_max, outfile); // GC: Calls the parent class' constructor
    // GC: Members unique to the Cauchy Lorentz function
    m_peakLoc = peakLoc_val;
    m_gamma = gamma_val;
      
}

cauchyLorentz::~cauchyLorentz(){
Gnuplot gp; //Set up gnuplot object
this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}


/*
###################
//Setters
###################
*/ 

// GC: Defining new setters for the members unique to the Cauchy Lorentz distribution
void cauchyLorentz::setPeakLoc(double peakLoc) {m_peakLoc = peakLoc;};
void cauchyLorentz::setGamma(double gamma) {m_gamma = gamma;};

/*
###################
//Getters
###################
*/ 

// GC: Defining new getters for the members unique to the Cauchy Lorentz distribution
double cauchyLorentz::peakLoc_val() {return m_peakLoc;};
double cauchyLorentz::gamma_val() {return m_gamma;};


/*
###################
//Integrating
###################
*/ 

double cauchyLorentz::integrate(int Ndiv){
    FiniteFunction::integrate(Ndiv);
}

double cauchyLorentz::integral(int Ndiv) { 
  FiniteFunction::integral(Ndiv); 

}


// Check path

void cauchyLorentz::checkPath(std::string outfile){
  FiniteFunction::checkPath(outfile);
}


/*
###################
//Plotting
###################
*/ 

// Plot function

void cauchyLorentz::plotFunction(){
  FiniteFunction::plotFunction();
}

// Generate Plot

void cauchyLorentz::plotData(std::vector<double> &points, int Nbins, bool isdata){
  FiniteFunction::plotData(points, Nbins, isdata);
}


void cauchyLorentz::generatePlot(Gnuplot &gp){
  FiniteFunction::generatePlot(gp);
}


/*
###################
//Function evals
###################
*/ 

// GC: Defining the Cauchy Lorentz distribution using new members
double cauchyLorentz::cauchyLorentzDistribution(double x) {return (1/(pi * m_gamma *(1 + pow((x-m_peakLoc)/m_gamma,2))));};
double cauchyLorentz::callFunction(double x) {return this->cauchyLorentzDistribution(x);}; //(overridable)


/*
###################
//Sampling
###################
*/ 

void cauchyLorentz::Sampling(){
  FiniteFunction::Sampling();
}
 

// GC: Print info to terminal

void cauchyLorentz::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "peakLoc " << m_peakLoc << std::endl;
  std::cout << "gamma: " << m_gamma << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
  std::cout << std::endl; 
}


/*
###################
// Negative Crystal Ball Distribution
###################
*/ 

//Empty constructor
negativeCrystalBall::negativeCrystalBall(){
  m_RMin = -5.0;
  m_RMax = 5.0;

// GC: Setting member values that best describe the mystery data
  m_mean = 0;
  m_stdev = 1;
  m_n = 1.5; // GC: Must be greater than 1
  m_alpha = 2; // GC: Must be greater than 0

  this->checkPath("NegativeCrystalBallFunction");
  m_Integral = NULL;
}

  negativeCrystalBall::negativeCrystalBall(double range_min, double range_max, std::string outfile, double mean_val, double stdev_val, double n_val, double alpha_val){ // GC: Overloaded constructor
      FiniteFunction(range_min, range_max, outfile);
      // GC: Members unique to the Cauchy Lorentz function
      m_mean = mean_val;
      m_stdev = stdev_val;
      m_n = n_val;
      m_alpha = alpha_val;
}

  negativeCrystalBall::~negativeCrystalBall(){
    Gnuplot gp; //Set up gnuplot object
    this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}


/*
###################
//Setters
###################
*/ 

// GC: Defining new setters for the members unique to the negative Crystal Ball distribution
void negativeCrystalBall::setMean(double mean) {m_mean = mean;};
void negativeCrystalBall::setStdev(double stdev) {m_stdev = stdev;};
void negativeCrystalBall::setN(double n) {m_n = n;};
void negativeCrystalBall::setAlpha(double alpha) {m_alpha = alpha;};



/*
###################
//Getters
###################
*/ 

// GC: Defining new getters for the members unique to the negative Crystal Ball distribution
double negativeCrystalBall::mean_val() {return m_mean;};
double negativeCrystalBall::stdev_val() {return m_stdev;};
double negativeCrystalBall::n_val() {return m_n;};
double negativeCrystalBall::alpha_val() {return m_alpha;};



/*
###################
//Integrating
###################
*/ 

double negativeCrystalBall::integrate(int Ndiv){
    FiniteFunction::integrate(Ndiv);
}

double negativeCrystalBall::integral(int Ndiv) { 
  FiniteFunction::integral(Ndiv); 

}



// Check path

void negativeCrystalBall::checkPath(std::string outfile){
  FiniteFunction::checkPath(outfile);
}


/*
###################
//Plotting
###################
*/ 

// GC: Plot function
void negativeCrystalBall::plotFunction(){
  FiniteFunction::plotFunction();
}


// GC: Plot Data
void negativeCrystalBall::plotData(std::vector<double> &points, int Nbins, bool isdata){
  FiniteFunction::plotData(points, Nbins, isdata);
}

// Generate Plot
void negativeCrystalBall::generatePlot(Gnuplot &gp){
  FiniteFunction::generatePlot(gp);
}

/*
###################
//Function evals
###################
*/ 

// GC: There are 5 functions (A,B,C,D,N) that are required to make up the negative Crystal Ball distribution

double negativeCrystalBall::A(){

  return pow(m_n/abs(m_alpha),m_n) * exp(-pow(abs(m_alpha),2) / 2);
}


double negativeCrystalBall::B(){

  return (m_n / abs(m_alpha)) - abs(m_alpha);
}


double negativeCrystalBall::C(){

  return (m_n / abs(m_alpha)) * (1/(m_n - 1)) * exp(-pow(abs(m_alpha),2) / 2) ;
}


double negativeCrystalBall::D(){

  return (sqrt((pi)/2)) * (1 + erf(abs(m_alpha)/sqrt(2))) ;
}

double negativeCrystalBall::N(){ // GC: Required for normalisation

  return 1 / (m_stdev * (C() * D()));
}


// GC: The form of the distribution is determined by the deviation of points from the mean divided by the standard deviation
double negativeCrystalBall::negativeCrystalBallDistribution(double x) {
  
  if ((x-m_mean)/m_stdev > -m_alpha){

    return N() * exp(-pow(x-m_mean,2) / 2*pow(m_stdev,2));

  }
  else{
    
  return N() * A() * pow((B() - (x-m_mean)/m_stdev),-m_n);

  }
  
  };


double negativeCrystalBall::callFunction(double x) {return this->negativeCrystalBallDistribution(x);}; //(overridable)

/*
###################
//Sampling
###################
*/ 

void negativeCrystalBall::Sampling(){
  FiniteFunction::Sampling();
}

// GC: Printing useful info to the terminal

void negativeCrystalBall::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "mean " << m_mean << std::endl;
  std::cout << "stdev: " << m_stdev << std::endl;
  std::cout << "n: " << m_n << std::endl;
  std::cout << "alpha: " << m_alpha << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
  std::cout << std::endl; 
}


