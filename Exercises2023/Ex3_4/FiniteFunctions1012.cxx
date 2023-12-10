// Author: Graeme Crawford
// Date: 9/12/23
// Finite Functions for the plotting of test functions and as well as executing Metropolis sampling


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <numbers>
#include "FiniteFunctions1012.h"
#include <filesystem> //To check extensions in a nice way
#include <random>
#include <vector>

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 

using std::filesystem::path;


// Function to read in the data

void Read(std::string file, std::vector<double>& Values){

    std::ifstream mysteryData(file); // Using ifstream to read in the inputs
   
    //mysteryData.ignore(numeric_limits<streamsize>::max(), '\n');
    if (mysteryData.fail()) { // if statement tells us if the file was not opened successfully
        std::cout << "File did not open" << std::endl;
        exit(1);

    }

    double Data;
   
    
    while(!mysteryData.eof() ) { 
        mysteryData>> Data;
        Values.push_back(Data);
       

        }

   sort(Values.begin(), Values.end());
   mysteryData.close(); // Closed after finishing

    }




//Empty constructor
FiniteFunction::FiniteFunction(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
  m_sampleData.empty();
}

//initialised constructor
FiniteFunction::FiniteFunction(double range_min, double range_max, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
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
double FiniteFunction::integrate(int Ndiv){ //private
 double sum = 0;
  double widths =  (this->rangeMax() - this->rangeMin()) / Ndiv ;
  for (int i = 0; i<Ndiv; i++) // we don't want last Ndiv as the last rectangle starts before this point
  {
    double xEvals = this->rangeMin() + (i*widths);
    // Now use callFunction to get the value of the rectangle at each step
    double heights = this->callFunction(xEvals);
    double rectArea = heights * widths;

    sum += rectArea;
    
  }
  std::cout << "widths are " << widths << std::endl;
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



void FiniteFunction::Sampling(){

   
    int random_num = 1;
    std::random_device rd;
    std::mt19937 mtEngine{rd()};
    std::uniform_real_distribution<float> rndNumber{m_RMin,m_RMax};
    for (int j=0; j < random_num; j++){
    double rndX = rndNumber(mtEngine);
    m_randomX[j] = rndX;

    }
     
    for (int i=0; i < 1000; i++){

    float m_mean_sampling = m_randomX[i];
    //std::cout << "Random mean:  " << m_mean_sampling << std::endl;
    float m_width = 1.0; // arbitrarily chose standard deviation
    std::normal_distribution<float> gaussianpdf{m_mean_sampling,m_width}; // change these variable names
    double rndY = gaussianpdf(mtEngine);
    m_randomY[i] = rndY;

    int m_Tmin = 0;
    int m_Tmax = 1;
    std::uniform_real_distribution<float> rndTNumber{m_Tmin,m_Tmax};
    double rndT = rndTNumber(mtEngine);
    m_randomT[i] = rndT;

    double A;

   if (callFunction(m_randomY[i])/callFunction(m_randomX[i]) < 1){

        //std::cout << "A is :  " << callFunction(m_randomY[i])/callFunction(m_randomX[i]) << std::endl;

        A = callFunction(m_randomY[i])/callFunction(m_randomX[i]);
  
    }

    else { 

        A = 1;
       
    }

    if (m_randomT[i] < A){

      m_randomX[i+1] = m_randomY[i]; // this should then set the next x component (i+1) to y in the next iteration of the loop
          
        }

    else{

      m_randomX[i+1] = m_randomX[i];        

    }

        m_sampleData.push_back(m_randomX[i+1]);
        sort(m_sampleData.begin(), m_sampleData.end());
        this->plotData(m_sampleData,75,false);
   

    }    
 }

//Print (overridable)
void FiniteFunction::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
  std::cout << "Function at x=0 is: " << callFunction(0) << std::endl;
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
    if (bindex<0 || bindex>Nbins){
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





// Normal distribution

  //Empty constructor
normalDistribution::normalDistribution(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  m_mean = 0;
  m_stdev = 0.5;
  this->checkPath("NormalFunction");
  m_Integral = NULL;
}

  normalDistribution::normalDistribution(double range_min, double range_max, std::string outfile, double mean_val, double stdev_val){ //Variable constructor
      FiniteFunction(range_min, range_max, outfile); //Variable constructor
      m_mean = mean_val;
      m_stdev = stdev_val;
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
void normalDistribution::setMean(double mean) {m_mean = mean;};
void normalDistribution::setStdev(double stdev) {m_stdev = stdev;};




/*
###################
//Getters
###################
*/ 
double normalDistribution::mean_val() {return m_mean;};
double normalDistribution::stdev_val() {return m_stdev;};



// integrate

double normalDistribution::integrate(int Ndiv){
    FiniteFunction::integrate(Ndiv);
}

double normalDistribution::integral(int Ndiv) { 
  FiniteFunction::integral(Ndiv); 

}

// Plot function

void normalDistribution::plotFunction(){
  FiniteFunction::plotFunction();
}

// Check path

void normalDistribution::checkPath(std::string outfile){
  FiniteFunction::checkPath(outfile);
}

// Generate Plot

void normalDistribution::plotData(std::vector<double> &points, int Nbins, bool isdata){
  FiniteFunction::plotData(points, Nbins, isdata);
}


void normalDistribution::generatePlot(Gnuplot &gp){
  FiniteFunction::generatePlot(gp);
}

// Sampling


void normalDistribution::Sampling(){
  FiniteFunction::Sampling();
}


// Function Evals
double normalDistribution::normDist(double x) {return (1/(m_stdev * sqrt(m_stdev * 3.14))) *  exp(pow((x-m_mean)/2,2)*(-0.5));};
double normalDistribution::callFunction(double x) {return this->normDist(x);}; //(overridable)

// exp(pow(((x-0)/2),2) * (-1/2))

void normalDistribution::printInfo(){

  std::cout << "meanVal: " << m_mean << std::endl;
  std::cout << "stdevVal: " << m_stdev << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
  std::cout << "Function at x=0 is: " << callFunction(0) << std::endl;
   std::cout << std::endl;
 

  
}




// Cauchy-Lorentz


  //Empty constructor
cauchyLorentz::cauchyLorentz(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  m_peakLoc = 0;
  m_gamma = 0.76;
  this->checkPath("CauchyLorentz");
  m_Integral = NULL;
  m_sampleData.empty();
  //m_randomX[0] = NULL;
  //m_randomY[0] = NULL;
}

  cauchyLorentz::cauchyLorentz(double range_min, double range_max, std::string outfile, double peakLoc_val, double gamma_val){ //Variable constructor
      FiniteFunction(range_min, range_max, outfile); //Variable constructor
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
void cauchyLorentz::setPeakLoc(double peakLoc) {m_peakLoc = peakLoc;};
void cauchyLorentz::setGamma(double gamma) {m_gamma = gamma;};
//void cauchyLorentz::setRandomX(int randomX) {m_randomX = randomX;};
//void cauchyLorentz::setRandomY(int randomY) {m_randomY = randomY;};
//void cauchyLorentz::setRandomT(int randomY) {m_randomT = randomT;};




/*
###################
//Getters
###################
*/ 
double cauchyLorentz::peakLoc_val() {return m_peakLoc;};
double cauchyLorentz::gamma_val() {return m_gamma;};

//double cauchyLorentz::randomT_val() {return m_randomT;};



// integrate

double cauchyLorentz::integrate(int Ndiv){
    FiniteFunction::integrate(Ndiv);
}

double cauchyLorentz::integral(int Ndiv) { 
  FiniteFunction::integral(Ndiv); 

}


// Function Evals
double cauchyLorentz::cauchyLorentzDistribution(double x) {return (1/(3.14 * m_gamma *(1+ pow((x-m_peakLoc)/m_gamma,2))));};
double cauchyLorentz::callFunction(double x) {return this->cauchyLorentzDistribution(x);}; //(overridable)



// Plot function

void cauchyLorentz::plotFunction(){
  FiniteFunction::plotFunction();
}

// Check path

void cauchyLorentz::checkPath(std::string outfile){
  FiniteFunction::checkPath(outfile);
}

// Generate Plot

void cauchyLorentz::plotData(std::vector<double> &points, int Nbins, bool isdata){
  FiniteFunction::plotData(points, Nbins, isdata);
}


void cauchyLorentz::generatePlot(Gnuplot &gp){
  FiniteFunction::generatePlot(gp);
}

// Sampling


void cauchyLorentz::Sampling(){
  FiniteFunction::Sampling();
}
 

void cauchyLorentz::printInfo(){

  std::cout << "peakLoc: " << m_peakLoc << std::endl;
  std::cout << "gamma: " << m_gamma << std::endl;
  std::cout << "Minimum Value: " << m_RMin << std::endl;
  std::cout << "Maximum Value: " << m_RMax << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
  std::cout << "Function at x=0 is: " << callFunction(0) << std::endl;
  std::cout << std::endl;

 
 

}


// Negative Crystal Ball Distribution

// Normal distribution

  //Empty constructor
negativeCrystalBall::negativeCrystalBall(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  m_mean = 0;
  m_stdev = 1;
  m_n = 1.5;
  m_alpha = 2;
  this->checkPath("NegativeCrystalBallFunction");
  m_Integral = NULL;
}

  negativeCrystalBall::negativeCrystalBall(double range_min, double range_max, std::string outfile, double mean_val, double stdev_val, double n_val, double alpha_val){ //Variable constructor
      FiniteFunction(range_min, range_max, outfile); //Variable constructor
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
void negativeCrystalBall::setMean(double mean) {m_mean = mean;};
void negativeCrystalBall::setStdev(double stdev) {m_stdev = stdev;};
void negativeCrystalBall::setN(double n) {m_n = n;};
void negativeCrystalBall::setAlpha(double alpha) {m_alpha = alpha;};



/*
###################
//Getters
###################
*/ 
double negativeCrystalBall::mean_val() {return m_mean;};
double negativeCrystalBall::stdev_val() {return m_stdev;};
double negativeCrystalBall::n_val() {return m_n;};
double negativeCrystalBall::alpha_val() {return m_alpha;};



// integrate

double negativeCrystalBall::integrate(int Ndiv){
    FiniteFunction::integrate(Ndiv);
}

double negativeCrystalBall::integral(int Ndiv) { 
  FiniteFunction::integral(Ndiv); 

}

// Plot function

void negativeCrystalBall::plotFunction(){
  FiniteFunction::plotFunction();
}

// Check path

void negativeCrystalBall::checkPath(std::string outfile){
  FiniteFunction::checkPath(outfile);
}

// Generate Plot

void negativeCrystalBall::plotData(std::vector<double> &points, int Nbins, bool isdata){
  FiniteFunction::plotData(points, Nbins, isdata);
}


void negativeCrystalBall::generatePlot(Gnuplot &gp){
  FiniteFunction::generatePlot(gp);
}

double negativeCrystalBall::A(){

  return pow(m_n/abs(m_alpha),m_n) * exp(-pow(abs(m_alpha),2) / 2);
}


double negativeCrystalBall::B(){

  return (m_n / abs(m_alpha)) - abs(m_alpha);
}


double negativeCrystalBall::N(){

  return 1 / (m_stdev * (C() * D()));
}

double negativeCrystalBall::C(){

  return (m_n / abs(m_alpha)) * (1/(m_n - 1)) * exp(-pow(abs(m_alpha),2) / 2) ;
}


double negativeCrystalBall::D(){

  return (sqrt((3.14)/2)) * (1 + erf(abs(m_alpha)/sqrt(2))) ;
}

// Function Evals
double negativeCrystalBall::negativeCrystalBallDistribution(double x) {
  
  if ((x-m_mean)/m_stdev > -m_alpha){

    return N() * exp(-pow(x-m_mean,2) / 2*pow(m_stdev,2));

  }
  else{
    
  return N() * A() * pow((B() - (x-m_mean)/m_stdev),-m_n);

  }
  
  };


double negativeCrystalBall::callFunction(double x) {return this->negativeCrystalBallDistribution(x);}; //(overridable)

// Sampling


void negativeCrystalBall::Sampling(){
  FiniteFunction::Sampling();
}

void negativeCrystalBall::printInfo(){

  std::cout << "meanVal: " << m_mean << std::endl;
  std::cout << "stdevVal: " << m_stdev << std::endl;
  std::cout << "nVal: " << m_n << std::endl;
  std::cout << "alphaVal: " << m_alpha << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
  std::cout << "Function at x=0 is:" << callFunction(0) << std::endl;
 

  
}