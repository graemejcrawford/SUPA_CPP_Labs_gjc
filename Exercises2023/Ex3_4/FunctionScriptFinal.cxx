// Author: Graeme Crawford
// Date: 9.12.2023

// Main function which plots the mystery data, the default function, the test functions and sampled data from Metropolis sampling.


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <vector>
#include <sstream>
#include <cmath>
#include "FiniteFunctionsFinal.h"



int main(){

   
    std::string filename = ("/workspaces/SUPA_CPP_Labs_gjc/Exercises2023/Ex3_4/Outputs/data/MysteryData13126.txt"); //GC: Path to mystery data file
    std::vector<double> plotValues; // GC: Will save mystery data as a vector
    Read(filename, plotValues);  // GC: Using read functionality to read in the mystery data
   
    
    FiniteFunction myFunction; // declared instance of the class for the default function
  
    
    myFunction.plotFunction(); // GC: Plots the default function
    myFunction.plotData(plotValues, 100, true); // GC: Using 100 bins to plot mystery data
    myFunction.Sampling(); // GC: Perform sampling on default function and plots sampled data
    myFunction.printInfo(); // GC: Printing useful information about the function

   



    normalDistribution myNormalDistribution; // declared instance of the class for the normal distribution


    myNormalDistribution.plotFunction();
    myNormalDistribution.plotData(plotValues, 100, true);
    myNormalDistribution.Sampling();
    myNormalDistribution.printInfo();
    


    cauchyLorentz myCauchyLorentz; // declared instance of the class for the Cauchy Lorentz distribution


    myCauchyLorentz.plotFunction();
    myCauchyLorentz.plotData(plotValues, 100, true);
    myCauchyLorentz.Sampling();
    myCauchyLorentz.printInfo();
    
  
    

    negativeCrystalBall myNegativeCrystalBall; // declared instance of the class for the negative Crystal Ball distribution


    myNegativeCrystalBall.plotFunction();
    myNegativeCrystalBall.plotData(plotValues, 100, true);
    myNegativeCrystalBall.Sampling();
    myNegativeCrystalBall.printInfo();
   


return 0;

}