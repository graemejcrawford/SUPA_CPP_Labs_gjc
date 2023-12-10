// Author: Graeme Crawford
// Date: 9.12.2023

// Plots the mystery data, the default function, the test functions and performs Metropolis sampling.
// 

// To-do: Command line arguments


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <vector>
#include <sstream>
#include <cmath>
#include "FiniteFunctions1012.h"
#include <random>

using namespace std;

int main(){

   

    std::string filename = ("/workspaces/SUPA_CPP_Labs_gjc/Exercises2023/Ex3_4/Outputs/data/MysteryData13126.txt");
    vector<double> plotValues;
    Read(filename, plotValues);


 // Create an instance for the FiniteFunction class

    string output = "/workspaces/SUPA_CPP_Labs_gjc/Exercises2023/Ex3_4/Outputs/test";
    FiniteFunction myFunction; // declared instance of the class
  
    myFunction.setRangeMin(-5.0);
    myFunction.setRangeMax(5.0);
    myFunction.setOutfile(output);
    myFunction.plotFunction();
    myFunction.plotData(plotValues, 100, true);
    myFunction.Sampling();
    myFunction.printInfo();



    normalDistribution myNormalDistribution;

    myNormalDistribution.setRangeMin(-5.0);
    myNormalDistribution.setRangeMax(5.0);
    myNormalDistribution.plotFunction();
    myNormalDistribution.plotData(plotValues, 100, true);
    myNormalDistribution.Sampling();
    myNormalDistribution.printInfo();


    cauchyLorentz myCauchyLorentz;

    myCauchyLorentz.setRangeMin(-5);
    myCauchyLorentz.setRangeMax(5.0);
    myCauchyLorentz.plotFunction();
    myCauchyLorentz.plotData(plotValues, 100, true);
    myCauchyLorentz.Sampling();
    myCauchyLorentz.printInfo();
  
    
    
  

    negativeCrystalBall myNegativeCrystalBall;

    myNegativeCrystalBall.setRangeMin(-5.0);
    myNegativeCrystalBall.setRangeMax(5.0);
    myNegativeCrystalBall.plotFunction();
    myNegativeCrystalBall.plotData(plotValues, 100, true);
    myNegativeCrystalBall.Sampling();
    myNegativeCrystalBall.printInfo();


/*

 int random_num = 1;
    int a = -5;
    int b = 5;
    vector<double> randomXValue;
    vector<double> randomYValue;
   
    random_device rd;
    mt19937 mtEngine{rd()};
    uniform_real_distribution<float> rndNumber{a,b};
    for  (int i=0; i < random_num; i++){

        int randomX = rndNumber(mtEngine);
        cout << "Random x value is: " << randomX << endl;
        randomXValue.push_back(randomX);

         float mean = randomX;
         float width = 1; // arbitrarily chose standard deviation
         normal_distribution<float> gaussianpdf{mean,width}; // change these variable names
         int randomY = gaussianpdf(mtEngine);
           cout << "Random y value is: " << randomY << endl;
           randomYValue.push_back(randomY);

}
*/


// define the function for alpha in custom functions
// want to be able to set the value for f(y) and f(xi) such that f(randomXValue) and f(randomYValue) in main




return 0;

}