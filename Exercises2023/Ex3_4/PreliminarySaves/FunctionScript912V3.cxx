// Author: Graeme Crawford
// Date: 1.12.2023

// Script to test the default code of FinteFunctions
// This script replaces the Test scripts in the Makefile

// To-do: Command line arguments


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <vector>
#include <sstream>
#include <cmath>
#include "FiniteFunctions912V3.h"
#include <random>

using namespace std;

int main(){

   

    string filename = ("/workspaces/SUPA_CPP_Labs_gjc/Exercises2023/Ex3_4/Outputs/data/MysteryData13126.txt");

    ifstream mysteryData(filename); // Using ifstream to read in the inputs
   
    //mysteryData.ignore(numeric_limits<streamsize>::max(), '\n');
    if (mysteryData.fail()) { // if statement tells us if the file was not opened successfully
        cout << "File did not open" << endl;
        exit(1);
    }

  
   
    double Data;
    vector<double> plotValues; 
   
    
    while(!mysteryData.eof() ) { 
        mysteryData>> Data;
        plotValues.push_back(Data);

    }


   sort(plotValues.begin(), plotValues.end());
   cout << plotValues.size() << endl;
   mysteryData.close(); // Closed after finishing
   


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
