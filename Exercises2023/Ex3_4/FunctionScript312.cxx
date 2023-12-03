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
#include "FiniteFunctions312.h"

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
    myFunction.printInfo();
  


    normalDistribution myNormalDistribution;

    myNormalDistribution.setRangeMin(-5.0);
    myNormalDistribution.setRangeMax(5.0);
    myNormalDistribution.plotFunction();
    myNormalDistribution.plotData(plotValues, 100, true);
    myNormalDistribution.printInfo();


return 0;

}