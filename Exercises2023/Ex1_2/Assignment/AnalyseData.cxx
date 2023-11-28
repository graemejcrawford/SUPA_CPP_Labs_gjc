// Author: Graeme Crawford
// Date: 26.11.23

// This file is the running script




//Including all relevant libraries
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <vector>
#include <sstream>
#include <cmath>
#include "CustomFunctions.h"



using namespace std; 




int main(){

    
    vector<double> xData,yData; // Define vectors to hold the x and y data
    vector<double> xErrors, yErrors; // Define vectors to hold the x and y errors
    string dataFile = "/workspaces/SUPA_CPP_Labs_gjc/Exercises2023/Ex1_2/input2D_float.txt"; // Define full file path for data
    string errorFile = "/workspaces/SUPA_CPP_Labs_gjc/Exercises2023/Ex1_2/error2D_float.txt"; // Define full file path for errors
    Read(dataFile, xData, yData); // Execute function to read in data and store in vectors
    Read(errorFile, xErrors, yErrors); // Execute function to read in errors and store in vectors

// Making a switch case to execute the function the user requires

int i; // switch works by users inserting an interger for the function they want to execute
bool go = true;

while(go){

cout << "\n Enter 1 to print the data, 2 for magnitudes, 3 for chi square, 4 for x raised to the power of y and any other interger to exit \n" << endl;
cin >> i;
switch(i){

case 1:{

// Integer 1 runs the function to print the x and y values

cout << "You have chosen to print the data \n"<< endl;
print(xData,yData);
    break;

}

case 2:{

    // Integer 2 runs the function to calculate vector magnitudes

    cout << "You have chosen magnitudes \n"<< endl;
    vector<double> vMag; // defining a vector to store all the magnitudes
    for (unsigned i = 0; i < xData.size(); ++i)
        Mag(xData[i], yData[i], vMag); // When Mag function is run all magnitudes are stored in vMag
        print(vMag);
        break;

    

}


 case 3:{

    // Integer 3 runs the function that creates a best fit line (least squares) and calculate chi squared fit value

    cout << "You have chosen chi squared fit \n" << endl;  
    double p,q; // where y = px + q
    double chisq; // Variable to store the value from the chi square test
    vector<double> errMag, yExp; // errMag stores the varaince of the y data, yExp are the expected y values
    for (unsigned i = 0; i < xData.size(); ++i)
        Lsq(xData, yData, errMag, yErrors[i], yExp, chisq, p, q);
         print(p,q,chisq); // printing line of best fit and chi square
         printBestFit(xData, yExp); // printing line of best fit values for plotting
         break;

    }   

case 4:{

    // Integer 4 runs the function that finds x^y, without loops or pow function

    cout << "You have chosen x raised to the power of y \n" << endl;
    vector<double> xY; // Vector to store x^y
    for (unsigned i = 0; i < xData.size(); ++i)
        Power(xData[i], yData[i], xY); // When Power function is run all values are stored in xY
        printxy(xY);
        break;

    }   

default:{

    // Used to exit the switch case when an invalid interger is given as input

    cout << "A valid integer was not selected, so no function was run." << endl;
    cout << "Exit" << endl;
    go = false;
    break;
    }

    }

}

    return 0;
}  
    


    



