// Graeme Crawford; 15.11.2023

// cxx file for executing the magnitude of a vector based on user input

#include <iostream> // using system libraries
#include <cmath>
#include "vecMagFnc.h"



int main(){
    

    std::cout << "Please enter a value for the x component of the vector: "; // prompt for user input
    double x_in; // vector components will include decimals
    std::cin >> x_in; // user input for x component
    std::cout << "Please enter a value for the y component of the vector: "; // prompt for user input
    double y_in; // vector components will include decimals
    std::cin >> y_in;  // user input for y component

    double magnitude = vectorMag(x_in,y_in);

    std::cout << "The magnitude of this 2D vector is " << magnitude << std::endl; 

    return 0;

}
