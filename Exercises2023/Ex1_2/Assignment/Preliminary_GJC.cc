// Graeme Crawford, 15.11.2023

// Preliminary exercises for assignment 1 (not marked)

// Bring in the relevant header files. We will need to define vectors in this exercise.

#include <iostream>
#include <cmath>




int main(){
    std::cout << "Hello World!" << std::endl; // printing out an initial 'Hello World!'

     
    double x = 2.3; // defining the x component of the 2D vector
    double y=4.5; // defining the y component of the 2D vector
    double mag = sqrt(pow(x,2) + pow(y,2)); // Taking the square root of the squared components
    

    std::cout << "The magnitude of this 2D vector is " << mag << std::endl; 

    std::cout << "Please enter a value for the x component of the vector: "; // prompt for user input
    double x_in; // vector components will include decimals
    std::cin >> x_in; // user input for x component
    std::cout << "Please enter a value for the y component of the vector: "; // prompt for user input
    double y_in; // vector components will include decimals
    std::cin >> y_in;  // user input for y component

    double vec_mag = sqrt(pow(x_in,2) + pow(y_in,2)); // Calculating the magnitude

    std::cout << "The magnitude of this 2D vector is " << vec_mag << std::endl; 

    return 0;

}


