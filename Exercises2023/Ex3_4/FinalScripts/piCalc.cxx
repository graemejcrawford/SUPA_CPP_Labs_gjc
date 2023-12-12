#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <numbers>
#include <random>
#include <vector>

// Imagine a circle within a rectangle, with the circle touching the midpoint of each side of the rectangle
// Given the area of quarter of a circle is (π r^2)/4 , the area of a quarter of the rectangle is r^2

/* If we can approximate the number of points within the circle to represent the area of the circle,
 the ratio of the number of points in the circle to the overall number of points (area of rectangle) will be
 π/4*/

// To check if a point lies inside of the circle, we need to calculate the magnitude of the (x,y) vector
// Function to calculate the magnitude 
void Magnitude(double xPosition, double yPosition, std::vector<double> &vecMag)
{

    double magCalc = sqrt(xPosition * xPosition + yPosition * yPosition); // Calculate magnitude
    vecMag.push_back(magCalc); // Storing magnitudes
};

int main()
{
    // The two inputs required are the radius of the circle and the number of random points/darts
    double radius;
    double n_random;

    // Prompt user to enter a value for the radius
    std::cout << "\n Enter a value for the radius \n"<< std::endl;
    std::cin >> radius;
     if (radius <= 0){
        std::cout << "Invlaid value, radius must be greater than 0 \n" << std::endl;
        exit(1);
    }

    // Prompt user to enter a value for the number of random points
    std::cout << std::endl;
    std::cout << "Enter a value for the number of random points \n" << std::endl;
    std::cin >> n_random;
        if (n_random <= 0){
            std::cout << "Invlaid value, number of values must be greater than 0 \n" << std::endl;
        exit(1);
    }
    std::cout << std::endl;

    double dartsCircle = 0; // Defining a variable that will count the number of points inside the circle

    std::vector<double> randomX; // Generate random x coordinates
    std::vector<double> randomY; // Generate random y coordinates
    std::random_device rd;
    std::mt19937 mtEngine{rd()};
    std::uniform_real_distribution<float> rndNumber{0, 1}; // Generate random numbers from a uniform distribution
    for (int j = 0; j < n_random; j++)
    {
        double rndX = rndNumber(mtEngine);
        double rndY = rndNumber(mtEngine);
        randomX.push_back(rndX); // Storing random X values 
        randomY.push_back(rndY); // Storing random Y values 
        
    }
        
    // Now calculate the magnitude to see if the coordinate lies in the circle
    std::vector<double> vMag; // defining a vector to store all the magnitudes
    for (unsigned i = 0; i < randomX.size(); ++i){

        // Using the magnitude function with generated random numbers as arguments
        Magnitude(randomX[i], randomY[i], vMag); 

    if (vMag[i] < 1){

        dartsCircle ++; // If coordinate lies within the circle increase the count by +1
          
        }
        
    }
        
    double pi = (4 * dartsCircle) / (n_random) ; // The ratio of the area of the circle and rectangle is π/4

    std::printf("%.10f\n", pi); // Printing value of π to 10 decimal places

}
