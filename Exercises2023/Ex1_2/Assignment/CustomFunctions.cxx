// Author: Graeme Crawford
// Date: 26.11.23

// This file contains all the custom functions



//Including all relevant libraries
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <vector>
#include <sstream>
#include <cmath>


using namespace std;

// Function to read in the data and using vectors (dynamic data structure) to store the values

void Read(string filename, vector<double>& xValues, vector<double>& yValues)
{

    ifstream inputData(filename); // Using ifstream to read in the inputs
    inputData.ignore(numeric_limits<streamsize>::max(), '\n'); // Ignoring the first line with column headers
    if (inputData.fail()) { // if statement tells us if the file was not opened successfully
        cout << "File did not open" << endl;
        exit(1);
    }

    int number_lines;

    printf("Please insert number of lines to be read in '%s'\n ", filename.c_str() ); // user inputs number of lines
    cin >> number_lines;

    if (number_lines <= 0){
        cout << "Must insert a positive integer" << endl;
        exit(1);
    }
   
    int presentLine = 0; // Used to keep track of what line in the file we are at
    double x, y;
    char c;
    
    while(inputData >> x >> c >> y  && c == ',') { // Character is used to deal with the values being comma separated
        ++ presentLine;
        xValues.push_back(x);
        yValues.push_back(y);
        if (presentLine == number_lines) break; // break once we have reached the required line or end of the file

      
    }
      if (number_lines > xValues.size()) // warning code if user inputs an invalid number of lines
    {
        cout << "Warning: Requested number of lines is larger than file." << endl;
        cout << "Maximum number of lines is: " << presentLine << endl;
        cout << "x,y values of first five lines are: \n" << endl;
         for (unsigned i = 0; i < 5; ++i)
            cout <<  xValues[i] << " " << yValues[i] << endl;
            exit(1);
        
        }
      

    inputData.close(); // Closed after finishing

   
}


// Function that calculates the magnitude of a vector (x,y) 
// Arguments are the x and y components of the vectors and a vector (dynamic data structure) to store the magnitudes

void Mag( double xData, double yData, vector<double>& vecMag) 
{
   
   double s = sqrt(xData*xData + yData*yData);
   vecMag.push_back(s); 
   
};



/* Function that calculates the best fit line using the least squares method
   and then calculates the reduced chi square value */ 


void Lsq( vector<double>& xData, vector<double>& yData, vector<double>& errMag, double yErrors, vector<double>& yExp, double& redChisq, double& p, double& q)
{
    
    int N = 25; // Number of observations, since assuming full dataset is used for the rest of the task


    double sumx = 0, sumx2 = 0, sumy = 0, sumxy = 0; // Variables to store the summations
    for (int i=0; i < N; i++){
    
        sumx += xData[i]; 
        sumy += yData[i];
        sumx2 += pow(xData[i],2);
        sumxy += xData[i]*yData[i];

    }
    
   
    // Calculation of p and q that will be used in the least squares equation

    p = (N*sumxy - sumx*sumy) / (N*sumx2 - sumx*sumx);
    q = (sumx2*sumy - sumxy*sumx) / (N*sumx2 - sumx*sumx);


    // Calculation of the variance of y

     double varY = yErrors*yErrors; 
            errMag.push_back(varY); 



    /* Reduced chi square = chi square / degrees of freedom 
        degrees of freedom = Number of measurements - number of parameters fitted */

    double rChisq = 0; // Rediced chi square variable
    for (int i=0; i < N; i++){

        double yFit = (p*xData[i] + q); // Using the least squares equation to find fitted y values
        yExp.push_back(yFit);
        rChisq += (((pow(yData[i] - yExp[i],2)) / errMag[i])) / (N-2); // reduced chi square formula provided

    }
        
        redChisq = rChisq;


    
    };




// Function to calculate x^y without using the pow function or a for loop

void Power(double xData, double yData, vector<double>& power){ 


yData = yData + 0.5; /* Add 0.5 to the y value as compiler rounds down */
int y_values = (int)yData; // converting doubles to integers


if (power.size() < 25){ // using recursion instead of for loop. Again assuming full dataset in use.

    double xY = exp (y_values * log(xData)); // Using Exp and Log to create own power function
    power.push_back(xY);

}

};
   


// Using function overloading to have single print function


// Print function for data 

void print(vector<double> xData, vector<double>yData){

cout << "x and y data is: \n" << endl;
for (unsigned i = 0; i < xData.size(); ++i){
    cout << xData[i] << "," << yData[i] << endl;

}

}

// Print function for Vector Magnitudes 

void print(vector<double> magnitudes){

  cout << "Vector Magnitudes are: \n" << endl;
  ofstream out("../Outputs/VMag_output.txt"); // use this to create an output .txt file stored in Outputs folder
  out << "Vector Magnitudes are: \n" << endl;
  for(int i = 0; i < magnitudes.size(); i++){
        printf("%.4f\n", magnitudes[i]);
        out.precision(4);  // 4 decimal places was the lowest in the dataset
        out << fixed << magnitudes[i] << endl; // precision and fixed are used to write to 4 decimal places 
        
  }

  out.close();

} 





// Print function for Least Squares and Reduced Chi Square

void print(double p, double q, double chisq){

   cout << "Least squares fit is " << "y = " <<p << "x + " << q<< endl;
   printf("Value for chi squared fit is %.4f\n", chisq); // using printf to print to 4 decimal places
   
   ofstream out("../Outputs/lsq_output.txt"); // creating ofstream to write outputs to a .txt file
        out << "Least squares fit is " << "y = " <<p << "x + " << q<< endl;
        out << "Chi squared value is: " << chisq << endl;
        out.close();
}



// Print function for x^y values, following same format as print function for magnitudes (so defined separate function)

void printxy(vector<double> XY){

    cout << "x to the power of y: \n" << endl;
    ofstream out("../Outputs/XpowerY_output.txt");
    out << "x to the power of y values: \n" << endl;
     for(int i = 0; i < XY.size(); i++){
        printf("%.4f\n", XY[i]);
        out.precision(4);
        out << fixed << XY[i] << endl;
        

}
out.close();

}


/* Print function that produces the values for the best fit line which will be plotted. 
   Updated plotdata.gp to plot the best fit line */

void printBestFit(vector<double> xData, vector<double> yExp){

    
    ofstream out("../Outputs/bestFit_output.txt");
    for (unsigned i = 0; i < xData.size(); ++i){
    out << xData[i] << "," << yExp[i] << endl;
    }

}





