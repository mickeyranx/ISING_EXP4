#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <set>
#include <tuple>
#include <iomanip>
#include "Simulation.h"
#include <random>
#include "test.h"
#include <string>
#include <fstream>
using namespace std;



// RNG with the box müller method
static double boxMueller(double mu, double sigma) {
    //Zufallszahlen aus GLeichverteilung
    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_real_distribution<double> unif_distr(0.0, 1.0);
    double m = unif_distr(gen);
    double n = unif_distr(gen);
    double r = sigma * sqrt(-2.0 * log(m));

    return r * cos(2.0 * M_PI * n) + mu;

}

// approximation of a certain gaussian Integral with MC Integration (A1 b)
static double approxGaussianIntegral(int sample_size) {
    double upper = 1.5;
    double lower = -1.5;
    double sigma = 1;
    double sum = 0;
    set<double> x_values;
    //calculate with box-müller
    for (int i = 0; i < sample_size; i++)
    {
        double x = boxMueller(0, sigma);
        while (x < lower || x > upper) { //only allow samples within boundaries
            x = boxMueller(0, sigma);
        }
        x_values.insert(x);
    }
    //cout << x_values.size() << endl;
    for (double x : x_values)
    {
        double n = exp(-pow(x, 2) / 2);
        sum += n;
    }

    //approximate with integration and summs
    return (upper - lower) * sum/sample_size;
    
}

static void approximatonOfPi(int sample_size) {
    random_device rand_dev; 
    mt19937 gen(rand_dev());
    uniform_real_distribution<double> unif_distr(-1.0, 1.0);
    int square = 0;
    int circle = 0;
    for (int i = 0; i < sample_size; i++)
    {
        double x = unif_distr(gen);
        double y = unif_distr(gen);
        if (pow(x,2) + pow(y,2)  <= 1) {
            circle++;
            
        }
        square++;

    }

    cout << 4 * circle / double(square) << endl;

}

//return a small list of 4 indices of Positions of the neighbors (top, right, bottom, left)

int main()
{

    /*TODO
    -implement external magnetic field coupling h 
    -implement write to file
    -read file in python and evaluate
    -start making measurements with changing params
    */




    approximatonOfPi(10000);
    //cout << approxGaussianIntegral(40000) << endl;


    clock_t start = clock();
    
    //controll variables
    double beta = 0.6;
    double h = 0;
    int lattice_length = 100;
    int therm_steps = 0;
    int draw_interval = 2;
    int number_of_draws = 400;
    //name of file to write data
    string filename = "test.txt";

    //create file
    ofstream data_file(filename);

    /*

    //to change the update algorithm go to init() and draw() in Simulation.cpp and switch via commenting
    // J is always set to 1

    //thermalize  
    Simulation::init(lattice_length,therm_steps, beta, h);
    //
    for (int i = 0; i < number_of_draws; i++)
    {
        tuple<double, double> mean_vals = Simulation::draw(beta, h, draw_interval);
        //cout << "E = " << fixed << setprecision(7) << get<0>(mean_vals) << ", M = " << get<1>(mean_vals) << endl;


        //write to data file with format E E_err kap kap_err M M_ERR
        data_file << fixed << setprecision(7) << get<1>(mean_vals) << "\n";



    }

    data_file.close();
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    printf("execution time: %.3f sec", elapsed);
    */
}

//sources
//measure execution time (code from https://levelup.gitconnected.com/8-ways-to-measure-execution-time-in-c-c-48634458d0f9#:~:text=The%20function%20clock()%20returns,has%20been%20running%2C%20in%20seconds)