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
#include <bitset>
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
//approximation of PI witch MC (A1 a)
static double approximatonOfPi(int sample_size) {
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

    return 4 * circle / double(square);

}

//backtrack algorithm to invoke all possible configurations of a LxL lattice 
static void backtrack(int L, int counter, vector<int> config, vector<vector<int>>& list_of_configs) {
    if (counter >= pow(L, 2)) {
        list_of_configs.push_back(config);
        return;
    }
    config.push_back(1);
    counter++;
    backtrack(L, counter, config, list_of_configs);

    counter--;
    config.pop_back();
    config.push_back(0);
    counter++;
    backtrack(L, counter, config, list_of_configs);
    return;

}

//calculates the mean energy pp, mean absolute magnetism pp 
//and mean magnetism pp explicitly with the partition function over all possible configurations
static void explicitIsing(int L) {
    //-----------------------------------
    //               setup
    //-----------------------------------
    //list of all possible configurations
    vector<vector<int>> list_of_configs = {};
    //create all configs with backtrack-algorithm
    backtrack(L, 0, {}, list_of_configs);
    cout << list_of_configs.size() << endl;
    //initialize list of betas
    vector<double> betas = {};
    for (double i = 0; i <= 1; i+= 0.05)
    {
        betas.push_back(i);
    }
    double K = (double) L * L;
    //-----------------------------------
    //      calculate observables
    //-----------------------------------
    ofstream File("explicit_vals_L=" + to_string(L) + ".txt");
    File << "beta" << "\t" << "<e>" << "\t" << "<m>" << "\t" << "<|m|>" << "\n";
    File << fixed << setprecision(7);
    for (double beta : betas) {
        //1.partition-function
        double Z = 0;
        vector<double> energies = {};
        for (vector<int> config : list_of_configs) {
            double H = Simulation::averageEnergy(K,L , 0, config);
            energies.push_back(H);
            Z += exp(-beta * H);
        }

        //2. mean energy pp, mean magnetism pp, mean absolute magnetism pp
        double mean_energy = 0;
        double mean_mag = 0;
        double mean_abs_mag = 0;
        int i = 0;
        for (vector<int> config : list_of_configs) {
            double H_i = energies[i];
            mean_energy += 1 / K * 1 / Z * exp(-beta * H_i) * H_i;
            double M_i = Simulation::averageMagnetisation(K, config);
            mean_mag += 1 / K * 1 / Z * M_i * exp(-beta * H_i);
            mean_abs_mag += 1 / K * 1 / Z * abs(M_i) * exp(-beta * H_i);
            i++;
        }

        File << beta << "\t" << mean_energy << "\t" << mean_mag << "\t" << mean_abs_mag << "\n";


    }
    
    File.close();





}



int main()
{

    clock_t start = clock();
    /*TODO
    -implement external magnetic field coupling h 
    -implement write to file
    -read file in python and evaluate
    -start making measurements with changing params
    */

    //--------------------------------------------
    //                 exercise 1
    //--------------------------------------------
    /*
    printf("approximation of PI: %.3f", approximatonOfPi(10000));
    printf("\napproximation of gaussian integral: %.4f", approxGaussianIntegral(50000));
    */

    //--------------------------------------------
    //                 exercise 2
    //--------------------------------------------
    /*
    explicitIsing(2);
    explicitIsing(3);
    explicitIsing(4);
    */
    //--------------------------------------------
    //                 exercise 3
    //--------------------------------------------




    
    
    //controll variables
    double beta = 0.6;
    double h = 0;
    int lattice_length = 100;
    int therm_steps = 0;
    int draw_interval = 2;
    int number_of_draws = 400;
    //name of file to write data
    //string filename = "test.txt";

    //create file
    //ofstream data_file(filename);

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
    */
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    printf("execution time: %.3f sec", elapsed);
    
}

//sources
//measure execution time (code from https://levelup.gitconnected.com/8-ways-to-measure-execution-time-in-c-c-48634458d0f9#:~:text=The%20function%20clock()%20returns,has%20been%20running%2C%20in%20seconds)