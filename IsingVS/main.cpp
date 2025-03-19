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
#include <string>
#include <fstream>
#include <numeric>
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
            double H = Simulation::averageEnergy(config, K,L , 0 );
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

//Simulation for exercises 3 and 4
static vector<double> startSimulation(int L, double beta, double h,int therm_steps, int N, int draw_interval) {
    //------------------------------
    //            setup
    //------------------------------
    Simulation::initRNG();

    vector<int> config = Simulation::initializeLatticeCold(L);
    //vector<int> config = Simulation::initializeLatticeHot(L);
    for (int i = 0; i < therm_steps; i++)
    {
        //Simulation::sweepHeatbath(config, L, beta, h);
        Simulation::sweepMetropolisMultihit(config, L, beta, h, 5);
    }

    //Observables
    vector<double> energies1(N, 0);
    vector<double> energies2(N, 0);
    vector<double> absmag1(N, 0);
    vector<double> absmag2(N, 0);
    vector<double> mag1(N, 0);
    vector<double> mag2(N, 0);
    //------------------------------
    //            sweeps
    //------------------------------
    double K = (double) L * L;
    for (int i = 0; i < N; i++)
    {
        Simulation::draw(config, L ,beta , h , draw_interval);
        double E = Simulation::averageEnergy(config, L*L, L, h);
        energies1[i] = E/K;
        energies2[i] = E/K * E/K;
        double M = Simulation::averageMagnetisation(L * L, config);
        double absmag_i = abs(M) / K;
        absmag1[i] = absmag_i;
        absmag2[i] = absmag_i * absmag_i;
        mag1[i] = M / K;
        mag2[i] = M / K * M / K;
    }
    //------------------------------
    //    calculate observables
    //------------------------------
    vector<double> vals = {};
    double e_mean = accumulate(energies1.begin(), energies1.end(), 0.0)/ N;
    double e2_mean = accumulate(energies2.begin(), energies2.end(), 0.0) / N;
    double m_mean = accumulate(mag1.begin(), mag1.end(), 0.0) / N;
    double m2_mean = accumulate(mag2.begin(), mag2.end(), 0.0) / N;
    double absm_mean = accumulate(absmag1.begin(), absmag1.end(), 0.0) / N;
    double absm2_mean = accumulate(absmag2.begin(), absmag2.end(), 0.0) / N;
    vals.push_back(e2_mean);
    vals.push_back(sqrt(abs(pow(e_mean, 2) - e2_mean)) / (N - 1));
    vals.push_back(m_mean);
    vals.push_back(sqrt(abs(pow(m_mean, 2) - m2_mean)) / (N - 1));
    vals.push_back(e2_mean);
    vals.push_back(sqrt(abs(pow(absm_mean, 2) - absm2_mean)) / (N - 1));
    return vals;
   
}


int main()
{
    clock_t start = clock();
    /*TODO
    
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
    //          exercises 3 & 4
    //--------------------------------------------
    
    /*
    vector<double> betas = {};
    for (double beta = 0; beta < 1.05; beta+= 0.05)
    {

    }
    */

    double beta = 0.6;
    double h = 0;
    int therm_steps = 0;
    int draw_interval = 1;
    
    int L = 32;
    //number of draws
    int N = 1000;

    vector<double> vals = startSimulation(L, beta, h, therm_steps, N, draw_interval);
    ofstream File("test.txt");
    File << fixed << setprecision(5);
    File << "beta" << "\t" << "<e>" << "\t" << "de" << "\t" << "<m>" << "\t" << "dm" << "\t" << "<|m|>" << "\t" << "d|m|" << "\n";
    File << beta << "\t" << vals[0] << "\t" << vals[1] << "\t" << vals[2] << "\t" << vals[3] << "\t" << vals[4] << "\t" << vals[5] << "\n";
    File.close();
    
    //controll variables
    
    //name of file to write data
    //string filename = "test.txt";

    //create file
    //ofstream data_file(filename);

    /*

    
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