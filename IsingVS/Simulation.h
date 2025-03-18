#pragma once
#include <vector>
#include <random>
using namespace std;
class Simulation {
public:
    //Lattice length
    static int L;

    //old configuration
    static vector<int> config_1;

    //new configuration
    static vector<int> config_2;

    static default_random_engine generator;
    static uniform_real_distribution<double> distr;


    //
    static vector<int> initializeLatticeCold(int L);

    static vector<int> initializeLatticeHot(int L);

    //return a small list of 4 indices of Positions of the neighbors (top, right, bottom, left)
    static vector<int> getNeighborPos(int i, int L);

    //calculate the change in energy with new spin and next 2D neighbors
    static int changeInEnergy(int pos, int s);

    static int changeInEnergy(int pos, int s_new, double h);

    
    static void sweepMetropolis(double beta, double h);

    static void sweepMetropolisMultihit(double beta, double h, int tries);

    static void sweepHeatbath(double beta, double h);

    static void init(int L, int therm_size, double beta, double h);

    static tuple<double, double> draw(double beta, double h, int sample_size_interval);

    static double averageEnergy();

    static double averageEnergy(int K, int L, double h, vector<int>& config);

    

    static double averageMagnetisation();

    static double averageMagnetisation(int M, vector<int>& config);

    

    static void printConfig(vector<int>& config);

};

