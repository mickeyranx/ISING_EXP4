#include <iostream>
#include <vector>
#include "Simulation.h"
#include <tuple>
using namespace std;





vector<int> Simulation::config_1 = {};
vector<int> Simulation::config_2 = {};
int Simulation::L = 0;
uniform_real_distribution<double> Simulation::distr = uniform_real_distribution<double>(0.0, 1.0);
default_random_engine Simulation::generator = default_random_engine();



vector<int> Simulation::initializeLatticeCold(int L) {
    vector<int> lattice(pow(L, 2), 1);
    return lattice;
}

vector<int> Simulation::initializeLatticeHot(int L)
{
    vector<int> lattice = {};
    //set current time as seed for rand()
    srand(time(0));
    for (int i = 0; i < pow(L, 2); i++)
    {
        //generate random number
        double p = double(rand() % 101);
        if (p < 50) {
            lattice.push_back(-1);
        }
        else
        {
            lattice.push_back(1);
        }

    }
    
    return lattice;
}


vector<int> Simulation::getNeighborPos(int i, int L) {
    //determine which row we are in
    double row = floor(double(i) / L);
    //cout << "row = " << row << endl;
    vector<int> neigbor_pos;


    //check top neighbor
    if (i - L < 0) {
        //add L^2 - (L - i) 
        neigbor_pos.push_back(pow(L, 2) - (L - i));
    }
    else {
        //add i - L
        neigbor_pos.push_back(i - L);
    }

    //check right neighbor
    if (i + 1 > ((row + 1) * L) - 1) {
        //add i - L - 1

        neigbor_pos.push_back(i - L + 1);
    }
    else {
        //add i+1

        neigbor_pos.push_back(i + 1);
    }

    //check bottom neighbor
    if (i + L > pow(L, 2) - 1) {
        //add i - (L*(L-1))

        neigbor_pos.push_back(i - (L * (L - 1)));
    }
    else {
        //add i+L

        neigbor_pos.push_back(i + L);
    }

    //check left neigbor
    if (i - 1 < ((row + 1) - 1) * L) {
        //i + (L - 1)

        neigbor_pos.push_back(i + (L - 1));
    }
    else {
        //add i-1

        neigbor_pos.push_back(i - 1);
    }
    //cout << "successfully calculated neighbors" << neigbor_pos.size() << endl;
    return neigbor_pos;

}


//change in energy in case of spin flip suggestion
int Simulation::changeInEnergy(int pos, int s) {
    vector<int> neighbors = getNeighborPos(pos, L);
    /*for (int index : neighbors) {
        cout << "index = " << index << endl;
    }*/
    int dH = 2 * s * (config_1[neighbors[0]] + config_1[neighbors[1]] + config_1[neighbors[2]] + config_1[neighbors[3]]);
    //cout << "calculated energy change successfully" << endl;
    return double(dH);
}

int Simulation::changeInEnergy(int pos, int s, double h) {
    vector<int> neighbors = getNeighborPos(pos, L);
    
    int dH = 2 * s * (config_1[neighbors[0]] + config_1[neighbors[1]] + config_1[neighbors[2]] + config_1[neighbors[3]] + h);

    return double(dH);
}


//Metropolis with spin flip suggestion
void Simulation::sweepMetropolis(double beta, double h) {
    for (int i = 0; i < config_1.size(); i++)
    {
        //current spin
        int s_i = config_1[i];
        //sugest spin flip 
        int s_i_new = -s_i;
        //calculate change
        int dH = changeInEnergy(i, s_i, h);

        if (dH < 0) {
            config_1[i] = s_i_new;
            
        }
        else {
            //generate random number between 0 and 1
            
            double rnd = distr(generator);
            //double prob = exp(-beta * dH) / (1 + exp(-beta * dH));
            double prob = exp(-beta * dH);

            //cout << "r must be smaller than " << exp(-beta * dH) << endl;
            //cout << "r = " << r << endl;

            double scale = pow(10, 10);
            //cout << "p = " << p << endl;
            //double b = round(p * scale) / scale;
            //double a = round(r * scale) / scale;

            //cout << "comparing " << a << " < " << b << endl;
            if (round(rnd*scale)/scale < round(prob*scale)/scale) {
                config_1[i] = s_i_new;
                //cout << "flip accepted randomly" << endl;
            }
            else {
                continue;
            }


        }

    }
}

void Simulation::sweepMetropolisMultihit(double beta, double h, int tries) {
    for (int i = 0; i < config_1.size(); i++)
    {
        int s_i = config_1[i];
        
               
        for (int j = 0; j < tries; j++)
        {
            double rnd = distr(generator);
            //select new random spin
            int s_i_new = (rnd < 0.5 ? 1 : -1);



            double dH = 0;
            //calculate change in energy
            if (s_i != s_i_new) {
                dH = changeInEnergy(i, s_i, h);
            }
            
            
            if (dH < 0) {
                config_1[i] = s_i_new;
                //cout << "flip accepted immideatly" << endl;
               
            }
            else {
                //generate random number between 0 and 1
                double rnd = distr(generator);
                double prob = exp(-beta * dH);

                //equalize precision to compare doubles
                double scale = pow(10, 10);
                
                if (round(rnd * scale) / scale < round(prob * scale) / scale) {
                    config_1[i] = s_i_new;
                    
                }
                else {
                    continue;
                }


            }




        }


    }


}


void Simulation::sweepHeatbath(double beta, double h) {
    for (int i = 0; i < config_1.size(); i++)
    {
        int s_i = config_1[i];

        vector<int> neighbors = getNeighborPos(i, L);

        double delta = (config_1[neighbors[0]] + config_1[neighbors[1]] + config_1[neighbors[2]] + config_1[neighbors[3]]);
        double k = beta * (delta + h);
        double z = 2 * cosh(k);
      
        //2 values to compare
        double q = exp(k) / z;
        double r = distr(generator);
        //transform r,p to same precision
        double scale = pow(10, 10);
        double prob = round(q * scale) / scale;
        double rnd = round(r * scale) / scale;

        //cout << "r must be smaller than " << prob << endl;
        //cout << "r = " << r << endl;

        (rnd < prob) ? config_1[i] = 1 : config_1[i] = -1;
        
    }

}

void Simulation::init(int lattice_length, int therm_size ,double beta, double h)
{
    L = lattice_length;
    //initialize Lattice and assign to global variables
    config_1 = initializeLatticeHot(L);
    //config_1 = initializeLatticeCold(L);
    //printConfig(config_1);
    //cout << "\n" << "next config -------" << endl;
    //copy the spins from config_1 to config_2
    config_2 = config_1;
    cout << averageEnergy() << endl;

    //generate seed with current time for rand() function
    srand(time(0));

    for (int i = 0; i < therm_size; i++)
    {
        //sweepMetropolis(beta,h);
        sweepHeatbath(beta, h);

        //printConfig(config_2);
        //copies the vector (not the reference)
        config_1 = config_2;
        
       // cout << "\n" << "next config -------" << endl;
    }



}

tuple<double, double> Simulation::draw(double beta, double h, int draw_interval) {
    for (int i = 0; i < draw_interval; i++)
    {
        //sweepMetropolis(beta, h);
        sweepMetropolisMultihit(beta, h, 1);
        //sweepHeatbath(beta, h);
        //config_1 = config_2;
    }
    //calculate Observables
    double E = averageEnergy();
    double M = averageMagnetisation();
    //return them 
    return make_tuple(E, M);
}

double Simulation::averageEnergy() {
    double sum = 0;
    for (int i = 0; i < pow(L,2); i++)
    {
        int s = config_1[i];
        vector<int> neigbors = getNeighborPos(i,L);
        sum += -double(s)*(config_1[neigbors[1]] + config_1[neigbors[2]]);
        
        
    }
    return 1/(pow(L,2)) * sum;
}

double Simulation::averageMagnetisation() {
    double sum = 0;
    for (int i = 0; i < pow(L, 2); i++)
    {
        sum += config_1[i];

    }
    return 1 / (pow(L, 2)) * abs(sum);

}


//for debugging
void Simulation::printConfig(vector<int> &config) {
    int row = 0;
    int L_square = config.size();
    for (int i = 0; i < L_square; i++)
    {
        int current_spin = config[i];
        if (floor(double(i) / sqrt(L_square)) > row) {
            row++;
            cout << "\n";
        }

        if(current_spin == -1){
            cout << current_spin << " ";
        }
        else {
            cout << " " << current_spin << " ";
        }

            
   
    }

}








