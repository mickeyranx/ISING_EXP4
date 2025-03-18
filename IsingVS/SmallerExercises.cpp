#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <set>
using namespace std;


class SmallerExercises {
public:
    static vector<int> Prng(int seed, int iterations);
    static double approxPi(int N);

};

vector<int> SmallerExercises::Prng(int seed, int iterations) {
    int seed_length = to_string(seed).length();

    vector<int> output;

    //string number = to_string(seed);
    int squared = pow(seed, 2);
    string number = to_string(squared);

    for (int i = 0; i < iterations; i++)
    {



        int current_number_length = number.length();
        //add a zero to the number if its a three digit number
        if (current_number_length < 2 * seed_length) {
            int diff = 2 * seed_length - current_number_length;
            for (int j = 0; j < diff; j++)
            {
                number = "1" + number;
            }

        }

        //remove first and last character
        string middle = number.substr(seed_length / 2, seed_length);

        //stoi converts string to int
        int new_number = pow(stoi(middle), 2);
        number = to_string(new_number);
        output.push_back(stoi(middle));
    }

    //transformiere Zahl zwischen 0 und 1

    return output;
}


//approximation of pi with MC-Algorithm (A1 a)
double SmallerExercises::approxPi(int N) {
    //number of generated numbers in sqaure and in circle
    int square = 0;
    int circle = 0;



    for (int i = 0; i < N; i++)
    {

        double a = double(rand() % 1000) / 1000;
        double b = double(rand() % 1000) / 1000;

        double d = pow(a, 2) + pow(b, 2);
        if (d <= 1) {
            circle++;
        }
        square++;
    }

    double pi = double(4 * circle) / square;



    return pi;
}



