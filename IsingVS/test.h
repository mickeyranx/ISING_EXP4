#pragma once
#include <random>
#include <iostream>
using namespace std;

class test
{
public:

    static default_random_engine generator;
    static uniform_real_distribution<double> distr;

    static void testrng();


};

