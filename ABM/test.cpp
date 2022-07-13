#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* exp */
#include <time.h>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <list>
#include <vector>
#include <filesystem>
#include "sampled_distribution.hpp"

#define PI 3.14159265359

using namespace std;

// g++ -std=c++17 des.cpp -o des.o
// g++ -std=c++17 -stdlib=libc++ test.cpp -o test.o; ./test.o; python3 test.py
// ./des.o 42 42

void distribution()
{
    auto logFunc = [](double x) {
        const double k = -1.0;
        const double m = 10;
        return x*(k*std::log(x) + m - k); // PDF(x) = k*log(x) + m
    };
    auto sinFunc = [](double x) {
        const double x0 = 0.0;
        const double L = 44.0;
        return x + 0.75 * L * std::sin( 2. * PI * ( x - x0) / L  - PI / 4.0 ) / ( 2.0 * PI );
      }; // PDF(x) = 1 + 0.75cos(2 pi (x-x0) / K - pi / 4 )

    std::mt19937 gen;
    //Sampled_distribution<> dist(logFunc, 1.0, 1e4);
    Sampled_distribution<> dist(sinFunc, 0.0, 44.);
    std::ofstream file("d.txt");
    for (int i = 0; i < 100000; i++) {
      file << dist(gen) << std::endl;
    }


void matrix_inverse()
{

}


int main()
{

}
