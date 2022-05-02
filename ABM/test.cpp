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

#define PI 3.14159265359

using namespace std;

// g++ -std=c++17 des.cpp -o des.o
// ./des.o 42 42
random_device rd;
mt19937 generator(rd());
mt19937 len_generator(rd());
mt19937 ang_generator(rd());

uniform_real_distribution<double> len_distribution(-0.025, 0.025);
uniform_real_distribution<double> ang_distribution(-0.005 * PI , 0.005 * PI);
uniform_real_distribution<double> uni_length_distribution(-0.0, 0.0);

int main (int argc, char* argv[]) {
  cout << ang_distribution(ang_generator) << '\n';
  cout << ang_distribution(ang_generator) << '\n';
  cout << ang_distribution(ang_generator) << '\n';
}
