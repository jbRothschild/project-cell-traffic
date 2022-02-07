#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* exp */
#include <time.h>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <random>
#include <list>
#include <vector>

using namespace std;

#define PI 3.14159265359

random_device rd;
mt19937 generator(rd());
uniform_real_distribution<double> len_distribution(-0.025, 0.025);
uniform_real_distribution<double> ang_distribution(-0.02 * PI , 0.02 * PI);

class ABMagent
{
  // Our main actor in the simulation is the agent-> Code any behaviours
  // that are independent of the environment here. Examples include
  // searching, resting, or moving.
  //
  // For simplicity and the smoothness of our graphical output our agents
  // will be coded to move in the world coordinates. These are doubleing
  // point numbers between -1.0 and 1.0, with { 0 , 0 } at the centre.
  //
  protected:
  public:


};

int main (int argc, char **argv) {
  int a = 5;
  int b = a;
  a = 4;
  cout << a << b;
  return 0;
}
