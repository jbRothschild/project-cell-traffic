#ifndef BACTERIA_H
#define BACTERIA_H

class EColi
{
  public:
    double radius = 0.92 / 2.0 + 0.05; // 0.65 / 2 + 0.1
    double max_length = 4.56; // 3.71
    double inertia = 5.0;
    double growth_rate = 0.0173;
    double length = 2.25;
};

class EColiA22
{
  public:
    double radius = 1.20 / 2.0 + 0.05; // 0.89 / 2 + 0.1
    double max_length = 4.10; // 3.05
    double inertia = 5.0;
    double growth_rate = 0.0173;
    double length = 1.25;
};

class BSubt
{
  public:
    double radius = 0.83 / 2.0 + 0.1;
    double max_length = 7.95;
    double inertia = 5.0;
    double growth_rate = 0.0039;
    double length = 3.00;
};

#endif
