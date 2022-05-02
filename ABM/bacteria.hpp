#ifndef BACTERIA_H
#define BACTERIA_H

class EColi
{
  public:
    double radius = 0.5;
    double max_length = 4.0;
    double inertia = 5.0;
    double growth_rate = 0.013;
    double length = 2.0;
};

class CBacillus
{
  public:
    double radius = 0.5;
    double max_length = 6.0;
    double inertia = 5.0;
    double growth_rate = 0.0065;
    double length = 3.0;
};

class EColiAlt
{
  public:
    double radius = 0.5;
    double max_length = 2.5;
    double inertia = 5.0;
    double growth_rate = 0.0065;
    double length = 1.25;
};

#endif
