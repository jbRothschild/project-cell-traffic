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
#define SPEED 0.005
#define NAGENTS 20000
#define XCHANNEL 45.0
#define YCHANNEL 12.0
#define MAXLENCELL 4.0

using namespace std;

// g++ -std=c++17 des.cpp -o des.o
// ./des.o 42 42
random_device rd;
mt19937 generator(rd());
mt19937 len_generator(rd());
mt19937 ang_generator(rd());
uniform_real_distribution<double> len_distribution(-0.025, 0.025);
uniform_real_distribution<double> ang_distribution(-0.005 * PI , 0.005 * PI);
uniform_real_distribution<double> uni_length_distribution(-0.05, 0.05);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void pnt2line (
                double A_x,
                double A_y,
                double B_x,
                double B_y,
                double pnt_x,
                double pnt_y,
                double &pnt_on_line_x,
                double &pnt_on_line_y,
                double &distance
              )
{

  double AB_x = A_x - B_x;
  double AB_y = A_y - B_y;
  double line_length = sqrt( pow(AB_x, 2.0) + pow(AB_y, 2.0) );
  double Apnt_x = pnt_x - B_x;
  double Apnt_y = pnt_y - B_y;
  double t = ( AB_x * Apnt_x + AB_y * Apnt_y ) / pow(line_length, 2.0);
  if(t < 0)
  {
    t = 0;
  }
  if(t > 1)
  {
    t = 1;
  }

  pnt_on_line_x = B_x + t * AB_x;
  pnt_on_line_y = B_y + t * AB_y;
  distance = sqrt( pow(pnt_x - pnt_on_line_x, 2.0)
                  + pow(pnt_y - pnt_on_line_y, 2.0) );
  /*
  cout << "\n A_x : " << A_x;
  cout << "\n B_x : " << B_x;
  cout << "\n pnt_x : " << pnt_x;
  cout << "\n A_y : " << A_y;
  cout << "\n B_y : " << B_y;
  cout << "\n pnt_y : " << pnt_y;
  cout << "\n t : " << t;
  cout << "\n pnt_on_line_x : " << pnt_on_line_x;
  cout << "\n pnt_on_line_y : " << pnt_on_line_y;
  cout << "\n distance : " << distance;
  */

};


class Environment;

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
  friend class environment;
  protected:
  public:
    Environment* environment;
    string label;
    double radius;
    double length;
    double growth_rate;
    double max_length;
    double split_length;
    double inertia;
    double x;
    double y;
    double angle;
    double vel_x;
    double vel_y;
    double acc_x;
    double acc_y;
    double acc_angle;
    double vel_angle;
    double force_x_prev = 0.0;
    double force_y_prev = 0.0;
    double torque_prev = 0.0;
    double force_x = 0.0;
    double force_y = 0.0;
    double torque = 0.0;
    ABMagent* prev_ = NULL;
    ABMagent* next_ = NULL;
    ABMagent
    (
      Environment* environment_,
      double x_,
      double y_,
      double radius_,
      double max_length_,
      double length_,
      double vel_x_,
      double vel_y_,
      double acc_x_,
      double acc_y_,
      double angle_,
      double inertia_,
      double vel_angle_,
      double growth_rate_,
      string label_,
      double force_x_prev_,
      double force_y_prev_,
      double torque_prev_
    );
    void addForce
    (
      double point_x,
      double point_y,
      double ext_force_x,
      double ext_force_y
    );
    void split();
    void print();
    void grow(double dt);
    void move(double dt, double damping_lin, double damping_tor);
    void move_old(double dt, double damping_lin, double damping_tor);

};


class Environment
{
  public:
    double dt;
    double save_time;
    /* full force
    double mass = 1.0 / 5.0;
    double beta = 0.8;
    double g = 1000.0 / 3.6;
    double k_n = 2000000.0; // 3333
    double gamma_n = 220.0; // 0.00
    double gamma_t = 220.0; // 0.00
    double damping_lin = 0.8 * k_n; // 2000
    double damping_tor = 0.8 * k_n; // 8000 * PI
    double mu_cc = 0.1;
    double mu_cw = 0.8;
    */
    double mass = 1.0 / 5.0;
    double beta = 0.8;
    double g = 1000.0 / 3.6;
    double k_n = 4000000.0; // 3333
    double gamma_n = 0.0; // 0.00
    double gamma_t = 0.0; // 0.00
    double damping_lin = 200.0 * 60; // 12000
    double damping_tor = 200.0 * 60; // 8000 * PI
    double mu_cc = 0.1;
    double mu_cw = 0.8;
    bool split_bool= false;
    static constexpr double CELL_SIZE = 4.0;
    static constexpr double CHANNEL_WIDTH = 44.0;
    static constexpr double CHANNEL_HEIGHT = 12.0;
    static const int NUM_CELLS_WIDTH = 2 + ceil( CHANNEL_WIDTH / CELL_SIZE);
    static const int NUM_CELLS_HEIGHT = 2 + ceil( CHANNEL_HEIGHT / CELL_SIZE);


    Environment(double dt_, double save_time_)
    {
      dt = dt_;
      save_time = save_time_;
      // Clear the environment.
      for (int x = 0; x < NUM_CELLS_WIDTH; x++)
      {
        for (int y = 0; y < NUM_CELLS_HEIGHT; y++)
        {
          grid_[x][y] = NULL;
        }
      }
    }


    Environment(string filename_paramaters, string filename_agents);
    void add(ABMagent* agent);
    void move(ABMagent* agent, double x_prev, double y_prev);
    void applyForceCell2Cell(ABMagent* agent, ABMagent* other,
                          double point_agent_x, double point_agent_y,
                          double point_other_x, double point_other_y,
                          double distance
                       );
   void applyForceCell2Wall(ABMagent* agent, double point_agent_x,
                          double point_agent_y, double point_wall_y,
                          double separation
                      );
    void handleAgent(ABMagent* agent, ABMagent* other);
    void handleAgentWall(ABMagent* agent, double top);
    void handleCell(int x, int y);
    void handleCellWall(int x, int y, double top);
    void handleCellMoveAndGrow(int x, int y);
    void handleCellSplit(int x, int y);
    void handleInteractions();
    int countNumberAgents();
    void save_data();
    void writeSimulationAgents(string filename);
    void writeSimulationParameters(string filename);
  private:
    ABMagent* grid_[NUM_CELLS_WIDTH][NUM_CELLS_HEIGHT];
};


void Environment::add(ABMagent* agent)
{
  // if cell gone from simulation, delete it
  if ( ( agent->x < 0.0 || agent->x > CHANNEL_WIDTH )
        || ( agent->y < 0.0 || agent->y > CHANNEL_HEIGHT ))
  {
    return;
  }
  // Determine which environment cell it's in.
  int cellX = floor(agent->x / Environment::CELL_SIZE) + 1;
  int cellY = floor(agent->y / Environment::CELL_SIZE) + 1;

  // Add to the front of list for the cell it's in.
  agent->prev_ = NULL;
  agent->next_ = grid_[cellX][cellY];
  grid_[cellX][cellY] = agent;

  if (agent->next_ != NULL)
  {
    agent->next_->prev_ = agent;
  }
}


ABMagent::ABMagent
(
  Environment* environment_,
  double x_,
  double y_,
  double radius_,
  double max_length_,
  double length_,
  double vel_x_,
  double vel_y_,
  double acc_x_,
  double acc_y_,
  double angle_,
  double inertia_,
  double vel_angle_,
  double growth_rate_,
  string label_,
  double force_x_prev_,
  double force_y_prev_,
  double torque_prev_
)
{
  environment = environment_;
  x = x_;
  y = y_;
  radius = radius_;
  max_length = max_length_;
  length = length_;
  vel_x = vel_x_;
  vel_y = vel_y_;
  acc_x = acc_x_;
  acc_y = acc_y_;
  angle = angle_;
  vel_angle = vel_angle_;
  growth_rate = growth_rate_;
  label = label_;
  inertia = inertia_;
  force_x_prev = force_x_prev_;
  force_y_prev = force_y_prev_;
  torque_prev = torque_prev_;

  split_length = max_length * ( 1.0 + uni_length_distribution(len_generator));
  environment->add(this);
}


void ABMagent::print()
{
  cout << "\n ===== Cell info =====";
  cout << "\n label  : " << label;
  cout << "\n length : " << length;
  cout << "\n x      : " << x;
  cout << "\n y      : " << y;
  cout << "\n angle  : " << angle;
  cout << "\n =====================";
}


void Environment::writeSimulationParameters(string filename)
{
  // Saves parameters of the simulation in this order for use later in reloading
  // of simulations.
  ofstream myfile(filename);
  string save_items = "dt, save_time, mass, beta, k_n, gamma_n, gamma_t, "
            "damping_lin, damping_tor, mu_cc, mu_cw, CELL_SIZE, CHANNEL_WIDTH, "
            "CHANNEL_HEIGHT\n";
  myfile << save_items;
  myfile <<  dt << ", " << save_time << ", " << mass << ", " << beta << ", "
          << k_n << ", " << gamma_n << ", " << gamma_t << ", " << damping_lin <<
          ", " << damping_tor << ", " << mu_cc << ", " << mu_cw << ", " <<
          CELL_SIZE << ", " << CHANNEL_WIDTH << ", " << CHANNEL_HEIGHT << "\n";
  myfile.close();
}


void Environment::writeSimulationAgents(string filename)
{
  // Skips line at the end to show new timestep is being saved. Constructor
  // of Environment can reload last saved timestep.
  ofstream fout;
  //ifstream fin;
  //fin.open(filename);
  fout.open(filename, std::ios_base::app);
  ABMagent* agent;
  //if (fin.is_open())
  //{
    for (int x = 1; x < NUM_CELLS_WIDTH; x++)
    {
      for (int y = 1; y < NUM_CELLS_HEIGHT; y++)
      {
        agent = grid_[x][y];
        while (agent != NULL)
        {
          fout << agent->label << ", " << agent->x << ", " << agent->y
                << ", " << agent->angle << ", " << agent->length
                << ", " << agent->radius << ", " << agent->growth_rate
                << ", " << agent->max_length << ", " << agent->split_length
                << ", " << agent->inertia << ", " << agent->vel_x
                << ", " << agent->vel_y << ", " << agent->vel_angle
                << ", " << agent->acc_x << ", " << agent->acc_y << "\n";
          agent = agent->next_;
        }
      }
    }
    fout << "\n";
  //}
  //fin.close();
  fout.close();
}


void Environment::move(ABMagent* agent, double x_prev, double y_prev)
{
  // See which cell it was in.
  int oldCellX = floor(x_prev / Environment::CELL_SIZE) + 1;
  int oldCellY = floor(y_prev / Environment::CELL_SIZE) + 1;

  // If cell leaves simulation, delete it
  if ( ( agent->x < 0.0 || agent->x > CHANNEL_WIDTH )
        || ( agent->y < 0.0 || agent->y > CHANNEL_HEIGHT ))
  {
    if (agent->prev_ != NULL)
    {
      agent->prev_->next_ = agent->next_;
    }

    if (agent->next_ != NULL)
    {
      agent->next_->prev_ = agent->prev_;
    }

    // If it's the head of a list, remove it.
    if (grid_[oldCellX][oldCellY] == agent)
    {
      grid_[oldCellX][oldCellY] = agent->next_;
    }
    return;
  }

  // See which cell it's moving to.
  int cellX = floor(agent->x / Environment::CELL_SIZE) + 1;
  int cellY = floor(agent->y / Environment::CELL_SIZE) + 1;

  // If it didn't change cells, we're done.
  if (oldCellX == cellX && oldCellY == cellY)
  {
    return;
  }

  // Unlink it from the list of its old cell.
  if (agent->prev_ != NULL)
  {
    agent->prev_->next_ = agent->next_;
  }

  if (agent->next_ != NULL)
  {
    agent->next_->prev_ = agent->prev_;
  }

  // If it's the head of a list, remove it.
  if (grid_[oldCellX][oldCellY] == agent)
  {
    grid_[oldCellX][oldCellY] = agent->next_;
  }
  // Add it back to the environment at its new cell.
  add(agent);
}

void ABMagent::addForce(
              double point_x,
              double point_y,
              double ext_force_x,
              double ext_force_y
            )
{
  force_x += ext_force_x;
  force_y += ext_force_y;
  torque += ( point_x - x ) * ext_force_y - ( point_y - y ) * ext_force_x;
}


int Environment::countNumberAgents()
{
  int num_agents = 0;
  ABMagent* agent;
  for (int x = 1; x < NUM_CELLS_WIDTH; x++)
  {
    for (int y = 1; y < NUM_CELLS_HEIGHT; y++)
    {
      // cout << "\n -----> GRID CELL : (i,j) = (" << x << "," << y << ")";
      agent = grid_[x][y];
      while (agent != NULL)
      {
        num_agents++;
        //cout << "\n counting : " << num_agents;
        //cout << "\nlabel : " << agent->label << ", (x,y) : (" << agent->x << "," << agent->y << ")";
        agent = agent->next_;
      }
    }
  }
  cout << "\nNumber of cells: " << num_agents;
  return num_agents;
}


void ABMagent::split()
{
  // Need to give back a bunch of stuff for the daughter cell
  string label_daugh = label;
  label.append("0");
  label_daugh.append("1");

  // location of center of daughter cell and new location of cell
  cout << "\n SPLITTING CELL -----------";
  double x_prev = x;
  double y_prev = y;
  double rand_angle = ang_distribution(ang_generator);
  double angle_daugh = angle - rand_angle;
  double rand_length = len_distribution(len_generator);
  double x_daugh = x - length * cos(angle) * ( 1.0 - 4.0 * rand_length ) / 4.0;
  x += length * cos(angle) * ( 1.0 + 4.0 * rand_length ) / 4.0;
  double y_daugh = y - length * sin(angle) * ( 1.0 - 4.0 * rand_length ) / 4.0;
  y += length * sin(angle) * ( 1.0 + 4.0 * rand_length ) / 4.0;
  // Update length shortened
  angle += rand_angle;
  length /= 2.0;
  double length_daugh = length - 2.0 * rand_length;
  length += 2.0 * rand_length;

  // new CoM velocities

  double vel_x_daugh = vel_x - vel_angle * ( y_daugh - y_prev );
  double vel_y_daugh = vel_y + vel_angle * ( x_daugh - x_prev );
  vel_x += vel_angle * ( y - y_prev );
  vel_y -= vel_angle * ( x - x_prev );
  double vel_angle_daugh = vel_angle;

  force_x_prev = 0.0;
  force_y_prev = 0.0;
  torque_prev = 0.0;

  ABMagent* daughter = new ABMagent(environment, x_daugh, y_daugh, radius, max_length,
                                    length_daugh, vel_x_daugh, vel_y_daugh,
                                    acc_x, acc_y, angle_daugh,
                                    inertia, vel_angle_daugh, growth_rate, label_daugh,
                                    force_x_prev, force_y_prev, torque_prev);

  environment->move(this, x_prev, y_prev);
  environment->split_bool= true;
  print();
  daughter->print();
  //int num_agent = environment->countNumberAgents();
};


void ABMagent::grow(double dt)
{
  length *= exp ( growth_rate * dt );
};


void ABMagent::move_old(double dt, double damping_lin, double damping_tor)
{
  //double reduced_mass = ( 2.0 * (length - 2 * radius ) )
  //                      / ( ( PI * radius ) + 1.0 );
  double reduced_mass = environment->mass * ( ( length - 2 * radius ) * 4
                        / ( 2 * environment->mass * radius) + 1.0);
  double inertia_dt = reduced_mass * pow(length, 2.0) / 12.;
  double vel_x_prev = vel_x;
  double vel_y_prev = vel_y;
  double vel_angle_prev = vel_angle;
  double acc_x_prev = acc_x;
  double acc_y_prev = acc_y;
  double acc_angle_prev = acc_angle;

  // accelerations, implicit linear acceleration method
  double A_1 = reduced_mass * ( 1.0 + dt * damping_lin / 2.0);
  double A_2 = reduced_mass * ( 1.0 - dt * damping_lin / 2.0);
  double A_3 = inertia_dt * ( 1.0 + dt * damping_tor / 2.0);
  double A_4 = inertia_dt * ( 1.0 - dt * damping_tor / 2.0);
  acc_x = ( A_2 * acc_x_prev + force_x - force_x_prev ) / A_1;
  acc_y = ( A_2 * acc_y_prev + force_y - force_y_prev ) / A_1;
  acc_angle = ( A_4 * acc_angle_prev + torque - torque_prev ) / A_3;

  // velocities calculated with damping_lin
  vel_x += dt * ( acc_x + acc_x_prev ) / 2.0;
  vel_y += dt * ( acc_y + acc_y_prev ) / 2.0;
  vel_angle += dt * ( acc_angle + acc_angle_prev ) / 2.0; //
  //vel_angle = dt * ( torque / ( damping_lin * inertia ) );

  // reposition cell and re-orient
  double x_prev = x;
  double y_prev = y;
  x += dt * vel_x_prev + pow(dt, 2.0) * ( acc_x + 2.0 * acc_x_prev ) / 6.0;
  y += dt * vel_y_prev + pow(dt, 2.0) * ( acc_y + 2.0 * acc_y_prev ) / 6.0;
  angle += dt * vel_angle + pow(dt, 2.0) * ( acc_angle + 2.0 * acc_angle_prev ) / 6.0;

  // reset forces to zero for next round
  force_x_prev = force_x;
  force_y_prev = force_y;
  torque_prev = torque;
  force_x = 0.0;
  force_y = 0.0;
  torque = 0.0;
  environment->move(this, x_prev, y_prev);
};


void ABMagent::move(double dt, double damping_lin, double damping_tor)
{
  //double reduced_mass = ( 2.0 * (length - 2 * radius ) )
  //                      / ( ( PI * radius ) + 1.0 );

  double reduced_mass = environment->mass * ( ( length - 2 * radius ) * 4
                        / ( 2 * environment->mass * radius) + 1.0);
  reduced_mass = 1.0;
  double inertia_dt = reduced_mass * pow(length, 2.0) / 12.;

  // reposition cell and re-orient
  double x_prev = x;
  double y_prev = y;
  x += dt * force_x / ( reduced_mass * damping_lin * length );
  y += dt * force_y / ( reduced_mass * damping_lin * length );
  angle += dt * torque / ( inertia_dt * damping_tor * length );

  // reset forces to zero for next round
  force_x = 0.0;
  force_y = 0.0;
  torque = 0.0;
  environment->move(this, x_prev, y_prev);
};


void Environment::applyForceCell2Cell
(
  ABMagent* agent,
  ABMagent* other,
  double point_agent_x,
  double point_agent_y,
  double point_other_x,
  double point_other_y,
  double distance
)
{
  double M_e = mass / 2.0;
  double delta = agent->radius + other->radius - distance;
  double ext_force_x = 0.0;
  double ext_force_y = 0.0;
  double normal_x;
  double normal_y;
  double vel_delta_x;
  double vel_delta_y;
  double vel_norm;
  double tngt_x = 0.0;
  double tngt_y = 0.0;
  double vel_tngt_x;
  double vel_tngt_y;
  double vel_tngt;
  double norm_force = 0.0;
  double tngt_force = 0.0; // tangential force
  double max_friction;
  double friction;

  normal_x = ( point_agent_x - point_other_x ) / distance;
  normal_y = ( point_agent_y - point_other_y ) / distance;

  vel_delta_x = agent->vel_x - other->vel_x +
                  ( agent->vel_angle * ( point_agent_y - agent->y )
                  - other->vel_angle * ( point_other_y - other->y ));
  vel_delta_y = agent->vel_y - other->vel_y +
                  ( - agent->vel_angle * ( point_agent_x - agent->x )
                  + other->vel_angle * ( point_other_x - other->x ));

  // Calculating the normal force
  vel_norm = vel_delta_x * normal_x + vel_delta_y * normal_y;
  norm_force = k_n * pow(delta, 1.5) - gamma_n * M_e * delta * vel_norm;

  // Calculating tangential force
  vel_tngt_x = vel_delta_x - vel_norm * normal_x;
  vel_tngt_y = vel_delta_y - vel_norm * normal_y;
  vel_tngt = sqrt( pow(vel_tngt_x, 2.0) + pow(vel_tngt_y, 2.0) );
  if ( vel_tngt != 0.0)
  {
    tngt_x = vel_tngt_x / vel_tngt;
    tngt_y = vel_tngt_x / vel_tngt;
    friction = gamma_t * M_e * pow(delta, 0.5) * vel_tngt;
    max_friction = mu_cc * norm_force;
    tngt_force = - min(max_friction, friction);
  }

  ext_force_x = norm_force * normal_x + tngt_force * tngt_x;
  ext_force_y = norm_force * normal_y + tngt_force * tngt_y;

  agent->addForce(point_agent_x, point_agent_y, ext_force_x, ext_force_y);
  other->addForce(point_other_x, point_other_y, -ext_force_x, -ext_force_y);
}


void Environment::applyForceCell2Wall
(
  ABMagent* agent,
  double point_agent_x,
  double point_agent_y,
  double point_wall_y,
  double separation
)
{
  double M_e = mass / 2.0;
  double delta = separation;
  double ext_force_x = 0.0;
  double ext_force_y = 0.0;
  double normal_x;
  double normal_y;
  double vel_delta_x;
  double vel_delta_y;
  double vel_norm;
  double tngt_x = 0.0;
  double tngt_y = 0.0;
  double vel_tngt_x;
  double vel_tngt_y;
  double vel_tngt;
  double norm_force = 0.0;
  double tngt_force = 0.0; // tangential force
  double max_friction;
  double friction;

  normal_x = 0.0;
  normal_y = ( point_agent_y - point_wall_y ) / separation;

  vel_delta_x = agent->vel_x + agent->vel_angle * ( point_agent_y - agent->y );
  vel_delta_y = agent->vel_y - agent->vel_angle * ( point_agent_x - agent->x );

  // Calculating the normal force
  vel_norm = vel_delta_x * normal_x + vel_delta_y * normal_y;
  norm_force = k_n * pow(delta, 1.5) - gamma_n * M_e * delta * vel_norm;

  // Calculating tangential force
  vel_tngt_x = vel_delta_x - vel_norm * normal_x;
  vel_tngt_y = vel_delta_y - vel_norm * normal_y;
  vel_tngt = sqrt( pow(vel_tngt_x, 2.0) + pow(vel_tngt_y, 2.0) );
  if ( vel_tngt != 0.0)
  {
    tngt_x = vel_tngt_x / vel_tngt;
    tngt_y = vel_tngt_x / vel_tngt;
    friction = gamma_t * M_e * pow(delta, 0.5) * vel_tngt;
    max_friction = mu_cw * norm_force;
    tngt_force = - min(max_friction, friction);
  }

  ext_force_x = norm_force * normal_x + tngt_force * tngt_x;
  ext_force_y = norm_force * normal_y + tngt_force * tngt_y;

  agent->addForce(point_agent_x, point_agent_y, ext_force_x, ext_force_y);
  // cout << "Walls definitely happen!";
}

void Environment::handleAgent(ABMagent* agent, ABMagent* other)
{
  while (other != NULL)
  {
  double distance_x = other->x - agent->x;
  double distance_y = other->y - agent->y;
  double distance_centers = sqrt( pow(distance_x, 2) + pow(distance_y, 2) );
    // Continue only if centers are a certain distance from eachother
    if ( 2.0 * distance_centers < agent->length + other->length )
    {
      //cout << "\n filter 1";
      // Continue only if there is potential overlap
      double vector_agent_x = ( agent->length - 2.0 * agent->radius ) * cos(agent->angle) / 2.0;
      double vector_agent_y = ( agent->length - 2.0 * agent->radius ) * sin(agent->angle) / 2.0;
      double vector_other_x = ( other->length - 2.0 * other->radius ) * cos(other->angle) / 2.0;
      double vector_other_y = ( other->length - 2.0 * other->radius ) * sin(other->angle) / 2.0;
      double project_agent;
      double project_other;
      project_agent = ( distance_x * vector_agent_x
                        + distance_y * vector_agent_y ) / distance_centers; // I put length distance_length before... why?
      project_other = ( distance_x * vector_other_x
                        + distance_y * vector_other_y ) / distance_centers;
      if (agent->radius + other->radius
                  > distance_centers - abs(project_agent) - abs(project_other) )
      {
        //cout << "\n filter 2";
        double point_agent_x;
        double point_agent_y;
        double point_other_x;
        double point_other_y;
        double line_seg_x1;
        double line_seg_y1;
        double line_seg_x2;
        double line_seg_y2;
        double distance;

        if ( abs(project_agent) >= abs(project_other) )
        {
          point_agent_x = agent->x + sgn(project_agent) * vector_agent_x; // I thought this had to be + ???
          point_agent_y = agent->y + sgn(project_agent) * vector_agent_y;
          line_seg_x1 = other->x - vector_other_x;
          line_seg_x2 = other->x + vector_other_x;
          line_seg_y1 = other->y - vector_other_y;
          line_seg_y2 = other->y + vector_other_y;
          pnt2line(line_seg_x1, line_seg_y1, line_seg_x2, line_seg_y2,
                    point_agent_x, point_agent_y, point_other_x,
                    point_other_y, distance);
        }
        if ( abs(project_agent) < abs(project_other)  )
        {
          point_other_x = other->x - sgn(project_other) * vector_other_x;
          point_other_y = other->y - sgn(project_other) * vector_other_y;
          line_seg_x1 = agent->x - vector_agent_x;
          line_seg_x2 = agent->x + vector_agent_x;
          line_seg_y1 = agent->y - vector_agent_y;
          line_seg_y2 = agent->y + vector_agent_y;
          pnt2line(line_seg_x1, line_seg_y1, line_seg_x2, line_seg_y2,
                    point_other_x, point_other_y, point_agent_x,
                    point_agent_y, distance);
        }
        if ( distance < agent->radius + other->radius)
        {
          //cout << "\n filter 3";
          applyForceCell2Cell(agent, other, point_agent_x, point_agent_y,
                                point_other_x, point_other_y, distance);
        }
      }
    }

    other = other->next_;
  }
}


void Environment::handleAgentWall(ABMagent* agent, double top)
{
  double vector_agent_y = ( agent->length - 2.0 * agent->radius )
                          * sin(agent->angle) / 2.0;
  double vector_agent_x = ( agent->length - 2.0 * agent->radius )
                          * cos(agent->angle) / 2.0;
  double point_agent_x;
  double point_agent_y;
  if ( top == 1.0 )
  {
    double sep_cell_wall = agent->y + abs(vector_agent_y) + agent->radius
                           - CHANNEL_HEIGHT;
    if ( sep_cell_wall > 0.0 )
      {
        point_agent_y = agent->y + abs(vector_agent_y);
        point_agent_x = agent->x + sgn(vector_agent_y)*vector_agent_x;
        applyForceCell2Wall(agent, point_agent_x, point_agent_y, CHANNEL_HEIGHT,
                           sep_cell_wall);
      }
  }
  if ( top == -1.0 )
  {
    double sep_cell_wall =  - agent->y + abs(vector_agent_y) + agent->radius;
    if ( sep_cell_wall > 0.0 )
      {
        point_agent_y = agent->y - abs(vector_agent_y);
        point_agent_x = agent->x - sgn(vector_agent_y)*vector_agent_x;
        applyForceCell2Wall(agent, point_agent_x, point_agent_y, 0.0,
                           sep_cell_wall);
      }
  }
  return;
}

void Environment::handleCell(int x, int y)
{
  ABMagent* agent = grid_[x][y];
  while (agent != NULL)
  {
    // Handle other agents in this cell.
    handleAgent(agent, agent->next_);

    // Hand agents in other cells.
    handleAgent(agent, grid_[x][y+1]);
    handleAgent(agent, grid_[x+1][y-1]);
    handleAgent(agent, grid_[x+1][y]);
    handleAgent(agent, grid_[x+1][y+1]);
    agent = agent->next_;
  }
}

void Environment::handleCellWall(int x, int y, double top)
{
  // top is 1.0 if top, -1.0 if bottom
  ABMagent* agent = grid_[x][y];
  while (agent != NULL)
  {
    // Handle other agents in this cell.
    handleAgentWall(agent, top);
    agent = agent->next_;
  }
}


void Environment::handleCellMoveAndGrow(int x, int y)
{
  ABMagent* agent = grid_[x][y];
  ABMagent* next_agent;
  while (agent != NULL)
  {
    next_agent = agent->next_;
    agent->move(dt, damping_lin, damping_tor);
    agent->grow(dt);
    agent = next_agent;
  }
}


void Environment::handleCellSplit(int x, int y)
{
  ABMagent* agent = grid_[x][y];
  ABMagent* next_agent;
  while (agent != NULL)
  {
    next_agent = agent->next_;
    if (agent->length > agent->split_length) {
      agent->split();
    };
    agent = next_agent;
  }
}


void Environment::handleInteractions()
{
  for (int x = 1; x < NUM_CELLS_WIDTH; x++)
  {
    for (int y = 1; y < NUM_CELLS_HEIGHT; y++)
    {
      handleCell(x, y);
    }
  }
  // Wall interactions to be added shortly.
  for (int x = 1; x < NUM_CELLS_WIDTH; x++)
  {
    handleCellWall(x, 1, -1.0);
    handleCellWall(x, NUM_CELLS_HEIGHT-2, 1.0);
  }

  // Growing and Mouvement
  for (int x = 1; x < NUM_CELLS_WIDTH; x++)
  {
    for (int y = 1; y < NUM_CELLS_HEIGHT; y++)
    {
      handleCellMoveAndGrow(x, y);
    }
  }

  for (int x = 1; x < NUM_CELLS_WIDTH; x++)
  {
    for (int y = 1; y < NUM_CELLS_HEIGHT; y++)
    {
      handleCellSplit(x, y);
    }
  }
}

/*
bool readDataFile (string fname, int layer) {
  int dimension = 1 + 2/environmentRES;
  int row[dimension];
  char delim = ',';
  fstream fin;
  fin.open(fname,ios::in);
  string temp, line, word;
  if(layer > NLAYERS) {
    cout << "Too many input files. Input at most: " << NLAYERS;
    cout << "input files.";
    return false;
  }

  while(fin >> temp) {
    for (int i = 0; i < dimension; i++) {
      getline(fin, line);
      stringstream s(line);
      for (int j = 0; j < dimension; j++) {
        getline(s, word, delim);
  if(word != "") enviro.push(j, i, layer, stoi(word,nullptr,10));
      }
    }
  }
  fin.close();
  return true;
}
*/

int main (int argc, char* argv[]) {
  /*
  if (argc == 1)   // no arguments were passed
  {
    // do things for no arguments

    // usually those variables are set here for a generic flow onwards
    char* default_args[3];
    default_args[0] = argv[0];
    default_args[1] = 0;
    default_args[2] = 0;
    argv = default_args;
    cout << *default_args[2] << "\n";
  }
  cout << "Hi\n";
  cout << "Hi\n";
  cout << argv[1] << " " << argv[2] << "\n";
  cout << argv[1] << " " << argv[2] << "\n";
  */

  // Setup: set the random number generator's seed, initialize our display
  // window.
  double dt = 0.000025; // in minutes
  double save_time = 5.0; // X minutes
  int num_sub_iter = save_time / dt;
  int num_save_iter = 12 * 60 / ( num_sub_iter * dt );
  int num_agents = 0;

  string datafolder = "./data";
  int EXP_NUM = atoi(argv[1]);
  int SIM_NUM = atoi(argv[2]);
  string simulation_name = datafolder + "/" + "c_exp_" + to_string(EXP_NUM) + "/";
  string simulation_param_file = simulation_name + "params.txt";
  string simulation_agent_file = simulation_name + "sim" + to_string(SIM_NUM) + ".txt";
  filesystem::create_directories(simulation_name);

  Environment enviro(dt, save_time);
  enviro.writeSimulationParameters(simulation_param_file);

  // Initialize the cells
  double x = 1.0;
  double y = 3.0;
  double radius = 0.5;
  double max_length = 4.0;
  double length = 2.0;
  //double angle = PI / 4.0;
  double angle = 0;
  double inertia = 5.0;
  double growth_rate = 0.013;

  ABMagent* newBacteriaPntr = new ABMagent(&enviro, 22.0, 6.0, radius,
                    max_length, length, 0.0, 0.0, 0.0, 0.0,
                    0.0, inertia, 0.0, growth_rate, to_string(0),
                    0.0, 0.0, 0.0);
  ABMagent* newBacteriaPntr2 = new ABMagent(&enviro, 23.5, 5.5, radius,
                    max_length, length, 0.0, 0.0, 0.0, 0.0,
                    PI/2.0, inertia, 0.0, growth_rate, to_string(1),
                    0.0, 0.0, 0.0);

  /*
  // reg
  ABMagent* newBacteriaPntr = new ABMagent(&enviro, 23.7014, 4.55602, radius,
                    max_length, length, 0.0, 0.0, 0.0, 0.0,
                    PI/2.0, inertia, 0.0, growth_rate, to_string(0),
                    0.0, 0.0, 0.0);
  ABMagent* newBacteriaPntr2 = new ABMagent(&enviro, 23.8783, 6.455671, radius,
                    max_length, length, 0.0, 0.0, 0.0, 0.0,
                    PI/2.0, inertia, 0.0, growth_rate, to_string(1),
                    0.0, 0.0, 0.0);

  // other
  ABMagent* newBacteriaPntr = new ABMagent(&enviro, 22.0, 6.0, radius,
                    max_length, length, 0.0, 0.0, 0.0, 0.0,
                    0.0, inertia, 0.0, growth_rate, to_string(0),
                    0.0, 0.0, 0.0);
  ABMagent* newBacteriaPntr2 = new ABMagent(&enviro, 23.5, 5.5, radius,
                    max_length, length, 0.0, 0.0, 0.0, 0.0,
                    PI/2.0, inertia, 0.0, growth_rate, to_string(1),
                    0.0, 0.0, 0.0);

  // next
  ABMagent* newBacteriaPntr = new ABMagent(&enviro, 10.0, 4.0, radius,
                    max_length, length, 0.0, 0.0, 0.0, 0.0,
                    PI/3.0, inertia, 0.0, growth_rate, to_string(0),
                    0.0, 0.0, 0.0);
  ABMagent* newBacteriaPntr2 = new ABMagent(&enviro, 31.0, 8.0, radius,
                    max_length, length, 0.0, 0.0, 0.0, 0.0,
                    PI/7.0, inertia, 0.0, growth_rate, to_string(1),
                    0.0, 0.0, 0.0);
  */

  enviro.writeSimulationAgents(simulation_agent_file);
  // Simulation loop
  //num_save_iter = 1;
  //num_sub_iter = 1;
  for (int i = 0 ; i < num_save_iter; i++)
  {
    for (int j = 0 ; j < num_sub_iter; j++)
    {
      enviro.handleInteractions();
    }
    cout << "\n\n-------------\n\n";
    cout << "Number of " << to_string(save_time) << " minutes runs: " << i + 1;
    num_agents = enviro.countNumberAgents();
    enviro.writeSimulationAgents(simulation_agent_file);
  }
  cout << "\n\n-------------\n\n";
  // Anything here won't be run until the window is closed. Want to cout some
  // stats or other information? Write a file? Do something else? put it here.
  return 0;
}
