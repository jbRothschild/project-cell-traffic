#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* exp */
#include <time.h>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <sstream>
#include <list>
#include <vector>
#include <filesystem>
#include <chrono>
#include "intersect.hpp"
#include "bacteria.hpp"
#include "sampled_distribution.hpp"

#define PI 3.14159265359

using namespace std;
using namespace std::chrono;

// desktop : g++ -std=c++17 des.cpp -o des.o
// niagara : module load gcc, g++ -std=c++17 -o des.o des.cpp -lstdc++fs
// laptop : g++-8 -std=c++17 des.cpp -o des.o -lstdc++fs
//
// ./des.o 42 42
random_device rd;  // could set this to a number for a fixed seed
//mt19937 generator(rd() OR 0);
mt19937 generator(rd());
mt19937 len_generator(rd());
mt19937 ang_generator(rd());

float ang = 0.0;
float len = 0.05;
float max_len = 0.1;
uniform_real_distribution<double> len_distribution(-len, len);
uniform_real_distribution<double> ang_distribution(-ang * PI , ang * PI);
uniform_real_distribution<double> uni_length_distribution(-max_len, max_len);

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
                double &dist
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
  dist = sqrt( pow(pnt_x - pnt_on_line_x, 2.0)
                  + pow(pnt_y - pnt_on_line_y, 2.0) );
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
    double mass = 1.0; // 1/5
    double beta = 0.8;
    double k_n = 45000000.0; // 4000000.0
    double gamma_n = 0.0; // 0.00
    double gamma_t = 0.0; // 0.00
    double damping_lin = 200.0 * 60; // 12000
    double damping_tor = 200.0 * 60; // 8000 * PI
    double mu_cc = 0.1;
    double mu_cw = 0.8;
    bool split_bool= false;
    string file_param;
    string file_agents;
    string file_data;
    int nbr_strains;
    static constexpr double CELL_SIZE = 8.0;

    // static constexpr double CHANNEL_WIDTH = 44.0;
    // static constexpr double CHANNEL_HEIGHT = 12.0;
    static constexpr double CHANNEL_WIDTH = 4.56 * 10;
    static constexpr double CHANNEL_HEIGHT = 1.10;

    static const int NUM_CELLS_WIDTH = 2 + ceil( CHANNEL_WIDTH / CELL_SIZE);
    static const int NUM_CELLS_HEIGHT = 2 + ceil( CHANNEL_HEIGHT / CELL_SIZE);
    ABMagent* grid_[NUM_CELLS_WIDTH][NUM_CELLS_HEIGHT];

    Environment(double dt_, double save_time_, string file_param_,
                string file_agents_, string file_data_)
    {
      dt = dt_;
      save_time = save_time_;
      file_param = file_param_;
      file_agents = file_agents_;
      file_data = file_data_;
      // Clear the environment.
      for (int x = 0; x < NUM_CELLS_WIDTH; x++)
      {
        for (int y = 0; y < NUM_CELLS_HEIGHT; y++)
        {
          grid_[x][y] = NULL;
        }
      }
    }

    void add(ABMagent* agent);
    void move(ABMagent* agent, double x_prev, double y_prev);
    void applyForceCell2Cell(ABMagent* agent, ABMagent* other,
                          double point_agent_x, double point_agent_y,
                          double point_other_x, double point_other_y,
                          double dist
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
    void writeSimulationAgents();
    bool writeSimulationData();
    void writeSimulationParameters();
  private:
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


void Environment::writeSimulationParameters()
{
  // Saves parameters of the simulation in this order for use later in reloading
  // of simulations.
  ofstream myfile(file_param);
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


void Environment::writeSimulationAgents()
{
  // Skips line at the end to show new timestep is being saved. Constructor
  // of Environment can reload last saved timestep.
  ofstream fout;
  //ifstream fin;
  //fin.open(filename);
  fout.open(file_agents, std::ios_base::app);
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


int countDigit(int n)
{
    if (n/10 == 0)
        return 1;
    return 1 + countDigit(n / 10);
}


bool Environment::writeSimulationData()
{
  // Skips line at the end to show new timestep is being saved. Constructor
  // of Environment can reload last saved timestep.
  ofstream fout;
  int digits = countDigit(nbr_strains);
  fout.open(file_data, std::ios_base::app);
  ABMagent* agent;
  double tot_length [nbr_strains] = {0}; // TODO :replace 2 with number starting strains
  int tot_count [nbr_strains] = {0};
  bool extinct = false;
  // Data acquisition
  for (int x = 1; x < NUM_CELLS_WIDTH; x++)
  {
    for (int y = 1; y < NUM_CELLS_HEIGHT; y++)
    {
      agent = grid_[x][y];
      while (agent != NULL)
      {
        tot_length[stoi( (agent->label).substr(0, digits) )] += agent->length;
        tot_count[stoi( (agent->label).substr(0, digits) )] += 1;
        agent = agent->next_;
      }
    }
  }
  string str_len = "[";
  for (int i = 0; i < nbr_strains; i++)
  {
    str_len += to_string(tot_length[i]);
    str_len += ",";
    if (tot_count[i] == 0)
    {
      extinct = true;
    }
  }
  str_len.pop_back();
  str_len += "]";
  fout << str_len << "\n";
  fout.close();

  return extinct;
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
  // cout << "\nNumber of cells: " << num_agents;  // comment on niagara
  return num_agents;
}


void ABMagent::split()
{
  // Need to give back a bunch of stuff for the daughter cell
  string label_daugh = label;
  label.append("0");
  label_daugh.append("1");

  // location of center of daughter cell and new location of cell
  // cout << "\n SPLITTING CELL -----------";
  double x_prev = x;
  double y_prev = y;
  double rand_angle = ang_distribution(ang_generator);
  double angle_daugh = angle - rand_angle;
  double rand_length = len_distribution(len_generator);
  double x_daugh = x - length * cos(angle) * ( 1.0 + 4.0 * rand_length ) / 4.0 ;
  x += length * cos(angle) * ( 1.0 - 4.0 * rand_length ) / 4.0;
  double y_daugh = y - length * sin(angle) * ( 1.0 + 4.0 * rand_length ) / 4.0 ;
  y += length * sin(angle) * ( 1.0 - 4.0 * rand_length ) / 4.0;
  // Update length shortened
  angle += rand_angle;
  length /= 2.0;
  double length_daugh = length * ( 1.0 - 4.0 * rand_length );
  length += 4.0 * rand_length * length;

  // new CoM velocities

  double vel_x_daugh = vel_x - vel_angle * ( y_daugh - y_prev );
  double vel_y_daugh = vel_y + vel_angle * ( x_daugh - x_prev );
  vel_x += vel_angle * ( y - y_prev );
  vel_y -= vel_angle * ( x - x_prev );
  double vel_angle_daugh = vel_angle;

  force_x_prev = 0.0;
  force_y_prev = 0.0;
  torque_prev = 0.0;

  environment->move(this, x_prev, y_prev);

  ABMagent* daughter = new ABMagent(environment, x_daugh, y_daugh, radius, max_length,
                                    length_daugh, vel_x_daugh, vel_y_daugh,
                                    acc_x, acc_y, angle_daugh,
                                    inertia, vel_angle_daugh, growth_rate, label_daugh,
                                    force_x_prev, force_y_prev, torque_prev);

  environment->split_bool= true;
  //print(); // Information about the cell when split
  //daughter->print();
  //int num_agent = environment->countNumberAgents(); // number of cells after split
};


void ABMagent::grow(double dt)
{
  length *= exp ( growth_rate * dt );
};


void ABMagent::move(double dt, double damping_lin, double damping_tor)
{
  double reduced_mass = environment->mass * radius * ( 2.0 * ( length - 2 * radius )
                        + ( PI * radius ) );
  reduced_mass = 1.0;

  // reposition cell and re-orient
  double x_prev = x;
  double y_prev = y;
  x += dt * force_x / ( reduced_mass * damping_lin * length );
  y += dt * force_y / ( reduced_mass * damping_lin * length );
  double inertia_dt = reduced_mass * pow(length, 2.0) / 12.;
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
  double dist
)
{
  double M_e = mass / 2.0;
  double delta = agent->radius + other->radius - dist;
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

  normal_x = ( point_agent_x - point_other_x ) / dist;
  normal_y = ( point_agent_y - point_other_y ) / dist;

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
  /*
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
  */
  ext_force_x = norm_force * normal_x; //+ tngt_force * tngt_x;
  ext_force_y = norm_force * normal_y; //+ tngt_force * tngt_y;

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
  /*
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
  */
  ext_force_x = norm_force * normal_x + tngt_force * tngt_x;
  ext_force_y = norm_force * normal_y + tngt_force * tngt_y;

  agent->addForce(point_agent_x, point_agent_y, ext_force_x, ext_force_y);
  // cout << "Walls definitely happen!";
}

void Environment::handleAgent(ABMagent* agent, ABMagent* other)
{
  while (other != NULL)
  {
  double dist_x = other->x - agent->x;
  double dist_y = other->y - agent->y;
  double dist_centers = sqrt( pow(dist_x, 2) + pow(dist_y, 2) );
    // Continue only if centers are a certain dist from eachother
    // cout << "\n filter 1";
    if ( 2.0 * dist_centers < agent->length + other->length )
    {
      // Continue only if there is potential overlap, that is to say if the
      // projection onto eachother
      // cout << "\n filter 2";
      double vector_agent_x = ( agent->length - 2.0 * agent->radius ) * cos(agent->angle) / 2.0;
      double vector_agent_y = ( agent->length - 2.0 * agent->radius ) * sin(agent->angle) / 2.0;
      double vector_other_x = ( other->length - 2.0 * other->radius ) * cos(other->angle) / 2.0;
      double vector_other_y = ( other->length - 2.0 * other->radius ) * sin(other->angle) / 2.0;
      double project_agent;
      double project_other;
      project_agent = ( dist_x * vector_agent_x
                        + dist_y * vector_agent_y ) / dist_centers; // I put length dist_length before... why?
      project_other = ( dist_x * vector_other_x
                        + dist_y * vector_other_y ) / dist_centers;
      if (agent->radius + other->radius
                  > dist_centers - abs(project_agent) - abs(project_other) )
      {
        double point_agent_x[4];
        double point_agent_y[4];
        double point_other_x[4];
        double point_other_y[4];
        double X [4] = {other->x + vector_other_x, other->x - vector_other_x,
                         agent->x + vector_agent_x, agent->x - vector_agent_x};
        double Y [4] = {other->y + vector_other_y, other->y - vector_other_y,
                        agent->y + vector_agent_y, agent->y - vector_agent_y};
        double dist[4];
        int idx_min;
        for (int i = 0; i<4; i++)
        {
          if ( i > 1 )
          {
            point_agent_x[i] = X[i];
            point_agent_y[i] = Y[i];
            pnt2line(X[0], Y[0], X[1], Y[1], point_agent_x[i], point_agent_y[i],
                     point_other_x[i], point_other_y[i], dist[i]);
          }
          else
          {
            point_other_x[i] = X[i];
            point_other_y[i] = Y[i];
            pnt2line(X[2], Y[2], X[3], Y[3], point_other_x[i], point_other_y[i],
                     point_agent_x[i], point_agent_y[i], dist[i]);
          }
        }
        idx_min = distance(begin(dist), min_element(begin(dist), end(dist)));

        if ( dist[idx_min] < agent->radius + other->radius)
        {

          //cout << "\n filter 3";
          applyForceCell2Cell(agent, other, point_agent_x[idx_min],
                              point_agent_y[idx_min], point_other_x[idx_min],
                              point_other_y[idx_min], dist[idx_min]);
        }
        /*
        if ( (agent->label == problem && other->label == problem2) ||
              (other->label == problem && agent->label == problem2) )
        {
          cout << "other.\n";
          cout << "Cell " << agent->label << ", p1 : ( " << X[2]<< " , "
              << Y[2] << " ) and p2 : ( " << X[3] << " , " << Y[3] << " ).\n";
          cout << "Cell " << other->label << ", p1 : ( " << X[0]<< " , "
              << Y[0] << " ) and p2 : ( " << X[1] << " , " << Y[1] << " ).\n";
          cout << "Contact " << agent->label <<  " : ( " << point_agent_x[idx_min]
              << " , " << point_agent_y[idx_min] << " ) and " << other->label
              <<  " : ( " << point_other_x[idx_min] << " , " << point_other_y[idx_min] << " ).\n";
          cout << dist[idx_min] << "\n";
        }
        */
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


int initialize_cells_load_relabel(Environment &enviro, string filename,
                                    int timepoint, int nbr_zeros)
{
  int n = 0;
  int i = 0;
  string new_label;
  string label_i;
  ifstream file(filename);
  string str;

  while ( getline(file, str) )
  {
    if (str.empty())
    {
      n+=1;
    }
    else if (timepoint==n)
    {
      // vector of the cell information, parse comma separated string
      vector<string> vect;
      stringstream ss(str);
      string word;
      while (getline(ss, word, ',')) {
        vect.push_back(word);
      }
      label_i = to_string(i);
      int len_label = label_i.length();
      new_label = string(nbr_zeros - std::min(nbr_zeros, len_label), '0')
                          + label_i;
      ABMagent* newBacteriaPntr = new ABMagent(&enviro,
                        stod(vect[1]),
                        stod(vect[2]),
                        stod(vect[5]),
                        stod(vect[7]),
                        stod(vect[4]),
                        0.0, 0.0, 0.0, 0.0,
                        stod(vect[3]),
                        stod(vect[8]),
                        0.0,
                        stod(vect[6]),
                        new_label,
                        0.0, 0.0, 0.0);
      i++;
    }
    if ( timepoint < n ){
      break;
    }
  }
  enviro.writeSimulationAgents();
  enviro.nbr_strains = i;
  enviro.writeSimulationData();
  return i;
}


int initialize_cells_load(Environment &enviro, string filename, int timepoint)
{
  int n = 0;
  ifstream file(filename);
  string str;
  int nbr_strains = 2;

  while ( getline(file, str) )
  {
    if (str.empty())
    {
      n+=1;
    }
    else if (timepoint==n)
    {
      // vector of the cell information, parse comma separated string
      vector<string> vect;
      stringstream ss(str);
      string word;
      while (getline(ss, word, ',')) {
        vect.push_back(word);
      }
      ABMagent* newBacteriaPntr = new ABMagent(&enviro,
                        stod(vect[1]),
                        stod(vect[2]),
                        stod(vect[5]),
                        stod(vect[7]),
                        stod(vect[4]),
                        0.0, 0.0, 0.0, 0.0,
                        stod(vect[3]),
                        stod(vect[8]),
                        0.0,
                        stod(vect[6]),
                        vect[0],
                        0.0, 0.0, 0.0);
    }
    if ( timepoint < n ){
      break;
    }
  }
  enviro.writeSimulationAgents();
  enviro.nbr_strains = nbr_strains;
  enviro.writeSimulationData();
  return nbr_strains;
}

int initialize_cells2(Environment &enviro, int SIM_NUM) {
  // initialize 2 different cells from different types of cells
  BSubt bacteria1;
  BSubt bacteria2;
  double length1 = bacteria1.max_length/2.0;
  double length2 = bacteria2.max_length/2.0;
  double x1, x2, y1, y2, angle1, angle2;
  double inertia = 5.0;
  mt19937 cell_placement(SIM_NUM);

  uniform_real_distribution<double> x_dist(0 , enviro.CHANNEL_WIDTH);
  uniform_real_distribution<double> y1_dist(length1/2.0,
                         enviro.CHANNEL_HEIGHT - length1/2.0);
  uniform_real_distribution<double> y2_dist(length2/2,
                         enviro.CHANNEL_HEIGHT - length2/2);
  uniform_real_distribution<double> init_angle(0, 2*PI);

  bool intersect = true;
  angle1 = init_angle(cell_placement);
  angle2 = init_angle(cell_placement);

  while (intersect)
  {
     x1 = x_dist(cell_placement);
     x2 = x_dist(cell_placement);

     y1 = y1_dist(cell_placement);
     y2 = y2_dist(cell_placement);
     intersect = pointsIntersect(x1 + length1 / 2 * cos(angle1),
                                 x1 - length1 / 2 * cos(angle1),
                                 x2 + length2 / 2 * cos(angle2),
                                 x2 - length2 / 2 * cos(angle2),
                                 y1 + length1 / 2 * sin(angle1),
                                 y1 - length1 / 2 * sin(angle1),
                                 y2 + length2 / 2 * sin(angle2),
                                 y2 - length2 / 2 * sin(angle2)
                               );
  }

  ABMagent* newBacteriaPntr = new ABMagent(&enviro, x1, y1,
                    bacteria1.radius,
                    bacteria1.max_length,
                    bacteria1.length,
                    0.0, 0.0, 0.0, 0.0,
                    angle1,
                    bacteria1.inertia,
                    0.0,
                    bacteria1.growth_rate,
                    to_string(0),
                    0.0, 0.0, 0.0);
  ABMagent* newBacteriaPntr2 = new ABMagent(&enviro, x2, y2,
                    bacteria2.radius,
                    bacteria2.max_length,
                    bacteria2.length,
                    0.0, 0.0, 0.0, 0.0,
                    angle2,
                    bacteria2.inertia,
                    0.0,
                    bacteria2.growth_rate,
                    to_string(1),
                    0.0, 0.0, 0.0);
  enviro.writeSimulationAgents();
  enviro.nbr_strains = 2;
  enviro.writeSimulationData();

  return 2;
}


int initialize_N_strains(Environment &enviro, int SIM_NUM, int nbr_strains_WT,
                         int nbr_strains_A22, int nbr_strains_bsub) {

  // Initialize the cells
  EColi ecoliWT;
  EColiA22 ecoliA22;
  BSubt bsubt;
  ABMagent* agent;
  double x, y, angle, length;
  bool intersect = true;
  bool check_intersect;
  mt19937 cell_placement(SIM_NUM);

  auto sinFunc = [](double x) {
      const double x0 = 0.0;
      const double L = 44.0;
      return x + 0.75 * L * std::sin( 2. * PI * ( x - x0) / L  - PI / 4.0 ) / ( 2.0 * PI );
    }; // PDF(x) = 1 + 0.75cos(2 pi (x-x0) / K - pi / 4 )

  //uniform_real_distribution<double> x_dist(0 , enviro.CHANNEL_WIDTH);
  Sampled_distribution<> x_dist(sinFunc, 0.0, enviro.CHANNEL_WIDTH);
  uniform_real_distribution<double> y_WT_dist(ecoliWT.max_length/2.0,
                         enviro.CHANNEL_HEIGHT - ecoliWT.max_length/2.0);
  uniform_real_distribution<double> y_A22_dist(ecoliA22.max_length/2,
                         enviro.CHANNEL_HEIGHT - ecoliA22.max_length/2);
  uniform_real_distribution<double> init_angle(0, 2*PI);


  // placing wild type ecoli
  for (int i = 0 ; i < nbr_strains_WT; i++){
    length = ecoliWT.max_length / 2.0;
    intersect = true;
    while (intersect)
    {
       intersect = false;
       angle = init_angle(cell_placement);
       x = x_dist(cell_placement);
       y = y_WT_dist(cell_placement);

       // loop to place cells on grid
       for (int x = 1; x < enviro.NUM_CELLS_WIDTH; x++)
       {
         for (int y = 1; y < enviro.NUM_CELLS_HEIGHT; y++)
         {
           agent = enviro.grid_[x][y];

           // check that cells don't intersect
           while (agent != NULL)
           {
             check_intersect= pointsIntersect(x + length / 2 * cos(angle),
                                              x - length / 2 * cos(angle),
                                              agent->x + agent->length / 2 * cos(agent->angle),
                                              agent->x - agent->length / 2 * cos(agent->angle),
                                              y + length / 2 * sin(angle),
                                              y - length / 2 * sin(angle),
                                              agent->y + agent->length / 2 * sin(agent->angle),
                                              agent->y - agent->length / 2 * sin(agent->angle)
                                             );
             intersect = (intersect || check_intersect);
             agent = agent->next_;
           }
         }
       }
    }

    // after checks are passed, actually initialize new cell
    ABMagent* newBacteriaPntr = new ABMagent(&enviro, x, y, ecoliWT.radius,
                                             ecoliWT.max_length,
                                             ecoliWT.length,
                                             0.0, 0.0, 0.0, 0.0,
                                             angle,
                                             ecoliWT.inertia,
                                             0.0,
                                             ecoliWT.growth_rate,
                                             to_string(i),
                                             0.0, 0.0, 0.0
                                            );
  }

  // placing A22 mutant of ecoli
  for (int i = 0 ; i < nbr_strains_A22; i++){
    length = 2.0; // ecoliA22.max_length/2.0; // have same init. vol. all cells
    intersect = true;
    while (intersect)
    {
       intersect = false;
       angle = init_angle(cell_placement);
       x = x_dist(cell_placement);
       y = y_A22_dist(cell_placement);
       for (int x = 1; x < enviro.NUM_CELLS_WIDTH; x++)
       {
         for (int y = 1; y < enviro.NUM_CELLS_HEIGHT; y++)
         {
           agent = enviro.grid_[x][y];
           while (agent != NULL)
           {
             check_intersect= pointsIntersect(x + length / 2 * cos(angle),
                                              x - length / 2 * cos(angle),
                                              agent->x + agent->length / 2 * cos(agent->angle),
                                              agent->x - agent->length / 2 * cos(agent->angle),
                                              y + length / 2 * sin(angle),
                                              y - length / 2 * sin(angle),
                                              agent->y + agent->length / 2 * sin(agent->angle),
                                              agent->y - agent->length / 2 * sin(agent->angle)
                                             );
             intersect = (intersect || check_intersect);
             agent = agent->next_;
           }
         }
       }
    }
    ABMagent* newBacteriaPntr = new ABMagent(&enviro, x, y, ecoliA22.radius,
                                             ecoliA22.max_length,
                                             ecoliA22.length,
                                             0.0, 0.0, 0.0, 0.0,
                                             angle,
                                             ecoliA22.inertia,
                                             0.0,
                                             ecoliA22.growth_rate,
                                             to_string(nbr_strains_WT + i),
                                             0.0, 0.0, 0.0
                                            );
  }

  enviro.writeSimulationAgents();
  enviro.nbr_strains = nbr_strains_WT + nbr_strains_A22 + nbr_strains_bsub;
  enviro.writeSimulationData();

  return enviro.nbr_strains;
}


int initialize_1boundary(Environment &enviro, int EXP_NUM, int bndry_nbr) {
  // EXP_NUM has to be nbr_cells + f * 10; f is the fraction cells on the left

  // Initialize the cells
  EColi bacteria;
  double length = bacteria.max_length;
  double x, y;
  double angle = 0.0;
  double inertia = 5.0;

  // count number
  int nbr_cells = enviro.CHANNEL_WIDTH / length;
  int boundary_pos = EXP_NUM - nbr_cells;
  int cell_boundary = nbr_cells * boundary_pos / 10;
  y = 0.55;

  // add cells
  for (int i = 1; i <=cell_boundary; i++)
  {
    x = length * (i - 0.5);
    ABMagent* newBacteriaPntr = new ABMagent(&enviro, x, y,
                      bacteria.radius,
                      bacteria.max_length,
                      length,
                      0.0, 0.0, 0.0, 0.0,
                      angle,
                      bacteria.inertia,
                      0.0,
                      bacteria.growth_rate,
                      to_string(0),
                      0.0, 0.0, 0.0);
  }

  for (int i = cell_boundary + 1; i <=nbr_cells; i++)
  {
    x = length * (i - 0.5);
    ABMagent* newBacteriaPntr = new ABMagent(&enviro, x, y,
                      bacteria.radius,
                      bacteria.max_length,
                      length,
                      0.0, 0.0, 0.0, 0.0,
                      angle,
                      bacteria.inertia,
                      0.0,
                      bacteria.growth_rate,
                      to_string(1),
                      0.0, 0.0, 0.0);
  }
  enviro.writeSimulationAgents();
  enviro.nbr_strains = 2;
  enviro.writeSimulationData();

  return 2;
}


int main (int argc, char* argv[]) {

  // setup simulation parameters
  double dt = 0.00025; // in minutes 0.000025
  double save_time = 5.0; // X minutes
  int nbr_hours = 12;

  // metadata of simulations
  string datafolder = "./data";
  int EXP_NUM = atoi(argv[1]);
  int SIM_NUM = atoi(argv[2]);
  string sim_name = datafolder + "/" + "c_exp_" + to_string(EXP_NUM) + "/";
  string sim_param_file = sim_name + "params.txt";
  string sim_agent_file = sim_name + "sim" + to_string(SIM_NUM) + ".txt";
  string sim_data_file = sim_name + "sim" + to_string(SIM_NUM) + "_data.txt";
  filesystem::create_directories(sim_name);

  // initialize environment
  Environment enviro(dt, save_time, sim_param_file, sim_agent_file, sim_data_file);
  enviro.writeSimulationParameters();

  // Josh stuff
  // enviro.nbr_strains = initialize_N_strains(enviro, SIM_NUM, 2, 0, 0); // WT, A22, Bsub

  // Boundary stuff
  enviro.nbr_strains = initialize_1boundary(enviro, EXP_NUM, SIM_NUM);

  // Load previous simulation time
  //enviro.nbr_strains = initialize_cells_load(enviro, datafolder + "/c_exp_0/sim1.txt", 11);
  //enviro.nbr_strains = initialize_cells_load_relabel(enviro, datafolder + "/c_exp_11/sim" + argv[2] + ".txt", 193, 3);

  // timing simulation run
  auto start = high_resolution_clock::now();

  // fixation check
  bool fixate = false;

  // simulation loop times
  int num_sub_iter = save_time / dt;
  int num_save_iter = nbr_hours * 60 / ( num_sub_iter * dt );

  // quick check of simulation
  //num_save_iter = 1;
  //num_sub_iter = 1;

  // simulation loop
  for (int i = 0 ; i < num_save_iter; i++)
  {
    for (int j = 0 ; j < num_sub_iter; j++)
    {
      enviro.handleInteractions();
    }
    // saving data
    enviro.writeSimulationAgents();
    fixate = enviro.writeSimulationData();

    // ends simulation when one species goes extinct
    if (fixate) {break;} // comment out for multispecies stuff

    /* //uncomment when not on niagara to track progress
    cout << "\n\n-------------\n\n";
    cout << "Number of " << to_string(save_time) << " minutes runs: " << i + 1;
    // num_agents = enviro.countNumberAgents();
    */
  }
  auto stop = high_resolution_clock::now();

  auto duration = duration_cast<seconds>(stop - start);

  cout << duration.count() << endl;

  return 0;
}
