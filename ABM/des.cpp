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

#define PI 3.14159265359
#define SPEED 0.005
#define NAGENTS 20000
#define XCHANNEL 45.0
#define YCHANNEL 12.0
#define MAXLENCELL 4.0

using namespace std;
random_device rd;
mt19937 generator(rd());
mt19937 len_generator(rd());
mt19937 ang_generator(rd());
uniform_real_distribution<double> len_distribution(-0.025, 0.025);
uniform_real_distribution<double> ang_distribution(-0.02 * PI , 0.02 * PI);
uniform_real_distribution<double> uni_length_distribution(-0.025, 0.025);

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
  double AB_x = B_x - A_x;
  double AB_y = B_y - A_y;
  double line_length = sqrt( pow(AB_x, 2.0) + pow(AB_y, 2.0) );
  double Apnt_x = pnt_x - A_x;
  double Apnt_y = pnt_y - A_y;
  double t = ( AB_x * Apnt_x + AB_y * Apnt_x ) / pow(line_length, 2.0);
  if(t < 0)
  {
    t = 0;
  }
  if(t > 1)
  {
    t = 1;
  }
  pnt_on_line_x = A_x + t * AB_x;
  pnt_on_line_y = A_y + t * AB_y;
  distance = sqrt( pow(pnt_x - pnt_on_line_x, 2.0)
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
    double vel_angle;
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
      double angle_,
      double inertia_,
      double vel_angle_,
      double growth_rate_,
      string label_
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
    void move(double dt, double damping);

};


class Environment
{
  public:
    double dt;
    double mass = 1.0;
    double beta = 0.8;
    double k_n = 333.0;
    double gamma_n = 0.0022;
    double gamma_t = 0.0022;
    double damping = 0.8;
    double mu_cc = 0.1;
    bool split_bool= false;
    static constexpr double CELL_SIZE = 4.0;
    static constexpr double CHANNEL_WIDTH = 44.0;
    static constexpr double CHANNEL_HEIGHT = 12.0;
    static const int NUM_CELLS_WIDTH = 2 + ceil( CHANNEL_WIDTH / CELL_SIZE);
    static const int NUM_CELLS_HEIGHT = 2 + ceil( CHANNEL_HEIGHT / CELL_SIZE);
    Environment(double dt_)
    { dt = dt_;
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
    void moveTest(ABMagent* agent, double x_prev, double y_prev);
    void applyForceCell2Cell(ABMagent* agent, ABMagent* other, double point_agent_x,
                          double point_agent_y, double point_other_x,
                          double point_other_y, double distance
                       );
    void handleAgent(ABMagent* agent, ABMagent* other);
    void handleCell(int x, int y);
    void handleCellMoveAndGrow(int x, int y);
    void handleCellSplit(int x, int y);
    void handleInteractions();
    int countNumberAgents();
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
  double angle_,
  double inertia_,
  double vel_angle_,
  double growth_rate_,
  string label_
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
  angle = angle_;
  vel_angle = vel_angle_;
  growth_rate = growth_rate_;
  label = label_;
  inertia = inertia_;

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


void Environment::moveTest(ABMagent* agent, double x_prev, double y_prev)
{
  // See which cell it was in.
  int oldCellX = floor(x_prev / Environment::CELL_SIZE) + 1;
  int oldCellY = floor(y_prev / Environment::CELL_SIZE) + 1;
  //agent->print();

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
  torque += ( point_x - x ) * ext_force_y - ( point_y - y ) * ext_force_y;
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
  double angle_daugh = angle + ang_distribution(ang_generator);
  double rand_length = len_distribution(len_generator);
  double x_daugh = x - length * cos(angle) * ( 1.0 - 4.0 * rand_length ) / 4.0;
  x += length * cos(angle) * ( 1.0 + 4.0 * rand_length ) / 4.0;
  double y_daugh = y - length * sin(angle) * ( 1.0 - 4.0 * rand_length ) / 4.0;
  y += length * sin(angle) * ( 1.0 + 4.0 * rand_length ) / 4.0;
  // Update length shortened
  length /= 2.0;
  double length_daugh = length - 2.0 * rand_length;
  length += 2.0 * rand_length;

  // new CoM velocities
  double vel_x_daugh = vel_x + vel_angle * ( y_daugh - y_prev );
  double vel_y_daugh = vel_y - vel_angle * ( x_daugh - x_prev) ;
  vel_x += vel_angle * ( y - y_prev );
  vel_y -= vel_angle * ( x - x_prev );

  ABMagent* daughter = new ABMagent(environment, x_daugh, y_daugh, radius, max_length,
                                    length_daugh, vel_x_daugh, vel_y_daugh, angle_daugh,
                                    inertia, vel_angle, growth_rate, label_daugh);

  //daughter->print();
  //print();
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


void ABMagent::move(double dt, double damping)
{
  double reduced_mass = ( 2.0 * (length - 2 * radius )) / ( ( PI * radius ) + 1.0 );
  // velocities calculated with damping
  vel_x += dt * ( force_x / reduced_mass - damping * vel_x );
  vel_y += dt * ( force_y / reduced_mass - damping * vel_y );
  vel_angle += dt * ( torque / inertia - damping * vel_angle );

  // reposition cell and re-orient
  double x_prev = x;
  double y_prev = y;
  x += dt * vel_x;
  y += dt * vel_y;
  angle += dt * vel_angle;

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


void Environment::handleAgent(ABMagent* agent, ABMagent* other)
{
  while (other != NULL)
  {
  double distance_x = other->x - agent->x;
  double distance_y = other->y - agent->y;
  double distance_centers = sqrt( pow(distance_x, 2) + pow(distance_y, 2) );
    // Continue only if centers are a certain distance from eachother
    if ( 2.0 * distance_centers <
        agent->length + other->length
        + 2.0 * agent->radius + 2.0 * other->radius )
    {
      // Continue only if there is potential overlap
      double vector_agent_x = agent->length * cos(agent->angle) / 2.0;
      double vector_agent_y = agent->length * sin(agent->angle) / 2.0;
      double vector_other_x = other->length * cos(other->angle) / 2.0;
      double vector_other_y = other->length * sin(other->angle) / 2.0;
      double project_agent;
      double project_other;
      project_agent = ( distance_x * vector_agent_x
                        + distance_y * vector_agent_y ) / distance_centers; // I put length distance_length before... why?
      project_other = ( distance_x * vector_other_x
                        + distance_y * vector_other_y ) / distance_centers;
      if (agent->radius + other->radius
                  < distance_centers - abs(project_agent) - abs(project_other) )
      {
        double point_agent_x;
        double point_agent_y;
        double point_other_x;
        double point_other_y;
        double line_seg_x1;
        double line_seg_y1;
        double line_seg_x2;
        double line_seg_y2;
        double distance;

        if (abs(project_agent) > abs(project_other))
        {
          point_agent_x = agent->x + signbit(project_agent) * vector_agent_x;
          point_agent_y = agent->y + signbit(project_agent) * vector_agent_y;
          line_seg_x1 = other->x - vector_other_x;
          line_seg_x2 = other->x + vector_other_x;
          line_seg_y1 = other->y - vector_other_y;
          line_seg_y2 = other->y + vector_other_y;
          pnt2line(line_seg_x1, line_seg_y1, line_seg_x2, line_seg_y2,
                    point_other_x, point_other_y, point_agent_x,
                    point_agent_y, distance);
        }
        else
        {
          point_other_x = other->x - signbit(project_other) * vector_other_x;
          point_other_y = other->y - signbit(project_other) * vector_other_y;
          line_seg_x1 = agent->x - vector_agent_x;
          line_seg_x2 = agent->x + vector_agent_x;
          line_seg_y1 = agent->y - vector_agent_y;
          line_seg_y2 = agent->y + vector_agent_y;
          pnt2line(line_seg_x1, line_seg_y1, line_seg_x2, line_seg_y2,
                    point_agent_x, point_agent_y, point_other_x,
                    point_other_y, distance);
        }
        if ( distance < agent->radius + other->radius)
        {
          applyForceCell2Cell(agent, other, point_agent_x, point_agent_y,
                                point_other_x, point_other_y, distance);
        }
      }
    }

    other = other->next_;
  }
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
    /*
    if (split_bool)
    {
      agent->print();
    }
    */
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
    agent->move(dt, damping);
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
  /* Wall interactions to be added shortly.
  for (int y = 1; y < NUM_CELLS_WIDTH; y++)
  {
    handleCellWall(grid_[1][y]);
    handleCellWall(grid_[NUM_CELLS_HEIGHT-1][y]);
    handleCellWall(grid_[NUM_CELLS_HEIGHT-2][y]);
  }
  */
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
int main (int argc, char **argv) {

  // Setup: set the random number generator's seed, initialize our display
  // window.
  double dt = 0.0001; // in minutes
  Environment enviro(dt);

  double save_time = 10.0;
  int num_sub_iter = save_time / dt;
  int num_save_iter = 24 * 60 / ( num_sub_iter * dt );

  double x = 1.0;
  double y = 3.0;
  double radius = 0.5;
  double max_length = 4.0;
  double length = 2.0;
  //double angle = PI / 4.0;
  double angle = 0;
  double inertia = 5.0;
  double growth_rate = 0.013;
  /*
  for (int i = 0; i < 10 ; i++)
  {
    for (int j = 0; j < 1 ; j++)
    {
      ABMagent* newBacteriaPntr = new ABMagent(&enviro, x + 1.0*i, y + 4.0*j, radius, max_length, length, 0.0, 0.0,
                        angle, inertia, 0.0, growth_rate, to_string(i + i*j));
    }
  }
  */
  BMagent* newBacteriaPntr = new ABMagent(&enviro, 8, 6, radius, max_length, length, 0.0, 0.0,
                    angle, inertia, 0.0, growth_rate, to_string(i + i*j));
  int num_agent = enviro.countNumberAgents();

  for (int i = 0 ; i < num_save_iter; i++)
  {
    for (int j = 0 ; j < num_sub_iter; j++)
    {
      enviro.handleInteractions();
      /*
      if (enviro.split_bool)
      {
        cout << "\n"<< j;
        cout << "\nnext";
      }
      */
    }

    cout << "\n\n-------------\n\n";
    cout << "Number of 10 minutes run: " << i;
    num_agent = enviro.countNumberAgents();
  }

  // Anything here won't be run until the window is closed. Want to cout some
  // stats or other information? Write a file? Do something else? put it here.
  return 0;
}
