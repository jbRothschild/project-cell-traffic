#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* exp */
#include <time.h>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <list>
#include <vector>

#define GRIDRES 0.03
#define PI 3.14159265359
#define SPEED 0.005
#define NAGENTS 20000
#define XCHANNEL 45.0
#define YCHANNEL 12.0
#define MAXLENCELL 4.0

using namespace std;


class ABMagent
{
  // Our main actor in the simulation is the agent. Code any behaviours
  // that are independent of the environment here. Examples include
  // searching, resting, or moving.
  //
  // For simplicity and the smoothness of our graphical output our agents
  // will be coded to move in the world coordinates. These are floating
  // point numbers between -1.0 and 1.0, with { 0 , 0 } at the centre.
  //
  friend class Grid;
  protected:
  public:
    Grid* grid_;
    bool active;
    string label;
    float radius;
    float length;
    float growth_rate;
    float split_length;
    float inertia;
    float x;
    float y;
    float angle;
    float vel_x;
    float vel_y;
    float vel_angle;
    float force_x(0.0);
    float force_y(0.0);
    float torque(0.0);
    ABMagent* prev_(NULL);
    ABMagent* next_(NULL);

    void split()
    {
      // Need to give back a bunch of stuff for the daughter cell
      // This will not be easy to integrate...
      string daughter_label = label;
      label.append('0');
      daughter_label.append('1');
    };


    void grow(float dt)
    {
      length *= exp ( growth_rate * dt )
      if (length > split_length) {
        split();
      };
    };


    void addForce(float point_x, float point_y, float ext_force_x, float ext_force_y)
    {
      force_x += ext_force_x;
      force_y += ext_force_y;
      torque += ( point_x - x ) * ext_force_y - ( point_y - y ) * ext_force_y;
    };


    void move(float dt, float damping)
    {
      float reduced_mass = ( 2.0 * length ) / ( ( PI * radius ) + 1.0 );
      // velocities calculated with damping
      vel_x += dt * ( force_x / reduced_mass - damping * vel_x );
      vel_y += dt * ( force_y / reduced_mass - damping * vel_y );
      vel_angle += dt * ( torque / inertia - damping * vel_angle );

      // reposition cell and re-orient
      float x_old = x;
      float y_old = y;
      x += dt * vel_x;
      y += dt * vel_y;
      angle += dt * vel_angle;

      // reset forces to zero for next round
      force_x = 0.0;
      force_y = 0.0;
      torque = 0.0;

      grid_->move(this, x_old, y_old);
    };
};

ABMagent::ABMagent (
                    Grid* grid,
                    float x_,
                    float y_,
                    float radius_,
                    float length_,
                    float vel_x_,
                    float vel_y_,
                    float angle_,
                    float inertia_,
                    float vel_angle_,
                    float growth_rate_,
                    float split_length_,
                    string label_
                  )
{
  grid_ = grid
  x = 0.0;
  y = 0.0;
  radius = 1,0;
  length = 1.0;
  vel_x = 0.0;
  vel_y = 0.0;
  vel_angle = 0.0;
  torque = 0.0;
  growth_rate = 1.0;
  split_length = 2.0;
  label = "0";
  inertia = 5.0;

  grid_->add(this)
}


class Grid
{
  public:
    Grid()
    {
      // Clear the grid.
      for (int x = 0; x < NUM_CELLS_WIDTH; x++)
      {
        for (int y = 0; y < NUM_CELLS_HEIGHT; y++)
        {
          cells_[x][y] = NULL;
        }
      }
    }
    static const float CELL_SIZE = 4.0;
    static const float CHANNEL_WIDTH = 45.0;
    static const float CHANNEL_HEIGHT = 12.0;
    static const int NUM_CELLS_WIDTH = 2 + ceil( CHANNEL_WIDTH / CELL_SIZE);
    static const int NUM_CELLS_HEIGHT = 2 + ceil( CHANNEL_HEIGHT / CELL_SIZE);
  private:
    ABMagent* cells_[NUM_CELLS_HEIGHT][NUM_CELLS_WIDTH];
};


void Grid::add(ABMagent* agent)
{
  // Determine which grid cell it's in.
  int cellX = ceil(agent->x / Grid::CELL_SIZE);
  int cellY = ceil(agent->y / Grid::CELL_SIZE);

  // Add to the front of list for the cell it's in.
  agent->prev_ = NULL;
  agent->next_ = cells_[cellX][cellY];
  cells_[cellX][cellY] = agent;

  if (agent->next_ != NULL)
  {
    agent->next_->prev_ = agent;
  }
}


void Grid::move(ABMagent* agent, double x_old, double y_old)
{
  // See which cell it was in.
  int oldCellX = ceil(agent->x_old / Grid::CELL_SIZE);
  int oldCellY = ceil(agent->y_old / Grid::CELL_SIZE);

  // See which cell it's moving to.
  int cellX = ceil(agent->x / Grid::CELL_SIZE);
  int cellY = ceil(agent->y / Grid::CELL_SIZE);

  agent->x_ = x;
  agent->y_ = y;

  // If it didn't change cells, we're done.
  if (oldCellX == cellX && oldCellY == cellY) return;

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
  if (cells_[oldCellX][oldCellY] == agent)
  {
    cells_[oldCellX][oldCellY] = agent->next_;
  }

  // Add it back to the grid at its new cell.
  add(agent);
}


void Grid::handleInteractions()
{
  for (int x = 1; x < NUM_CELLS_HEIGHT; x++)
  {
    for (int y = 1; y < NUM_CELLS_WIDTH; y++)
    {
      handleCell(cells_[x][y]);
    }
  }

  for (int y = 1; y < NUM_CELLS_WIDTH; y++)
  {
    handleCellWall(cells_[1][y]);
    handleCellWall(cells_[NUM_CELLS_HEIGHT-1][y]);
    handleCellWall(cells_[NUM_CELLS_HEIGHT-2][y]);
  }
}


void Grid::handleCell(int x, int y)
{
  agent* agent = cells_[x][y];
  while (agent != NULL)
  {
    // Handle other agents in this cell.
    handleAgent(agent, agent->next_);

    // Hand agents in other cells.
    handleAgent(agent, cells_[x][y+1]);
    handleAgent(agent, cells_[x+1][y-1]);
    handleAgent(agent, cells_[x+1][y]);
    handleAgent(agent, cells_[x+1][y+1]);

    agent = agent->next_;
  }
}

void Grid::handleCellWall(int x, int y)
{
  agent* agent = cells_[x][y];
  while (agent != NULL)
  {
    // Handle other agents in this cell.
    handleAgentWall(agent, );
    agent = agent->next_;
  }
}


void Grid::handleAgent(ABMagent* agent, ABMagent* other)
{
  while (other != NULL)
  {
  float distance_x = agent.x - other.x;
  float distance_y = agent.y - other.y;
  float distance_length = sqrt( pow(distance_x, 2) + pow(distance_y, 2) );
    // Continue only if centers are a certain distance from eachother
    if ( 2 * distance_length < agent.length + other.length + agent.radius + other.radius)
    {
      // Continue only if there is potential overlap
      float project_agent;
      float project_other;
      project_agent = abs( distance_x * agent.x + distance_y * agent.y ) / distance_length;
      project_other = abs( distance_x * other.x + distance_y * other.y ) / distance_length;
      if (agent.radius + other.radius < distance_length - project_agent - project_other )
      {
        applyForces(agent, other);
      }
    }

    other = other->next_;
  }
}


void Grid::applyForces()
{

}


bool readDataFile (string fname, int layer) {
  int dimension = 1 + 2/GRIDRES;
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

int main (int argc, char **argv) {

  // Setup: set the random number generator's seed, initialize our display
  // window.
  srand (time(NULL));

  // Initialize our window.
  glutInitDisplayMode(GLUT_SINGLE);
  glutInitWindowSize(800, 800);
  glutInitWindowPosition(100, 100);
  int win = glutCreateWindow("C++ Agent Based Model");

  // Register OpenGL functional callbacks for handling keyboard inputs,
  // updating the display, and our controller function.
  glutKeyboardFunc(keyboard);
  glutDisplayFunc(drawAll);
  glutTimerFunc(25, updateTime, 0);

  // Fire up OpenGL.
  glutMainLoop();

  // Anything here won't be run until the window is closed. Want to cout some
  // stats or other information? Write a file? Do something else? put it here.
  return 0;
}
