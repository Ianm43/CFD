#ifndef INCOMPRESSIBLE_INVISCID
#define INCOMPRESSIBLE_INVISCID

#include <stdint.h>
#include <cmath>
#include "Bitmap.h"

struct cell
{  
    float u = 0; // x velocity in m/s
    float v = 0; // y velocity in m/s
    float density = 1;
};

static const float dt = 0.1; // time step in sec
static const float CELL_SIZE = 0.1; // size of each cell in m
static const size_t MESH_width = 100;
static const size_t MESH_height = 100;
static cell MESH[MESH_width][MESH_height];



int main();
void ExternalForce();
void Divergence( uint8_t iterations, float Relaxation );
void Advection();

#endif