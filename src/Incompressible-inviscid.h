#ifndef INCOMPRESSIBLE_INVISCID
#define INCOMPRESSIBLE_INVISCID

#include <stdint.h>

struct cell
{  
    float u = 0; // x velocity
    float v = 0; // y velocity
    float density = 1;
};

static const size_t MESH_width = 100;
static const size_t MESH_height = 100;
static cell MESH[MESH_width][MESH_height];

static float dt;

int main();
void ExternalForce();
void Divergence( uint8_t iterations, float Relaxation );
void Advection();

#endif