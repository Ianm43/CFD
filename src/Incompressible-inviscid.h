#ifndef INCOMPRESSIBLE_INVISCID
#define INCOMPRESSIBLE_INVISCID

#include <stdint.h>

struct cell
{  
    float u = 0;
    float v = 0;
    float density = 1;
};

static cell MESH[100][100];

static float dt;

int main();
void ExternalForce();
void Divergence( uint8_t iterations );
void Advection();

#endif