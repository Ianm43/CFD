#include "Incompressible-inviscid.h"

using namespace std;

int main()
{

    



    return 0;
}

void ExternalForce()
{



}

void Divergence( uint8_t iterations, float Relaxation )
{

    float Div = 0;

    for( uint8_t i = 0; i <= iterations; ++i )
    {

        for( size_t row = 0; row < MESH_height - 1; ++row )
        {
            for( size_t col = 0; col < MESH_width - 1; ++col  )
            {

                Div = (MESH[row][col].u + MESH[row][col].v - MESH[row+1][col].v - MESH[row][col+1].u) * Relaxation;

                if( row == 0 ) // bottom side treatment
                {

                    MESH[row][col].u -= Div/3.0;
                    MESH[row][col].v = 0;
                    MESH[row+1][col].v -= Div/3.0;
                    MESH[row][col+1].u -= Div/3.0;

                    continue;
                }

                if( row == MESH_height - 1 ) //Top side treatment
                {
                    MESH[row+1][col].v = 0; //wall condition

                    MESH[row][col].u -= Div/3.0;
                    MESH[row][col].v -= Div/3.0;
                    MESH[row][col+1].u -= Div/3.0;

                    continue;
                }

                if( col == 0 ) // left side treatment
                {

                    MESH[row][col].u = 5; // inflow condition

                    MESH[row][col].v -= Div/3.0;
                    MESH[row+1][col].v -= Div/3.0;
                    MESH[row][col+1].u -= Div/3.0;


                    continue;
                }

                if( col == MESH_width - 1 ) //right side treatment
                {
                    MESH[row][col+1].u = 5; // outflow condition

                    MESH[row][col].u -= Div/3.0;
                    MESH[row][col].v -= Div/3.0;
                    MESH[row+1][col].v -= Div/3.0;

                    continue;
                }

                

                MESH[row][col].u -= Div/4.0;
                MESH[row][col].v -= Div/4.0;
                MESH[row+1][col].v -= Div/4.0;
                MESH[row][col+1].u -= Div/4.0;
                
            }
        }

    }

}

void Advection()
{

    // estimate the previous position of the velocity on each edge, then update the new velocity on this edge
    // with the velocity at that previous position

    for( size_t row = 0; row < MESH_height - 1; ++row )
    {
        for( size_t col = 0; col < MESH_width - 1; ++col )
        {



        }
    }

}

