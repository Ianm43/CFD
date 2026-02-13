#include "Incompressible-inviscid.h"

using namespace std;

int main()
{
    float t = 0;

    while( t < 1000 )
    {
        
        ExternalForce();

        Divergence( 10, 1.9 );

        Advection();

        t += dt;
    }

    return 0;
}

void ExternalForce()
{

    for( auto &row : MESH )
    {
        for( auto &CELL : row )
        {
            CELL.v -= 9.8 * dt;
        }
    }

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

    float V_new, U_new, X_vert, Y_vert, X_horz, Y_horz;
    size_t r_vert, c_vert, r_horz, c_horz;

    for( size_t r = 0; r < MESH_height - 1; ++r )
    {
        for( size_t c = 0; c < MESH_width - 1; ++c )
        {

            V_new = (MESH[r][c].v + MESH[r][c-1].v - MESH[r+1][c-1].v - MESH[r+1][c].v )/4;
            U_new = (MESH[r-1][c].u - MESH[r-1][c+1].u + MESH[r][c].v - MESH[r][c+1].v )/4;

            X_horz = -MESH[r][c].u * dt;
            Y_horz = -U_new * dt;

            X_vert = -V_new * dt;
            Y_horz = -MESH[r][c].v * dt;

            // backtracked cells
            c_horz = ceil(( X_horz ) / MESH_width);
            r_horz = ceil(( Y_horz ) / MESH_height);

            r_vert = ceil(( X_vert ) / MESH_height);
            c_vert = ceil(( Y_vert ) / MESH_width);

            //weighted velocities
            

        }
    }

}

