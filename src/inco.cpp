#include "inco.h"

using namespace std;

int main()
{
    INCO_opts solver_opts;
    Graphics_opts graph_opts;

    graph_opts.PRINT_VEL = 1;
    graph_opts.PRINT_P = 0;
    graph_opts.PRINT_U = 0;
    graph_opts.PRINT_V = 0;
    graph_opts.PRINT_DEN = 0;
    graph_opts.PRINT_DIV = 0;

    bool DRAW_ELLIPSE = 1;

    double FLOW_VEL_u = 1;
    double FLOW_VEL_v = 0;

    //total mesh size (Height*Width) should be an odd number :)
    solver_opts.MESH_HEIGHT = 401;
    solver_opts.MESH_WIDTH = 1000;
    solver_opts.CELL_SIZE = 1;

    solver_opts.ITERATIONS = 40;
    solver_opts.REPORT_INTERVAL = 1;
    solver_opts.TIMESTEPS = 0; 
    solver_opts.OVER_RELAXATION = 1.9;

    solver_opts.TIMESTEP = solver_opts.CELL_SIZE/sqrt(FLOW_VEL_u*FLOW_VEL_u+FLOW_VEL_v*FLOW_VEL_v)*0.25;

    Profile P;
    
    Force Lift_Report( 0, 1, "Wing", "Lift", 1, 0 );



    // start with a converged solution
    //Divergence( 1000, solver_opts.OVER_RELAXATION );

    std::vector<cell> MESH( solver_opts.MESH_HEIGHT*solver_opts.MESH_WIDTH );

    // make it easier to populate mesh values by making row column indecies 
    size_t row;
    size_t col;

    // super dumb way of making a skewed ellipse 
    //for help defining varibles: https://www.desmos.com/calculator/fgdc56e9nm
        float X, Y, R, SLOPE, HORZ_SHIFT, VERT_SHIFT, VERT_SQUISH;
        SLOPE = 0;
        VERT_SHIFT = solver_opts.MESH_HEIGHT/2;
        HORZ_SHIFT = solver_opts.MESH_WIDTH/5;
        VERT_SQUISH = 1;
        R = (float)solver_opts.MESH_HEIGHT*solver_opts.MESH_HEIGHT / 100; //semi major axis

    for( size_t index = 0; index < MESH.size(); ++index )
    {
        col = index % solver_opts.MESH_WIDTH;
        row = index / solver_opts.MESH_WIDTH;

        MESH[index].x = index % solver_opts.MESH_WIDTH * solver_opts.CELL_SIZE;
        MESH[index].y = index / solver_opts.MESH_WIDTH * solver_opts.CELL_SIZE;

        //bottom side
        if( row == 0 )
        {
            MESH[index].Fluid = 0;
            MESH[index].u = FLOW_VEL_u;
            MESH[index].v = FLOW_VEL_v;
            MESH[index+solver_opts.MESH_WIDTH].v = FLOW_VEL_v;
        }

        //top side
        if( row == solver_opts.MESH_HEIGHT - 1 )
        {
            MESH[index].Fluid = 0;
            MESH[index].u = FLOW_VEL_u;
            MESH[index].v = FLOW_VEL_v;
        }

        // left side
        if( col == 0 )
        {
            MESH[index].Fluid = 0;
            MESH[index].u = FLOW_VEL_u;
            MESH[index+1].u = FLOW_VEL_u;
            MESH[index].v = FLOW_VEL_v;
            
        }

        // right side
        if( col == solver_opts.MESH_WIDTH - 1 )
        {
            MESH[index].Fluid = 0;
            MESH[index].u = FLOW_VEL_u;
            MESH[index].v = FLOW_VEL_v;
        }

        if( (row > (solver_opts.MESH_HEIGHT/2 - 20)) && (row < (solver_opts.MESH_HEIGHT/2 + 20))  && (col == 0) )
        {

            MESH[index].density = 0.5;

        }

        Y = row + SLOPE * (col - HORZ_SHIFT) - VERT_SHIFT;

        X = col + SLOPE * (col - HORZ_SHIFT) - HORZ_SHIFT;

        if (DRAW_ELLIPSE && round((Y) * (Y)*VERT_SQUISH + (X) * (X)) <= (R))
        {
            MESH[index].u = 0;
            MESH[index].v = 0;
            MESH[index].Fluid = 0;
        }
    }

    INCO_SOLVER solver( solver_opts, MESH );

    solver.Initilaize();

    Graphics graph( solver, graph_opts );

    graph.PrintBmp();

    for( size_t step = 0; step < solver_opts.TIMESTEPS; ++step )
    {

        solver.Solve();

        if( step % solver_opts.REPORT_INTERVAL == 0 )
            graph.PrintBmp();

    }
    
    return 0;
}