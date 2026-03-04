#include "inco.h"

using namespace std;

int main()
{
    INCO_opts solver_opts;
    MESH_opts mesh_options;
    Graphics_opts graph_opts;

    graph_opts.PRINT_VEL = 0;
    graph_opts.PRINT_P = 0;
    graph_opts.PRINT_U = 0;
    graph_opts.PRINT_V = 1;
    graph_opts.PRINT_DEN = 0;
    graph_opts.PRINT_DIV = 0;

    double FLOW_VEL_u = 1;
    double FLOW_VEL_v = 0;

    //total mesh size (Height*Width) should be an odd number :)
    solver_opts.MESH_HEIGHT = 401;
    solver_opts.MESH_WIDTH = 1000;
    solver_opts.CELL_SIZE = 1;

    solver_opts.ITERATIONS = 40;
    solver_opts.REPORT_INTERVAL = 500;
    solver_opts.TIMESTEPS = 2001; 
    solver_opts.OVER_RELAXATION = 1.9;

    solver_opts.TIMESTEP = solver_opts.CELL_SIZE/sqrt(FLOW_VEL_u*FLOW_VEL_u+FLOW_VEL_v*FLOW_VEL_v)*0.25;


    // make it easier to populate mesh values by making row column indecies 
    
    INCO_SOLVER solver( solver_opts, mesh_options );

    solver.Make_MESH();

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