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
    graph_opts.PRINT_V = 0;
    graph_opts.PRINT_DEN = 1;
    graph_opts.PRINT_DIV = 0;
    graph_opts.PRINT_CURL = 0;

    graph_opts.PixelsPerCell = 2;

    //total mesh size (Height*Width) should be an odd number :)
    solver_opts.MESH_HEIGHT = 400;
    solver_opts.MESH_WIDTH = 1000;
    solver_opts.CELL_SIZE = 0.005;

    mesh_options.DRAW_ELLIPSE = true;
    
    mesh_options.FLOW_VEL_u = 1;
    mesh_options.FLOW_VEL_v = 0;

    solver_opts.ITERATIONS = 30;
    solver_opts.OVER_RELAXATION = 1.9;
    solver_opts.CFL = 1;
    solver_opts.GRAVITY = 0;
    // this setting is broken
    solver_opts.VC = 0;


    
    INCO_SOLVER solver( solver_opts, mesh_options );

    Graphics graph( solver, graph_opts );

    solver.Make_MESH();

    solver.Initilaize();
    
    graph.PrintBmp();

    
    while( solver.GetStep() <= 10000 )
    {
        solver.Solve();
        if( solver.GetStep() % 100 == 0 )
        {
            solver.CalculateCurl();
            graph.PrintBmp();
        }
    }
    
    
    return 0;
}