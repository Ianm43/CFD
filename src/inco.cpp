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


    //total mesh size (Height*Width) should be an odd number :)
    solver_opts.MESH_HEIGHT = 400;
    solver_opts.MESH_WIDTH = 1000;
    solver_opts.CELL_SIZE = 0.05;

    solver_opts.ITERATIONS = 100;
    solver_opts.OVER_RELAXATION = 1.9;
    solver_opts.CFL = 0.5;
    solver_opts.GRAVITY = -2;


    
    INCO_SOLVER solver( solver_opts, mesh_options );

    solver.Make_MESH();

    solver.Initilaize();

    Graphics graph( solver, graph_opts );


    while( solver.GetStep() <= 3000 )
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