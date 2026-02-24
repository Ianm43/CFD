#include "inco.h"

using namespace std;

int main()
{
    INCO_opts solver_opts;
    Graphics_opts graph_opts;

    bool DRAW_ELLIPSE = 1;
    bool INITIAL = true;

    double FLOW_VEL_u = 1;
    double FLOW_VEL_v = 0;

    solver_opts.MESH_HEIGHT = 200;
    solver_opts.MESH_WIDTH = 200;
    solver_opts.CELL_SIZE = 0.005;

    solver_opts.ITERATIONS = 40;
    solver_opts.REPORT_INTERVAL = 1;
    solver_opts.TIMESTEPS = 5; 
    solver_opts.OVER_RELAXATION = 1.9;

    solver_opts.TIMESTEP = solver_opts.CELL_SIZE/sqrt(FLOW_VEL_u*FLOW_VEL_u+FLOW_VEL_v*FLOW_VEL_v)*0.25;

    Profile P;
    
    Force Lift_Report( 0, 1, "Wing", "Lift", 1, 0 );

    std::vector< std::vector< cell > > MESH( solver_opts.MESH_HEIGHT, std::vector<cell>(solver_opts.MESH_WIDTH) );;

    cout << "Starting " << to_string( solver_opts.MESH_WIDTH ) << " x " << to_string( solver_opts.MESH_HEIGHT ) << " Simulation with:\n" 
         << "Timestep:   " << to_string( solver_opts.TIMESTEP ) << "\n"
         << "Iterations: " <<to_string( solver_opts.ITERATIONS ) << "\n"
         << "Cell Size:  " << to_string( solver_opts.CELL_SIZE ) << "\n"
         << "For: " << to_string( solver_opts.TIMESTEPS ) << " steps, and reporting every: " << to_string( solver_opts.REPORT_INTERVAL ) << " steps" << endl;


    // super dumb way of making a skewed ellipse 
    //for help defining varibles: https://www.desmos.com/calculator/fgdc56e9nm
        float X, Y, R, SLOPE, HORZ_SHIFT, VERT_SHIFT, VERT_SQUISH;
        SLOPE = 0;
        VERT_SHIFT = solver_opts.MESH_HEIGHT/2.0;
        HORZ_SHIFT = solver_opts.MESH_WIDTH/2.0;
        VERT_SQUISH = 1;
        R = solver_opts.MESH_HEIGHT*solver_opts.MESH_HEIGHT / 100; //semi major axis

    for( size_t r = 0; r < solver_opts.MESH_HEIGHT; ++r )
    {
        //Left side treatment
        MESH[r][0].Fluid = 0;
        MESH[r][1].u = FLOW_VEL_u;
        MESH[r][0].u = FLOW_VEL_u;
        MESH[r][0].v = FLOW_VEL_v;
        MESH[r][1].v = FLOW_VEL_v;

        //right side treatment
        MESH[r][solver_opts.MESH_WIDTH-1].Fluid = 0;
        MESH[r][solver_opts.MESH_WIDTH-1].u = FLOW_VEL_u;
        MESH[r][solver_opts.MESH_WIDTH-1].v = FLOW_VEL_v;
    

        
        for( size_t c = 0; c < solver_opts.MESH_WIDTH; ++c )
        {

            if( c <= 200   &&  r > solver_opts.MESH_HEIGHT/2 - 20 && r < solver_opts.MESH_HEIGHT/2 + 20 )
            {
                //MESH[r][c].density = 0.5;
            }
            MESH[r][c].x = c * solver_opts.CELL_SIZE;
            MESH[r][c].y = r * solver_opts.CELL_SIZE;
            Y = r + SLOPE*( c - HORZ_SHIFT ) - VERT_SHIFT ;

            X = c + SLOPE*( c - HORZ_SHIFT ) - HORZ_SHIFT;

            if( DRAW_ELLIPSE && ( (Y)*(Y)*VERT_SQUISH + (X)*(X) ) <= (R) )
            {
                
                MESH[r][c].u = 0;
                MESH[r][c].v = 0;
                MESH[r][c].Fluid = 0;
                MESH[r][c].Tag = "Wing";
                continue;
            }
            

            MESH[0][c].Fluid = 0; // bottom treatment
            MESH[0][c].u = 0;//FLOW_VEL_u; 
            MESH[0][c].v = 0;//FLOW_VEL_v;
            MESH[solver_opts.MESH_HEIGHT-1][c].Fluid = 0; // top treatment
            MESH[solver_opts.MESH_HEIGHT-1][c].u = 0;//FLOW_VEL_u; 
            MESH[solver_opts.MESH_HEIGHT-1][c].v = 0; //FLOW_VEL_v;
        }

    }
    
    // initial solution "best" guess
    for( size_t r = 1; r < solver_opts.MESH_HEIGHT - 1 && INITIAL; ++r )
    {
        for( size_t c = 1; c < solver_opts.MESH_WIDTH - 1; ++c )
        {
            // MAKE SURE THAT SOLID EDGES DONT GET ASSIGNED VELOCITIES THAT NEVER CHANGE
            if( MESH[r][c-1].Fluid && MESH[r][c].Fluid )
            {
                MESH[r][c].u = FLOW_VEL_u;
                MESH[r][c].v = FLOW_VEL_v;
            }
        }
    }



    // start with a converged solution
    //Divergence( 1000, solver_opts.OVER_RELAXATION );

    std::vector<cell> NEW_MESH( solver_opts.MESH_HEIGHT*solver_opts.MESH_WIDTH );

    for( std::vector<cell> &ROW : MESH )
    {
        for( cell &CELL : ROW )
        {
            NEW_MESH.insert( NEW_MESH.begin(), CELL );
        }
    }
    


    INCO_SOLVER solver( solver_opts, NEW_MESH ); 

    Graphics graph( solver );

    for( size_t step = 0; step < solver_opts.TIMESTEPS; ++step )
    {

        solver.Solve();

        if( step % solver_opts.REPORT_INTERVAL == 0 )
            graph.PrintBmp();

    }
    /* 

    std::string name;
    size_t t = 0;


    while( t <= TIMESTEPS )
    {
        //ExternalForce();  // fuck gravity
        P.Start();
        Divergence( ITERATIONS, OVER_RELAXATION );
        P.End();
        cout << "Divergence step took: " << to_string( P.Elapsed_Seconds() ) << endl;

        P.Start();
        Advection(); //aprox material derivative 
        P.End();
        cout << "Advection step took: " << to_string( P.Elapsed_Seconds() ) << endl;
        
        if( t % REPORT_INTERVAL == 0 )
        {
            //Lift_Report.REPORT( MESH, CELL_SIZE );
            name = "Solution_Data_1/frame_" + to_string(t) +".bmp";

            Bitmap( MESH, name );
        }

        ++t;
    }
 */
    return 0;
}

/* 
void ExternalForce()
{
    for( size_t r = 1; r < MESH_height; ++r )
    {
        for( size_t c = 1; c < MESH_width; ++c )
        {
            //update velocites unless those velocies belong to a non fluid cell 
            MESH[r][c].v += dt * GRAVITY * ( MESH[r][c].Fluid * MESH[r-1][c].Fluid );
        }
    }
}

void Divergence( size_t iterations, double Relaxation )
{
    //Profile P;
    size_t r,c;
    #pragma omp parallel
    {
        for( size_t i = 0; i <= iterations; ++i )
        {
            
            //for( size_t row = 1; row < MESH_height - 1; ++row )
            //{
            //    for( size_t col = 1; col < MESH_width - 1; ++col )
            //    {
            //
            //        CellDivergence( row, col, Relaxation );
            //        
           //    }
            //} 

            #pragma omp for
           for( size_t index = 0; index < (MESH_height-1)*(MESH_width-1); index+=2 )
           {
                r = index / (MESH_width-1);
                c = index % (MESH_width-1);
                assert( r<MESH_height && c < MESH_width );
                CellDivergence( index/(MESH_width-1), index%(MESH_width-1), Relaxation );

           }
            #pragma omp barrier

            #pragma omp for
            for( size_t index = 1; index < (MESH_height-1)*(MESH_width-1); index+=2 )
            {

                CellDivergence( index/(MESH_width-1), index%(MESH_width-1), Relaxation );

            }

            #pragma omp barrier
        }
    }
}

void CellDivergence( const size_t &row, const size_t &col, const double &Relaxation )
{
    //skip non fluid cells
    if ( !MESH[row][col].Fluid ){ return; }

    uint8_t neighbors = 0;
    double Div = 0;

    neighbors = MESH[row - 1][col].Fluid + MESH[row + 1][col].Fluid + MESH[row][col + 1].Fluid + MESH[row][col - 1].Fluid;
    Div = (-MESH[row][col].u - MESH[row][col].v + MESH[row + 1][col].v + MESH[row][col + 1].u);

    // store the divergence values for visualization
    MESH[row][col].divergence = Div;

    Div *= Relaxation;

    // make the divergence of this cell 0 while leaving out boundary edges
    MESH[row][col].u += Div * MESH[row][col - 1].Fluid / neighbors;
    MESH[row][col].v += Div * MESH[row - 1][col].Fluid / neighbors;
    MESH[row + 1][col].v -= Div * MESH[row + 1][col].Fluid / neighbors;
    MESH[row][col + 1].u -= Div * MESH[row][col + 1].Fluid / neighbors;

    // estimate kinetic pressure
    MESH[row][col].p -= Div / neighbors * density * CELL_SIZE / dt;
}

void Advection()
{
    TEMP = MESH;

    double X_vert, Y_vert, X_horz, Y_horz;
    size_t r_vert, c_vert, r_horz, c_horz, r_den, c_den;

    for( size_t r = 0; r < MESH_height; ++r )
    {
        for( size_t c = 0; c < MESH_width; ++c )
        {
            //skip non fluid cells
            if( !MESH[r][c].Fluid ){ continue; }

            //get average velocities around each edge and look back

                X_vert = MESH[r][c].x - ( (MESH[r][c].u + MESH[r-1][c].u + MESH[r][c+1].u + MESH[r-1][c+1].u) / 4.0 ) * dt;
                X_vert *= ( X_vert > 0 );
                Y_vert = MESH[r][c].y - MESH[r][c].v * dt;
                Y_vert *= ( Y_vert > 0 );

                X_horz = MESH[r][c].x - MESH[r][c].u * dt;
                X_horz *= ( X_horz > 0 );
                Y_horz = MESH[r][c].y - ( (MESH[r][c].v + MESH[r][c-1].v + MESH[r+1][c].v + MESH[r+1][c-1].v) / 4.0 ) * dt;
                Y_horz *= ( Y_horz > 0 );


            // predict last location of "particle"

                c_vert = ( X_vert / CELL_SIZE );
                if( c_vert > MESH_width - 2 ){ c_vert = MESH_width - 2; }

                r_vert =  ( Y_vert / CELL_SIZE );
                if( r_vert > MESH_height - 2 ){ r_vert = MESH_height - 2; }

                c_horz =  ( X_horz / CELL_SIZE );
                if( c_horz > MESH_width - 2 ){ c_horz = MESH_width - 2; }

                r_horz = ( Y_horz / CELL_SIZE );
                if( r_horz > MESH_height - 2 ){ r_horz = MESH_height - 2; }

                c_den = ( X_horz / CELL_SIZE );
                if( c_den > MESH_width - 2 ){ c_den = MESH_width - 2; }

                r_den = ( Y_vert / CELL_SIZE );
                if( r_den > MESH_height - 2 ){ r_den = MESH_height - 2; }

            if( c_vert < 0 || c_vert > MESH_width - 2 )
            {
                cout << "c_vert is out of range: " << to_string(c_vert) << endl;
                cout << "X_vert: " << to_string(X_vert);
            }
            
            assert( c_vert <= MESH_width - 2 && c_vert >= 0);

            // get average velocities around position

            TEMP[r][c].u = GetWeightedCellAverage( r_horz, c_horz, ( X_horz - MESH[r_horz][c_horz].x ) / CELL_SIZE, ( Y_horz - MESH[r_horz][c_horz].y ) / CELL_SIZE ).u;

            TEMP[r][c].v = GetWeightedCellAverage( r_vert, c_vert, ( X_vert - MESH[r_vert][c_vert].x ) / CELL_SIZE, ( Y_vert - MESH[r_vert][c_vert].y ) / CELL_SIZE ).v;

            TEMP[r][c].density = GetWeightedCellAverage( r_den, c_den, ( X_horz - MESH[r_den][c_den].x ) / CELL_SIZE, ( Y_vert - MESH[r_den][c_den].y ) / CELL_SIZE).density;

        }
    }
    MESH = TEMP;
}

cell GetWeightedCellAverage( size_t r, size_t c, double x, double y )// x and y are assumed to be between 0 and 1
{
    cell AveragedCell;

    double BL = (1 - x)*(1-y);
    double BR = (x)*(1-y);
    double TL = (1-x)*(y);
    double TR = x*y;

    AveragedCell.v = ( BL*MESH[r][c].v + TL*MESH[r+1][c].v + BR*MESH[r][c+1].v + TR*MESH[r+1][c+1].v );
    AveragedCell.u = ( BL*MESH[r][c].u + TL*MESH[r+1][c].u + BR*MESH[r][c+1].u + TR*MESH[r+1][c+1].u );
    AveragedCell.density = ( BL*MESH[r][c].density + TL*MESH[r+1][c].density + BR*MESH[r][c+1].density + TR*MESH[r+1][c+1].density );

    return AveragedCell;

}

void Bitmap( vector< std::vector< cell > > & Img_Data, std::string &filename )
{
    cell MAX = MAX_CELL();
    cell MIN = MIN_CELL();
    
    MAX.u -= MIN.u * (MIN.u < 0);
    MAX.v -= MIN.v * (MIN.v < 0);
    MAX.p -= MIN.p * (MIN.p < 0);
    MAX.density -= MIN.density * (MIN.density < 0) ;
    MAX.divergence -= MIN.divergence * (MIN.divergence < 0);
    
    assert( MAX.divergence < 1 && "Divergence should not be this high" );

    double BIGGEST_VEL = sqrt(MAX.u*MAX.u + MAX.v*MAX.v);
    double BIGGEST = ( MAX.p * PRINT_P + MAX.u * PRINT_U + MAX.v * PRINT_V + MAX.divergence * PRINT_DIV + BIGGEST_VEL * PRINT_VEL + PRINT_DEN );
    double VEL = 0;

    double PRINT_VAL;

    BmpHeader header;
    IBmpInfoHeader infoHeader;

    infoHeader.Height = Img_Data.size();
    infoHeader.Width = Img_Data[0].size(); 
    infoHeader.horizontalResolution = 2 * abs(infoHeader.Width); // fits the image to 50 cm
    infoHeader.VerticalResolution = 2 * abs(infoHeader.Height);

    header.bfSize = 54 + abs(infoHeader.Width * infoHeader.Height * 3);

    ofstream file( filename, ios::binary);

    header.HeaderWrite( file );
    infoHeader.InfoHeaderWrite( file );

    Pixel pixel;
    
    for( auto const &ROW : Img_Data )
    {
        for( cell const &VALUE : ROW )
        {
            // Assign color based on tile value

            VEL = sqrt( VALUE.u*VALUE.u + VALUE.v*VALUE.v );
            PRINT_VAL = ( (VALUE.p-MIN.p*(MIN.p<0)) * PRINT_P + (VALUE.u-MIN.u*(MIN.p<0)) * PRINT_U + (VALUE.v-MIN.v*(MIN.v<0)) * PRINT_V + (VALUE.divergence-MIN.divergence*(MIN.divergence<0)) * PRINT_DIV + VEL * PRINT_VEL )/BIGGEST;

            if( (PRINT_VAL * 255) > 255 || (PRINT_VAL * 255) < 0 )
            {
                cout << "PRINT_VAL: " + to_string(PRINT_VAL) << endl;
            }

            if( !PRINT_DEN )
            {
                pixel.Red   = (uint8_t)( ( PRINT_VAL * 255 ) );
                pixel.Green = (uint8_t)( ( PRINT_VAL * 0 )  );
                pixel.Blue  = (uint8_t)( ( 255 - PRINT_VAL * 255 ) );
            }else
            {
                pixel.Red   = (uint8_t)( ( VALUE.density * 255 ) );
                pixel.Green = (uint8_t)( ( VALUE.density * 255 ) );
                pixel.Blue  = (uint8_t)( ( VALUE.density * 255 ) );
            }

            if( !VALUE.Fluid )
            {
                pixel.Red = 0; pixel.Blue = 0; pixel.Green = 0;
            }


            pixel.PixelWrite( file );

        }
        //extra "padding" cuz rows must be a multiple of 4 (WHYYY?)
        file.write( (char *)0, infoHeader.Width % 4 ); // pro tip width != Height
            
    }

    cout << "Saving: " << filename <<";  divergence: " << to_string( MAX.divergence ) << endl;
    
    file.close();    
}

cell MIN_CELL()
{
    cell min;

    for( auto ROW : MESH )
    {
        for( auto CELL : ROW )
        {
            if( CELL.p < min.p )
                min.p = CELL.p;
            if( CELL.u < min.u )
                min.u = CELL.u;
            if( CELL.v < min.v )
                min.v = CELL.v;
            if( CELL.density < min.density )
                min.density = CELL.density;
            if( CELL.divergence < min.divergence )
                min.divergence = CELL.divergence;
        }
    }
    return min;
}

cell MAX_CELL()
{
    cell max;
    for( const vector<cell> &ROW : MESH )
    {
        for( const cell& CELL: ROW )
        {
            if( CELL.p > max.p )
                max.p = CELL.p;
            if( CELL.u > max.u )
                max.u = CELL.u;
            if( CELL.v > max.v )
                max.v = CELL.v;
            if( CELL.density > max.density )
                max.density = CELL.density;
            if( CELL.divergence > max.divergence )
                max.divergence = CELL.divergence;
        }
    }
    return max;
}
 */