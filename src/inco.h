#ifndef INCO_H
#define INCO_H

#include <cmath>

#include <cstdint>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <chrono>


int main();

struct Profile
{
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_time;

    void Start()
    {
        start_time = std::chrono::high_resolution_clock::now();
    }

    void End()
    {
        end_time = std::chrono::high_resolution_clock::now();
    }

    double Elapsed_Seconds()
    {
        return std::chrono::duration<double>(end_time - start_time).count();
    }
};

//random vector structure for stuff :)
struct vec
{
    double x;
    double y;
};

// fundamental MESH unit
struct cell
{
    double x = 0; // x coordinate of the bottom left vertex
    double y = 0;
    double u = 0;       // x velocity in m/s
    double v = 0;       // y velocity in m/s
    double density = 0; // SMOKE DENSITY NOT PHYSICAL DENSITY
    double divergence = 0;
    double time_residual = 0;
    double curl = 0;
    double p = 0;   // estimated kinetic pressure for visulaization
    bool Fluid = 1; // pretty self explanatory
    std::string Tag = "Fluid"; // specific tag for Force Report
};

//option structs, maybe not the best way to consolidate settings, but MATLAB does it

    struct INCO_opts
    {
        size_t ITERATIONS = 40;
        double CELL_SIZE = 0.01;
        size_t MESH_WIDTH = 100;
        size_t MESH_HEIGHT = 100;
        double OVER_RELAXATION = 1.9;
        double GRAVITY = 0;
        double REF_density = 1.225;
        double CFL = 0.1;
        double VC = 0.01;
    };

    struct Graphics_opts
    {
        bool PRINT_P = 1;
        bool PRINT_V = 0;
        bool PRINT_U = 0;
        bool PRINT_DEN = 0;
        bool PRINT_DIV = 0;
        bool PRINT_VEL = 0;
        bool PRINT_CURL = 0;
        uint8_t PixelsPerCell = 1;
    };

    struct MESH_opts
    {

        bool DRAW_ELLIPSE = true;
        double SLOPE = 0;
        double VERT_SQUISH = 1;
        double FLOW_VEL_u = 1;
        double FLOW_VEL_v = 0;

    };

//structure for extracting and displaying force due to pressure
// I really don't think this is accurate since im not actually
// calculating pressure but oh well.
struct Force
{
    std::string Tag;
    double x_norm;
    double y_norm;
    std::string Name;
    bool _PrintToFile = 0;
    bool _PrintToTerm = 0;

    Force(double x, double y, std::string CELL_TAG, std::string name, bool term, bool file) : 
    Tag(CELL_TAG), x_norm(x), y_norm(y), Name(name), _PrintToFile(file), _PrintToTerm(term) {};

    double REPORT(const std::vector<std::vector<cell>> &MESH, double CELL_SIZE);
};

// extras copy and pasted from another project to print bmp images.
    struct BmpHeader // 14 bytes
    {
        char bfType[2] = {'B', 'M'};
        uint32_t bfSize;
        uint32_t ReservedBytes = 0;
        uint32_t Headeroffset = 54;

        void HeaderWrite(std::ofstream &file)
        {
            file.write((char *)&this->bfType, 2);
            file.write((char *)&this->bfSize, sizeof(uint32_t));
            file.write((char *)&this->ReservedBytes, sizeof(uint32_t));
            file.write((char *)&this->Headeroffset, sizeof(uint32_t));
        }
    };

    struct IBmpInfoHeader // 40 bytes
    {
        uint32_t Size = 40;
        int32_t Width;  // pixels
        int32_t Height; // pixels
        uint16_t Planes = 1;
        uint16_t colorDepth = 24;
        uint32_t Compression = 0;
        uint32_t SizeImage = 0;
        int32_t horizontalResolution; // pixels per meter
        int32_t VerticalResolution;   // pixels per meter
        uint32_t ClrTable = 0;
        uint32_t ClrImportant = 0;

        void InfoHeaderWrite(std::ofstream &file);
    };

    struct Pixel
    {
        uint8_t Red = 0;
        uint8_t Blue = 0;
        uint8_t Green = 0;

        void PixelWrite(std::ofstream &file);
    };


//main solver class
class INCO_SOLVER
{
    private:
        INCO_opts _opts;
        MESH_opts _mesh_opts;
        std::vector<cell> _MESH;
        std::vector<cell> _TEMP;
        double dt = 0.01;
        std::size_t _Step = 0; // curent time step
        double h;

        void ExternalForces();
        void CELLDIVERGENCE(std::size_t index);
        void ODDCELLS();
        void EVENCELLS();
        void Divergence();
        cell InterpolateCell( double x, double y);
        void StableTimeStep();
        vec Advect_Step_Back( double x, double y, double dx, double dy );
        void Advect();


        double dudy( double x, double y, double delta );
        double dvdx( double x, double y, double delta );
        double Curl( double x, double y );

        std::size_t CoordToCell(double x, double y);

    public:
        
        INCO_SOLVER(INCO_opts Options, MESH_opts M_Options ) : 
        _opts(Options), _mesh_opts( M_Options ), _MESH( _opts.MESH_HEIGHT*_opts.MESH_WIDTH ), _TEMP(_MESH), h( _opts.CELL_SIZE / 2 ) {};
        void Solve();
        void Make_MESH();
        void Initilaize();
        void CalculateCurl();
        const cell GetMax();
        const cell GetMin();
        const cell &GetCell(size_t index);
        std::vector<cell> InterpolateMesh( size_t Pix_Per_Cell );
        

        const size_t &GetMeshHeight() { return _opts.MESH_HEIGHT; }
        const size_t &GetMeshWidth() { return _opts.MESH_WIDTH; }
        const size_t &GetStep() { return _Step; }
};

// seperate graphics class, just prints a bitmap for the moment, but the plan
//is to integrate this with a GUI framework at some point.
class Graphics
{
    Graphics_opts _opts;
    INCO_SOLVER &_SOLVER;

public:
    Graphics( INCO_SOLVER &solver, Graphics_opts &opts ): _opts(opts), _SOLVER( solver ){};
    void PrintBmp();
};

//INCO_SOLVER METHODS: 

//main solver loop
void INCO_SOLVER::Solve()
{
    _Step++;

    StableTimeStep();

    ExternalForces();
        
    Divergence();
        
    Advect();       
}

//returns the index of the cell which contains the point (x,y)
std::size_t INCO_SOLVER::CoordToCell(double x, double y)
{
    // round negative values to 0
    x *= (x >= 0);
    y *= (y >= 0);

    // round values greater than mesh dimensions to the edge
    if (x > ( _opts.MESH_WIDTH ) * _opts.CELL_SIZE)
    {
        x = ( _opts.MESH_WIDTH ) * _opts.CELL_SIZE;
    }
    if (y > ( _opts.MESH_HEIGHT ) * _opts.CELL_SIZE)
    {
        y = ( _opts.MESH_HEIGHT ) * _opts.CELL_SIZE;
    }
    
    // convert from x,y coordinates to the cell that contains those coordinates
    // increase each coordinate by one thousandth of half a cell because floating point error...
    // there's got to be a better way to do this, but oh well

    return floor( (y+h/1000) / _opts.CELL_SIZE ) * _opts.MESH_WIDTH + floor((x+h/1000) / _opts.CELL_SIZE);
}

//basically let the divergence converge before we add the time steping
// maybe not necissary, or could be replaced with some low order scheme 
// to predict an initial solution based on boundary vals. 
void INCO_SOLVER::Initilaize()
{
    Divergence();
    for( size_t i = 0; i < 10000 && std::abs(GetMax().divergence) > 0.02; i += _opts.ITERATIONS )
        Divergence();

}

//populate the _MESH vector with a rectangular mesh
// this might be the biggest limitation of this solver, not even structured mesh
// excuseively rectangular cells stacked up next to each other.
void INCO_SOLVER::Make_MESH()
{
    size_t row;
    size_t col;

    double X, Y;

    // super dumb way of making a skewed ellipse 
    //for help defining varibles: https://www.desmos.com/calculator/fgdc56e9nm

    double R = (float)_opts.MESH_HEIGHT * _opts.MESH_HEIGHT / 200; // semi major axis
    double VERT_SHIFT = _opts.MESH_HEIGHT / 2;
    double HORZ_SHIFT = _opts.MESH_WIDTH / 5;

    for( size_t index = 0; index < _MESH.size(); ++index )
    {
        col = index % _opts.MESH_WIDTH;
        row = index / _opts.MESH_WIDTH;

       _MESH[index].x = col * _opts.CELL_SIZE;
       _MESH[index].y = row * _opts.CELL_SIZE;

        //bottom side
        if( row == 0 )
        {
           _MESH[index].Fluid = 0;
           _MESH[index].u = _mesh_opts.FLOW_VEL_u;
           _MESH[index].v = _mesh_opts.FLOW_VEL_v;
           _MESH[index+_opts.MESH_WIDTH].v = _mesh_opts.FLOW_VEL_v;
        }

        //top side
        if( row == _opts.MESH_HEIGHT - 1 )
        {
           _MESH[index].Fluid = 0;
           _MESH[index].u = _mesh_opts.FLOW_VEL_u;
           _MESH[index].v = _mesh_opts.FLOW_VEL_v;
        }

        // left side
        if( col == 0 )
        {
           _MESH[index].Fluid = 0;
           _MESH[index].u = _mesh_opts.FLOW_VEL_u;
           _MESH[index+1].u = _mesh_opts.FLOW_VEL_u;
           _MESH[index].v = _mesh_opts.FLOW_VEL_v;
            
        }

        // right side
        if( col == _opts.MESH_WIDTH - 1 )
        {
           _MESH[index].Fluid = 0;
           _MESH[index].u = _mesh_opts.FLOW_VEL_u;
           _MESH[index].v = _mesh_opts.FLOW_VEL_v;
        }

        if( (row > (_opts.MESH_HEIGHT/2 - 20)) && (row < (_opts.MESH_HEIGHT/2 + 20))  && (col == 0) )
        {

           _MESH[index].density = 0.5;
           //_MESH[index+1].density = 0.5;

        }

        Y = row + _mesh_opts.SLOPE * (col - HORZ_SHIFT) - VERT_SHIFT;

        X = col + _mesh_opts.SLOPE * (col - HORZ_SHIFT) - HORZ_SHIFT;

        if( _mesh_opts.DRAW_ELLIPSE && round((Y) * (Y)*_mesh_opts.VERT_SQUISH + (X) * (X)) <= (R))
        {
           _MESH[index].u = 0;
           _MESH[index].v = 0;
           _MESH[index].Fluid = 0;
        }
    }
    _TEMP = _MESH;

}

//basically just add gravity
// iv'e had trouble making this make any sense, but who needs gravity anyway
void INCO_SOLVER::ExternalForces()
{
    if( _opts.GRAVITY == 0 ){ return; }
    #pragma omp parallel
    {
        #pragma omp for
        for (std::size_t index = _opts.MESH_WIDTH; index < (_opts.MESH_HEIGHT - 1) * (_opts.MESH_WIDTH - 1); ++index)
        {
            _MESH[index].v += _opts.GRAVITY * dt * (_MESH[index].Fluid && _MESH[index-_opts.MESH_WIDTH].Fluid);
        }
        #pragma omp barrier
    }
}

//Steady state solver: 
// just sets divergence of the velocity field to zero for now
void INCO_SOLVER::Divergence()
{
    #pragma omp parallel
    {
        for (std::size_t i = 0; i <= _opts.ITERATIONS; ++i)
        {
            
            ODDCELLS();
            
            EVENCELLS();
            
        }
    }
}

// helper function for Divergence(), makes the divergence zero for the cell at index: index
void INCO_SOLVER::CELLDIVERGENCE(std::size_t index)
{
    // skip non fluid cells
    // also protects MESH bounds
    if ( !_MESH[index].Fluid || index % _opts.MESH_WIDTH == _opts.MESH_WIDTH - 1 )
    {
        return;
    }


    // neighboring cell indecies
    size_t TOP = index + _opts.MESH_WIDTH;
    size_t BOT = index - _opts.MESH_WIDTH;
    size_t LEFT = index - 1;
    size_t RIGHT = index + 1;

    assert( TOP < _MESH.size() );
    assert( BOT > 0 );
    

    uint8_t neighbors = _MESH[TOP].Fluid + _MESH[LEFT].Fluid + _MESH[BOT].Fluid + _MESH[RIGHT].Fluid;
    double Div = - _MESH[index].u - _MESH[index].v + _MESH[TOP].v + _MESH[RIGHT].u;

    // store the divergence values for visualization
    _MESH[index].divergence = Div;

    //vorticity confinement*
    // deffinetly doesn't work right, should be in the momentum term.
    if( _opts.VC != 0 )
        Div += _opts.VC * std::fabs( Curl( _MESH[index].x+_opts.CELL_SIZE/2, _MESH[index].y + _opts.CELL_SIZE/2 ) );

    // Over-Relax divergence values
    Div *= _opts.OVER_RELAXATION;

    // make the divergence of this cell 0 while leaving out boundary edges
    _MESH[index].u += Div * _MESH[LEFT].Fluid / neighbors;
    _MESH[index].v += Div * _MESH[BOT].Fluid / neighbors;
    _MESH[RIGHT].u -= Div * _MESH[RIGHT].Fluid / neighbors;
    _MESH[TOP].v -= Div * _MESH[TOP].Fluid / neighbors;

    //estimate kinetic pressure
    _MESH[index].p -= Div / neighbors * _opts.REF_density * _opts.CELL_SIZE / dt;
}


//These two functions just split up the Divegence() into a stagerd grid
// to avoid race conditions in parallel
    void INCO_SOLVER::ODDCELLS()
    {

        #pragma omp for
        for( std::size_t row = 0; row < _opts.MESH_HEIGHT; ++row )
        {
            for( std::size_t col = (row%2==1) ; col < _opts.MESH_WIDTH; col += 2 )
            {

                CELLDIVERGENCE( row*_opts.MESH_WIDTH + col );

            }
        }

    }

    void INCO_SOLVER::EVENCELLS()
    {

    #pragma omp for
        for( std::size_t row = 0; row < _opts.MESH_HEIGHT; ++row )
        {
            for( std::size_t col = (row%2==0); col < _opts.MESH_WIDTH; col  += 2)
            {

                CELLDIVERGENCE( row*_opts.MESH_WIDTH + col );

            }
        }
    }

//bi-linear cell interpolation to get mesh properties at any point (x,y) within the domain
cell INCO_SOLVER::InterpolateCell( double x, double y )
{
    //this would be much easer if everyting was stored in the cell centroid

    //find the cell that contains this point
    size_t index = CoordToCell( x, y );

    //cant interpolate edge cells, so we just throw them back lol
    // this treats edge cells as constant, theoretically you could 
    // use right difference scheme for left boundary, left differencing
    // for right boundary, etc. this would create weird deadzones that 
    // don't line up, so no thank you.
    if( !_TEMP[index].Fluid )
    {
        return _TEMP[index];
    }

    cell CELL;

    assert( index +_opts.MESH_WIDTH < _MESH.size() );
    assert( index - _opts.MESH_WIDTH >= 0 );
    
    // obtain 0-1 weights from x and y
    size_t U_SHIFT = 0;
    size_t V_SHIFT = 0;

    // correct for the veloceties being stored at cell edge centers:
    //sets the correct bottom left cell index for interpolation grid
    if( (y - _TEMP[index].y ) < h )
    {
        U_SHIFT = -_opts.MESH_WIDTH;
    }
    
    if( (x - _TEMP[index].x ) < h )
    {
        V_SHIFT = -1;
    }


    // horizontal edge velocities
    double x_U = (x - _TEMP[ index +U_SHIFT].x ) / _opts.CELL_SIZE;
    double y_U = (y - h - _TEMP[ index +U_SHIFT].y) / _opts.CELL_SIZE;

    double BL_U = (1 - x_U) * (1 - y_U);
    double BR_U = (x_U) * (1 - y_U);
    double TL_U = (1 - x_U) * y_U;
    double TR_U = x_U * y_U;

    CELL.u = _TEMP[ index +U_SHIFT].u * (BL_U) + _TEMP[ index +U_SHIFT + 1].u * (BR_U) + _TEMP[ index +U_SHIFT + _opts.MESH_WIDTH].u * (TL_U) + _TEMP[ index +U_SHIFT + _opts.MESH_WIDTH + 1].u * (TR_U);

    // vertical edge velocities
    double x_V = (x - h - _TEMP[ index + V_SHIFT].x ) / _opts.CELL_SIZE;
    double y_V = (y - _TEMP[ index + V_SHIFT].y) / _opts.CELL_SIZE;

    double BL_V = (1 - x_V) * (1 - y_V);
    double BR_V = (x_V) * (1 - y_V);
    double TL_V = (1 - x_V) * y_V;
    double TR_V = x_V * y_V;

    CELL.v = _TEMP[ index +V_SHIFT].v * (BL_V) + _TEMP[ index +V_SHIFT + 1].v * (BR_V) + _TEMP[ index +V_SHIFT + _opts.MESH_WIDTH].v * (TL_V) + _TEMP[ index +V_SHIFT + _opts.MESH_WIDTH + 1].v * (TR_V);

    x = (x - h - _TEMP[ index + V_SHIFT + U_SHIFT].x ) / _opts.CELL_SIZE;
    y = (y - h - _TEMP[ index + U_SHIFT + V_SHIFT].y ) / _opts.CELL_SIZE;

    double BL = (1 - x) * (1 - y);
    double BR = (x) * (1 - y);
    double TL = (1 - x) * y;
    double TR = x * y;

    CELL.density = _TEMP[ index + U_SHIFT + V_SHIFT ].density * BL + _TEMP[ index + U_SHIFT + V_SHIFT + 1].density * BR + _TEMP[ index + U_SHIFT + V_SHIFT + _opts.MESH_WIDTH].density * TL + _TEMP[ index + U_SHIFT + V_SHIFT + _opts.MESH_WIDTH + 1].density * TR;
    
    return CELL;
}

//calculates the time step needed to achive the couraunt number 
// specified in the solver options.
void INCO_SOLVER::StableTimeStep()
{
    //maybe unecissary for semi-lagrangian advection, but i like to have the option
    
    double max_speed = 0;
    double temp_speed = 0;
    for( const cell &CELL : _MESH )
    {
        temp_speed = std::sqrt( ( CELL.u*CELL.u + CELL.v*CELL.v ) );
        if( max_speed < temp_speed )
            max_speed = temp_speed;
    }
    dt = ( _opts.CFL * _opts.CELL_SIZE ) / max_speed;
}

//approximate curl field of the velocity
void INCO_SOLVER::CalculateCurl()
{
    for( cell &CELL : _MESH )
    {
        if( CELL.y < (_opts.MESH_HEIGHT-2)* _opts.CELL_SIZE )
            CELL.curl = Curl( CELL.x + h, CELL.y + h );
    }
}

double INCO_SOLVER::Curl( double x, double y )
{
    // curl(F) = dv/dx - du/dy
    return dvdx( x, y, h/2 ) - dudy( x, y, h/2 );
    
}

// central differencing approx
    double INCO_SOLVER::dudy( double x, double y, double delta )
    {

        return ( InterpolateCell(x,y+delta).u - InterpolateCell(x,y-delta).u ) / (2*delta);

    }

    double INCO_SOLVER::dvdx( double x, double y, double delta )
    {

        return ( InterpolateCell(x+delta,y).v - InterpolateCell(x-delta,y).v ) / (2*delta);

    }



// Semi-Lagrangian scheme to approximate time step:
void INCO_SOLVER::Advect()
{
    _TEMP = _MESH;
    #pragma omp parallel
    {
        #pragma omp for
        for (std::size_t index = 0; index < (_opts.MESH_HEIGHT-1) * (_opts.MESH_WIDTH-1); ++index)
        {

            // skip non-fluid cells
            if( !_MESH[index].Fluid )
            {
                continue;
            }

            //different Departure points for different storage locations within each cell
            vec vert, horz, central;

            // update cell centroids
            central = Advect_Step_Back( _TEMP[index].x + h, _TEMP[index].y + h, InterpolateCell( _TEMP[index].x + h, _TEMP[index].y + h ).u * dt, InterpolateCell( _TEMP[index].x + h, _TEMP[index].y + h ).v * dt );

            _MESH[index].density = InterpolateCell( central.x, central.y ).density;

            //why does this cause asymetry?
            // if( !_MESH[index-1].Fluid || !_MESH[index-_opts.MESH_WIDTH].Fluid )
            //     continue;

            //update only if the edge only touches fluid cells
            if( _MESH[index-_opts.MESH_WIDTH].Fluid )
            {
                vert = Advect_Step_Back( _TEMP[index].x + h, _TEMP[index].y, InterpolateCell(_TEMP[index].x + h,_TEMP[index].y).u*dt, _TEMP[index].v*dt );
                _MESH[index].v = InterpolateCell( vert.x, vert.y ).v;
            }

            if( _MESH[index-1].Fluid )
            {
                horz = Advect_Step_Back( _TEMP[index].x, _TEMP[index].y + h, _TEMP[index].u*dt, InterpolateCell(_TEMP[index].x,_TEMP[index].y + h).v*dt );
                _MESH[index].u = InterpolateCell( horz.x, horz.y ).u;

                _MESH[index].time_residual = _MESH[index].x - horz.x - _MESH[index].u * dt;
            }

        }
    }

}

// finds the departure point
vec INCO_SOLVER::Advect_Step_Back( double x, double y, double dx, double dy )
{
    //need local variables foe parallelization 
    vec VAL;
    vec LastVAL;
    cell TEMP_CELL;
    VAL.x = dx;
    VAL.y = dy;

    // iteratively search for the departure point
    for( size_t i = 0; i < 500; ++i )
    {
        LastVAL = VAL;
        TEMP_CELL = InterpolateCell( x - VAL.x, y - VAL.y );
        VAL.x = dt * TEMP_CELL.u;
        VAL.y = dt * TEMP_CELL.v;
        if( std::fabs(LastVAL.x - VAL.x) < 0.00005 && std::fabs(LastVAL.y - VAL.y) < 0.00005 )
            break; 
    }

    VAL.x = x - VAL.x;
    VAL.y = y - VAL.y;

    return VAL;
}

const cell INCO_SOLVER::GetMax()
{
    cell max;
    for (const cell &CELL : _MESH)
    {
        if (CELL.p > max.p)
            max.p = CELL.p;
        if (CELL.u > max.u)
            max.u = CELL.u;
        if (CELL.v > max.v)
            max.v = CELL.v;
        if (CELL.density > max.density)
            max.density = CELL.density;
        if (CELL.divergence > max.divergence)
            max.divergence = CELL.divergence;
        if (CELL.curl > max.curl)
            max.curl = CELL.curl;
        if( CELL.time_residual > max.time_residual )
            max.time_residual = CELL.time_residual;
    }
    return max;
}

const cell INCO_SOLVER::GetMin()
{
    cell min;
    for (const cell &CELL : _MESH)
    {
        if (CELL.p < min.p)
            min.p = CELL.p;
        if (CELL.u < min.u)
            min.u = CELL.u;
        if (CELL.v < min.v)
            min.v = CELL.v;
        if (CELL.density < min.density)
            min.density = CELL.density;
        if (CELL.divergence < min.divergence)
            min.divergence = CELL.divergence;
        if ( CELL.curl < min.curl ) 
            min.curl = CELL.curl;
    }
    return min;
}

const cell &INCO_SOLVER::GetCell(size_t index)
{
    return _MESH[index];
}

// creates a smoothed grid for visualization
std::vector<cell> INCO_SOLVER::InterpolateMesh( size_t Pix_Per_Cell )
{
    size_t index = 0;
    std::vector<cell> NEW_MESH( _MESH.size() * Pix_Per_Cell * Pix_Per_Cell );

    _TEMP = _MESH;

    double step = _opts.CELL_SIZE / Pix_Per_Cell;

    for( double y = step/2; y <= (_opts.MESH_HEIGHT) * _opts.CELL_SIZE; y += step )
    {
        for ( double x = step/2; x <= (_opts.MESH_WIDTH) * _opts.CELL_SIZE; x += step )
        {
            //if( index == NEW_MESH.size() )
                //return NEW_MESH;

            NEW_MESH[ index ] = InterpolateCell( x, y );

            ++index;
        }
    }

    return NEW_MESH;

}

void Graphics::PrintBmp()
{
    std::string filename = "Solution_Data_2/Iteration" + std::to_string(_SOLVER.GetStep()) + ".bmp";
    
    cell MAX = _SOLVER.GetMax();
    cell MIN = _SOLVER.GetMin();

    //std::cout << "v: [" << std::to_string( MIN.v ) << ", " << std::to_string( MAX.v ) << "]" << std::endl;
    
    MAX.u -= MIN.u * (MIN.u < 0);
    MAX.v -= MIN.v * (MIN.v < 0);
    MAX.p -= MIN.p * (MIN.p < 0);
    MAX.density -= MIN.density * (MIN.density < 0);
    MAX.divergence -= MIN.divergence * (MIN.divergence < 0);
    MAX.curl -= MIN.curl * (MIN.curl < 0);

    //assert(MAX.divergence < 1 && "Divergence should not be this high");

    double BIGGEST_VEL = sqrt(MAX.u * MAX.u + MAX.v * MAX.v);
    double BIGGEST = (MAX.p * _opts.PRINT_P + MAX.u * _opts.PRINT_U + MAX.v * _opts.PRINT_V + MAX.divergence * _opts.PRINT_DIV + MAX.curl * _opts.PRINT_CURL + BIGGEST_VEL * _opts.PRINT_VEL + _opts.PRINT_DEN );
    double VEL = 0; 

    //assert( BIGGEST > 0 );

    double PRINT_VAL;

    BmpHeader header;
    IBmpInfoHeader infoHeader;

    infoHeader.Height = _SOLVER.GetMeshHeight() * _opts.PixelsPerCell;
    infoHeader.Width = _SOLVER.GetMeshWidth() * _opts.PixelsPerCell;
    infoHeader.horizontalResolution = 2 * abs(infoHeader.Width); // fits the image to 50 cm
    infoHeader.VerticalResolution = 2 * abs(infoHeader.Height);

    std::string Padding( infoHeader.Width % 4 , '0');

    header.bfSize = 54 + infoHeader.Width * infoHeader.Height * 3;

    std::ofstream file(filename, std::ios::binary);

    header.HeaderWrite(file);
    infoHeader.InfoHeaderWrite(file);

    Pixel pixel;
    
    auto InterpolatedMESH = _SOLVER.InterpolateMesh( _opts.PixelsPerCell );

    for (size_t index = 0; index < InterpolatedMESH.size(); ++index)
    {

        const cell &VALUE = InterpolatedMESH[index];//_SOLVER.GetCell(index);
        // Assign color based on tile value

        VEL = sqrt(VALUE.u * VALUE.u + VALUE.v * VALUE.v);



        PRINT_VAL = ((VALUE.p - MIN.p * (MIN.p < 0)) * _opts.PRINT_P + 
                    (VALUE.u - MIN.u * (MIN.u < 0)) * _opts.PRINT_U + 
                    (VALUE.v - MIN.v * (MIN.v < 0)) * _opts.PRINT_V + 
                    (VALUE.divergence - MIN.divergence * (MIN.divergence < 0)) * _opts.PRINT_DIV + 
                    (VALUE.curl - MIN.curl * (MIN.curl < 0)) * _opts.PRINT_CURL +
                    VEL * _opts.PRINT_VEL) / BIGGEST;


        // if ( PRINT_VAL > 1 || (PRINT_VAL) < 0)
        // {
        //     std::cout << "PRINT_VAL: " + std::to_string(PRINT_VAL) << std::endl;
        // }

        if (!_opts.PRINT_DEN)
        {
            pixel.Red = (uint8_t)((PRINT_VAL * 255));
            pixel.Green = (uint8_t)((PRINT_VAL * 0));
            pixel.Blue = (uint8_t)(255 - (PRINT_VAL * 255));
        }
        else
        {
            pixel.Red = (uint8_t)((VALUE.density * 255));
            pixel.Green = (uint8_t)((VALUE.density * 255));
            pixel.Blue = (uint8_t)((VALUE.density * 255)); 
        }

        if (!VALUE.Fluid)
        {
            pixel.Red = 0;
            pixel.Blue = 0;
            pixel.Green = 0;
        }

        pixel.PixelWrite(file);
        
        // extra "padding" cuz rows must be a multiple of 4 (WHYYY?)
        if(index % infoHeader.Width == (size_t)infoHeader.Width-1 )
            file.write((char *)&Padding, infoHeader.Width % 4); // pro tip width != Height

    }

    std::cout << "Saving: " << filename << ";  divergence: " << std::to_string(MAX.divergence) << "; time Residual: " << std::to_string( MAX.time_residual ) << std::endl;

    double min_u = 0;
    double max_u = 0;

    for( cell &CELL : InterpolatedMESH )
    {
        
        if( CELL.u < min_u )
            min_u = CELL.u;
        if( CELL.u > max_u )
            max_u = CELL.u;

    }


    std::cout << "maximum u: " << std::to_string( max_u ) << "; minimum u: " << std::to_string( min_u ) << std::endl;

    std::cout << "maximum u: " << std::to_string( MAX.u ) << "; minimum u: " << std::to_string( MIN.u ) << std::endl;

    file.close();
}

void IBmpInfoHeader::InfoHeaderWrite(std::ofstream &file)
{
    file.write((char *)&this->Size, sizeof(uint32_t));
    file.write((char *)&this->Width, sizeof(int32_t));
    file.write((char *)&this->Height, sizeof(int32_t));
    file.write((char *)&this->Planes, sizeof(uint16_t));
    file.write((char *)&this->colorDepth, sizeof(uint16_t));
    file.write((char *)&this->Compression, sizeof(uint32_t));
    file.write((char *)&this->SizeImage, sizeof(uint32_t));
    file.write((char *)&this->horizontalResolution, sizeof(int32_t));
    file.write((char *)&this->VerticalResolution, sizeof(int32_t));
    file.write((char *)&this->ClrTable, sizeof(uint32_t));
    file.write((char *)&this->ClrImportant, sizeof(uint32_t));
}

void Pixel::PixelWrite(std::ofstream &file)
{
    file.write((char *)&this->Blue, 1);
    file.write((char *)&this->Green, 1);
    file.write((char *)&this->Red, 1);
}

double Force::REPORT(const std::vector<std::vector<cell>> &MESH, double CELL_SIZE)
{
    double Report_val = 0;
    for (size_t r = 1; r < MESH.size() - 1; ++r)
    {
        for (size_t c = 1; c < MESH[0].size() - 1; ++c)
        {
            Report_val += MESH[r][c].p * CELL_SIZE * (MESH[r][c + 1].Tag == this->Tag);
            Report_val -= MESH[r][c].p * CELL_SIZE * (MESH[r][c - 1].Tag == this->Tag);
        }
    }
    if (_PrintToTerm)
    {

        std::cout << Name << ": " << std::to_string(Report_val) << std::endl;
    }
    if (_PrintToFile)
    {

        // Make a report file output structure
    }
    return Report_val;
}

#endif