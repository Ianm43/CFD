#ifndef INCO_H
#define INCO_H

#include <stdint.h>
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

struct cell
{
    double x = 0; // x coordinate of the bottom left vertex
    double y = 0;
    double u = 0;       // x velocity in m/s
    double v = 0;       // y velocity in m/s
    double density = 0; // SMOKE DENSITY NOT PHYSICAL DENSITY
    double divergence = 0;
    double p = 0;   // estimated kinetic pressure for visulaization
    bool Fluid = 1; // pretty self explanatory
    std::string Tag = "Fluid";
};

struct INCO_opts
{
    size_t TIMESTEPS = 10;
    size_t ITERATIONS = 40;
    size_t REPORT_INTERVAL = 100;
    double CELL_SIZE = 0.01;
    size_t MESH_WIDTH = 100;
    size_t MESH_HEIGHT = 100;
    double OVER_RELAXATION = 1.9;
    double GRAVITY = 0;
    double REF_density = 1.225;
    double TIMESTEP = 0.01;
};

struct Force
{
    std::string Tag;
    double x_norm;
    double y_norm;
    std::string Name;
    bool _PrintToFile = 0;
    bool _PrintToTerm = 0;

    Force(double x, double y, std::string CELL_TAG, std::string name, bool term, bool file) : Tag(CELL_TAG), x_norm(x), y_norm(y), Name(name), _PrintToFile(file), _PrintToTerm(term) {};

    double REPORT(const std::vector<std::vector<cell>> &MESH, double CELL_SIZE);
};

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

struct Graphics_opts
{
    bool PRINT_P = 1;
    bool PRINT_V = 0;
    bool PRINT_U = 0;
    bool PRINT_DEN = 0;
    bool PRINT_DIV = 0;
    bool PRINT_VEL = 0;
};

struct MESH_opts
{

    bool DRAW_ELLIPSE = true;
    double SLOPE = 0;
    double VERT_SQUISH = 1;
    double FLOW_VEL_u = 1;
    double FLOW_VEL_v = 0;

};

class INCO_SOLVER
{
    private:
        INCO_opts _opts;
        MESH_opts _mesh_opts;
        std::vector<cell> _MESH;
        std::vector<cell> _TEMP;
        std::size_t _Step; // curent time step

        void ExternalForces();
        void CELLDIVERGENCE(std::size_t index);
        void ODDCELLS();
        void EVENCELLS();
        void Divergence();
        cell InterpolateCell(std::size_t index, double x, double y);
        void Advect();

        std::size_t CoordToCell( double x, double y )
        {
            // convert from x,y coordinates to the cell that contains those coordinates
            return floor(y / _opts.CELL_SIZE) * _opts.MESH_WIDTH + (x / _opts.CELL_SIZE);
        }

    public:
        // I get the feelling that an object shouldn't be passing a reference to itself to one of it's members but oh well
        INCO_SOLVER(INCO_opts Options, MESH_opts M_Options ) : 
        _opts(Options), _mesh_opts( M_Options ), _MESH( _opts.MESH_HEIGHT*_opts.MESH_WIDTH ) , _Step(0) {};
        void Solve();
        void Make_MESH();
        void Initilaize();
        const cell GetMax();
        const cell GetMin();
        const cell &GetCell(size_t index);

        const size_t &GetMeshHeight() { return _opts.MESH_HEIGHT; }
        const size_t &GetMeshWidth() { return _opts.MESH_WIDTH; }
        const size_t &GetStep() { return _Step; }
};

class Graphics
{
    Graphics_opts _opts;
    INCO_SOLVER &_SOLVER;

public:
    Graphics( INCO_SOLVER &solver, Graphics_opts &opts ): _opts(opts), _SOLVER( solver ){};
    void PrintBmp();
};

void INCO_SOLVER::Initilaize()
{
    Divergence();
    for( size_t i = 0; i < 10000 && std::abs(GetMax().divergence) > 0.02; i += _opts.ITERATIONS )
        Divergence();

}

void INCO_SOLVER::Make_MESH()
{
    size_t row;
    size_t col;

    double X, Y;

    // super dumb way of making a skewed ellipse 
    //for help defining varibles: https://www.desmos.com/calculator/fgdc56e9nm

    double R = (float)_opts.MESH_HEIGHT * _opts.MESH_HEIGHT / 100; // semi major axis
    double VERT_SHIFT = _opts.MESH_HEIGHT / 2;
    double HORZ_SHIFT = _opts.MESH_WIDTH / 5;

    for( size_t index = 0; index <_MESH.size(); ++index )
    {
        col = index % _opts.MESH_WIDTH;
        row = index / _opts.MESH_WIDTH;

       _MESH[index].x = index % _opts.MESH_WIDTH * _opts.CELL_SIZE;
       _MESH[index].y = index / _opts.MESH_WIDTH * _opts.CELL_SIZE;

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

}

void INCO_SOLVER::ExternalForces()
{
    //#pragma omp parallel
    {
        //#pragma omp for
        for (std::size_t index = 0; index < (_opts.MESH_HEIGHT - 1) * (_opts.MESH_WIDTH - 1); ++index)
        {
            _MESH[index].v += _opts.GRAVITY * _opts.TIMESTEP * _MESH[index].Fluid;
        }
        //#pragma omp barrier
    }
}

void INCO_SOLVER::CELLDIVERGENCE(std::size_t index)
{
    // skip non fluid cells
    // also protects MESH bounds
    if (!_MESH[index].Fluid)
    {
        return;
    }

    // neighboring cell indecies
    size_t TOP = index + _opts.MESH_WIDTH;
    size_t BOT = index - _opts.MESH_WIDTH;
    size_t LEFT = index - 1;
    size_t RIGHT = index + 1;

    uint8_t neighbors = _MESH[TOP].Fluid + _MESH[LEFT].Fluid + _MESH[BOT].Fluid + _MESH[RIGHT].Fluid;
    double Div = - _MESH[index].u - _MESH[index].v + _MESH[TOP].v + _MESH[RIGHT].u;

    // store the divergence values for visualization
    _MESH[index].divergence = Div;

    // Over-Relax divergence values
    Div *= _opts.OVER_RELAXATION;

    // make the divergence of this cell 0 while leaving out boundary edges
    _MESH[index].u += Div * _MESH[LEFT].Fluid / neighbors;
    _MESH[index].v += Div * _MESH[BOT].Fluid / neighbors;
    _MESH[RIGHT].u -= Div * _MESH[RIGHT].Fluid / neighbors;
    _MESH[TOP].v -= Div * _MESH[TOP].Fluid / neighbors;

    // estimate kinetic pressure
    _MESH[index].p -= Div / neighbors * _opts.REF_density * _opts.CELL_SIZE / _opts.TIMESTEP;
}

void INCO_SOLVER::Divergence()
{
    #pragma omp parallel
    {
        for (std::size_t i = 0; i <= _opts.ITERATIONS; ++i)
        {
            
            /* for( size_t index = 0; index < _MESH.size(); index++ )
            {

                CELLDIVERGENCE( index );

            } */
            
            ODDCELLS();
            
            EVENCELLS();
            
        }
    }
}

void INCO_SOLVER::ODDCELLS()
{
    #pragma omp for
    for (std::size_t index = 1; index < (_opts.MESH_HEIGHT) * (_opts.MESH_WIDTH); index += 2)
    {
        CELLDIVERGENCE(index);
        //_MESH[index].v = 1;
    }
}

void INCO_SOLVER::EVENCELLS()
{
    #pragma omp for
    for (std::size_t index = 0; index < (_opts.MESH_HEIGHT) * (_opts.MESH_WIDTH); index += 2)
    {
        CELLDIVERGENCE(index);
        //_MESH[index].v = -1;
    }
}

cell INCO_SOLVER::InterpolateCell(std::size_t index, double x, double y)
{
    // bi-linear cell interpolation
    cell CELL;

    if (index + _opts.MESH_WIDTH >= _MESH.size())
    {
        return CELL;
    }

    double BL = (1 - x) * (1 - y);
    double BR = (x) * (1 - y);
    double TL = (1 - x) * (y);
    double TR = x * y;

    CELL.u = _MESH[index].u * BL + _MESH[index + 1].u * BR + _MESH[index + _opts.MESH_WIDTH].u * TL + _MESH[index + _opts.MESH_WIDTH + 1].u * TR;
    CELL.v = _MESH[index].v * BL + _MESH[index + 1].v * BR + _MESH[index + _opts.MESH_WIDTH].v * TL + _MESH[index + _opts.MESH_WIDTH + 1].v * TR;
    CELL.density = _MESH[index].density * BL + _MESH[index + 1].density * BR + _MESH[index + _opts.MESH_WIDTH].density * TL + _MESH[index + _opts.MESH_WIDTH + 1].density * TR;

    return CELL;
}

void INCO_SOLVER::Advect()
{
    _TEMP = _MESH;
//#pragma omp parallel
{
    //#pragma omp for
        for (std::size_t index = 0; index < (_opts.MESH_HEIGHT) * (_opts.MESH_WIDTH); ++index)
        {

            // ignore non-fluid cells
            if (!_MESH[index].Fluid)
            {
                continue;
            }

            double X_vert, Y_vert, X_horz, Y_horz;
            std::size_t density_index, u_index, v_index;

            // get average velocities around each edge and look back

            X_vert = _MESH[index].x - (_MESH[index].u + _MESH[index + 1].u + _MESH[index + _opts.MESH_WIDTH].u + _MESH[index + _opts.MESH_WIDTH + 1].u) / 4.0 * _opts.TIMESTEP;

            Y_vert = _MESH[index].y - _MESH[index].v * _opts.TIMESTEP;

            X_horz = _MESH[index].x - _MESH[index].u * _opts.TIMESTEP;

            Y_horz = _MESH[index].y - (_MESH[index].v + _MESH[index - 1].v + _MESH[index + _opts.MESH_WIDTH].v + _MESH[index + _opts.MESH_WIDTH - 1].v) / 4.0 * _opts.TIMESTEP;

            // round negative values to 0
            X_vert *= (X_vert > 0);
            Y_vert *= (Y_vert > 0);
            X_horz *= (X_horz > 0);
            Y_horz *= (Y_horz > 0);

            // round values greater than meshh dimensions to the edge
            if (X_vert >= (_opts.MESH_WIDTH) * _opts.CELL_SIZE)
            {
                X_vert = (_opts.MESH_WIDTH-1) * _opts.CELL_SIZE;
            }
            if (X_horz >= (_opts.MESH_WIDTH) * _opts.CELL_SIZE)
            {
                X_horz = (_opts.MESH_WIDTH-1) * _opts.CELL_SIZE;
            }
            if (Y_horz >= (_opts.MESH_HEIGHT) * _opts.CELL_SIZE)
            {
                Y_horz = (_opts.MESH_HEIGHT-1) * _opts.CELL_SIZE;
            }
            if (Y_vert >= (_opts.MESH_HEIGHT) * _opts.CELL_SIZE)
            {
                Y_vert = (_opts.MESH_HEIGHT-1) * _opts.CELL_SIZE;
            }


            assert( X_vert >= 0 && X_vert < _opts.MESH_WIDTH * _opts.CELL_SIZE );

            // temporarily store cell indecies
            u_index = CoordToCell(X_horz, Y_horz);
            v_index = CoordToCell(X_vert, Y_vert);
            density_index = CoordToCell(X_horz, Y_vert);


            if( !(v_index >= 0 && v_index < _opts.MESH_HEIGHT*_opts.MESH_WIDTH) )
            {
                std::cout << std::to_string( X_vert ) << ", " << std::to_string( Y_vert ) << std::endl;
            }

            assert( u_index >= 0 && u_index < _opts.MESH_HEIGHT*_opts.MESH_WIDTH );

            assert( density_index >= 0 && density_index < _opts.MESH_HEIGHT*_opts.MESH_WIDTH );

            // get average velocities around position
            _TEMP[index].u = InterpolateCell(u_index, (X_horz - _MESH[u_index].x) / _opts.CELL_SIZE, (Y_horz - _MESH[u_index].y) / _opts.CELL_SIZE).u;
            _TEMP[index].v = InterpolateCell(v_index, (X_vert - _MESH[v_index].x) / _opts.CELL_SIZE, (Y_vert - _MESH[v_index].y) / _opts.CELL_SIZE).v;
            _TEMP[index].density = InterpolateCell(density_index, (X_horz - _MESH[density_index].x) / _opts.CELL_SIZE, (Y_vert - _MESH[density_index].y) / _opts.CELL_SIZE).density;
        }
    //#pragma omp barrier
}
    _MESH = _TEMP;
}

void INCO_SOLVER::Solve()
{
    _Step++;

    ExternalForces();
        
    Divergence();
        
    Advect();       
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
    }
    return min;
}

const cell &INCO_SOLVER::GetCell(size_t index)
{
    return _MESH[index];
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

    //assert(MAX.divergence < 1 && "Divergence should not be this high");

    double BIGGEST_VEL = sqrt(MAX.u * MAX.u + MAX.v * MAX.v);
    double BIGGEST = (MAX.p * _opts.PRINT_P + MAX.u * _opts.PRINT_U + MAX.v * _opts.PRINT_V + MAX.divergence * _opts.PRINT_DIV + BIGGEST_VEL * _opts.PRINT_VEL + _opts.PRINT_DEN);
    double VEL = 0; 

    //assert( BIGGEST > 0 );

    double PRINT_VAL;

    BmpHeader header;
    IBmpInfoHeader infoHeader;

    infoHeader.Height = _SOLVER.GetMeshHeight();
    infoHeader.Width = _SOLVER.GetMeshWidth();
    infoHeader.horizontalResolution = 2 * abs(infoHeader.Width); // fits the image to 50 cm
    infoHeader.VerticalResolution = 2 * abs(infoHeader.Height);

    std::string Padding( infoHeader.Width % 4 , '0');

    header.bfSize = 54 + infoHeader.Width * infoHeader.Height * 3;

    std::ofstream file(filename, std::ios::binary);

    header.HeaderWrite(file);
    infoHeader.InfoHeaderWrite(file);

    Pixel pixel;

    for (size_t index = 0; index < (size_t)(infoHeader.Height * infoHeader.Width); ++index)
    {
        
        const cell &VALUE = _SOLVER.GetCell(index);
        // Assign color based on tile value

        VEL = sqrt(VALUE.u * VALUE.u + VALUE.v * VALUE.v);



        PRINT_VAL = ((VALUE.p - MIN.p * (MIN.p < 0)) * _opts.PRINT_P + 
                    (VALUE.u - MIN.u * (MIN.u < 0)) * _opts.PRINT_U + 
                    (VALUE.v - MIN.v * (MIN.v < 0)) * _opts.PRINT_V + 
                    (VALUE.divergence - MIN.divergence * (MIN.divergence < 0)) * _opts.PRINT_DIV + 
                    VEL * _opts.PRINT_VEL) / BIGGEST;


        if ( PRINT_VAL > 1 || (PRINT_VAL) < 0)
        {
            std::cout << "PRINT_VAL: " + std::to_string(PRINT_VAL) << std::endl;
        }

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

    std::cout << "Saving: " << filename << ";  divergence: " << std::to_string(MAX.divergence) << std::endl;

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