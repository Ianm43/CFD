#include "inco.h"

static const size_t MESH_width = 800;
static const size_t MESH_height = 400;
static const double MESH_SIZE = 0.005; // size of each cell in m
static const double FLOW_VEL = 2;
static const double density = 1.225; // density in kg / m^3

static const double dt = 0.05/3 * 0.2 ; // time step in sec should be < MESH_SIZE / Flow Velocity

static const bool PRINT_P = 0;
static const bool PRINT_V = 1;
static const bool PRINT_U = 0;
static const bool PRINT_DEN = 0;
static const bool PRINT_DIV = 0;
static const bool PRINT_VEL = 0;

static std::vector< std::vector< cell > > MESH( MESH_height, std::vector<cell>(MESH_width) );;
static auto TEMP = MESH; 

using namespace std;

int main()
{
    for( size_t r = 0; r < MESH_height; ++r )
    {
        if( r > MESH_height/2 - 30 && r < MESH_height/2 + 30 )
        {
            MESH[r][0].density=0.5;
        }

        MESH[r][0].Fluid = 0;
        MESH[r][MESH_width-1].Fluid = 0;
        
        for( size_t c = 0; c < MESH_width; ++c )
        {
            if( ( ((double)r-MESH_height/2+1)*((double)r-MESH_height/2+1) + ((double)c-MESH_width/4)*((double)c-MESH_width/4) ) < ( MESH_height*MESH_height / 80) )
            {
                MESH[r][c].Fluid = 0;
                continue;
            }
            MESH[r][c].u = FLOW_VEL;

            

            MESH[0][c].Fluid = 0; // bottom treatment
            MESH[0][c].u = 0;
            MESH[MESH_height-1][c].Fluid = 0; // top treatment
            MESH[MESH_height-1][c].u = 0;
        }

    }

    Divergence( 1000, 1.9 );

    size_t t = 0;

    while( t < 10000 )
    {
        //ExternalForce();

        Divergence( 100, 1.9 );

        Advection();

        if( t % 100 == 0 )
        {
            std::string name = "frame_" + to_string(t) +".bmp";

            Bitmap( MESH, name );
        }

        ++t;
    }

    return 0;
}

void ExternalForce()
{

    for( auto &row : MESH )
    {
        for( auto &CELL : row )
        {
            CELL.v -= 9.8 * dt * ( CELL.Fluid );
        }
    }

}

void Divergence( size_t iterations, double Relaxation )
{
    uint8_t neighbors = 0;
    double Div = 0;

    for( size_t i = 0; i <= iterations; ++i )
    {
        for( size_t row = 1; row < MESH_height - 1; ++row )
        {
            for( size_t col = 1; col < MESH_width - 1; ++col  )
            {
                if( !MESH[row][col].Fluid ){ continue; }
                neighbors = MESH[row-1][col].Fluid + MESH[row+1][col].Fluid + MESH[row][col+1].Fluid + MESH[row][col-1].Fluid;
                Div = (-MESH[row][col].u - MESH[row][col].v + MESH[row+1][col].v + MESH[row][col+1].u);
                
                MESH[row][col].divergence = Div;

                Div *= Relaxation;

                MESH[row][col].u += Div * (double)MESH[row-1][col].Fluid / (double)neighbors;
                MESH[row][col].v += Div * (double)MESH[row][col-1].Fluid / (double)neighbors;
                MESH[row+1][col].v -= Div * (double)MESH[row+1][col].Fluid / (double)neighbors;
                MESH[row][col+1].u -= Div * (double)MESH[row][col+1].Fluid / (double)neighbors;

                MESH[row][col].p += Div/(double)neighbors * (density * MESH_SIZE )/dt;
                
            }
        }

    }

}

void Advection()
{
    // Update the edge velocities for each cell
    TEMP = MESH;

    double X_vert, Y_vert, X_horz, Y_horz;
    long long r_vert, c_vert, r_horz, c_horz, r_den, c_den;

    for( size_t r = 1; r < MESH_height - 1; ++r )
    {
        for( size_t c = 1; c < MESH_width - 1; ++c )
        {
            if( !MESH[r][c].Fluid ){ continue; }
            //get average velocities around each edge and look back

                X_vert = (double)c*MESH_SIZE - ( (MESH[r][c].u + MESH[r-1][c].u + MESH[r][c+1].u + MESH[r-1][c+1].u) / 4.0 ) * dt;
                X_vert *= !( X_vert < 0 );
                Y_vert = (double)r*MESH_SIZE - MESH[r][c].v * dt;
                Y_vert *= !( Y_vert < 0 );

                X_horz = (double)c*MESH_SIZE - MESH[r][c].u * dt;
                X_horz *= !( X_horz < 0 );
                Y_horz = (double)r*MESH_SIZE - ( (MESH[r][c].v + MESH[r][c-1].v + MESH[r+1][c].v + MESH[r+1][c-1].v) / 4.0 ) * dt;
                Y_horz *= !( Y_horz < 0 );

            // predict last location of "particle"

                c_vert = ( X_vert / MESH_SIZE );
                //c_vert *= !( c_vert < 0 ); 
                if( c_vert > (long long)MESH_width - 2 ){ c_vert = MESH_width - 2; }
                assert( c_vert <= (long long)MESH_width - 2 && c_vert >= 0);

                r_vert =  ( Y_vert / MESH_SIZE );
                //r_vert *= !( r_vert < 0 );
                if( r_vert > (long long)MESH_height - 2 ){ r_vert = MESH_height - 2; }

                c_horz =  ( X_horz / MESH_SIZE );
                //c_horz *= !( c_horz < 0 );
                if( c_horz > (long long)MESH_width - 2 ){ c_horz = MESH_width - 2; }

                r_horz = ( Y_horz / MESH_SIZE );
                //r_horz *= !( r_horz < 0 );
                if( r_horz > (long long)MESH_height - 2 ){ r_horz = MESH_height - 2; }

                c_den = ( X_horz / MESH_SIZE );
                //c_den *= !( c_den < 0 );
                if( c_den > (long long)MESH_width - 2 ){ c_den = MESH_width - 2; }

                r_den = ( Y_vert / MESH_SIZE );
                //r_den *= !( r_den < 0 );
                if( r_den > (long long)MESH_height - 2 ){ r_den = MESH_height - 2; }

            // get average velocities around position
            /*
                if (X_horz - (double)c_horz * MESH_SIZE < 0)
                {
                    cout << to_string(X_horz - (double)c_horz * MESH_SIZE) << endl;
                    cout << "r: " << to_string(r_horz) << "c: " << to_string(c_horz) << endl;
                    cout << "X_horz: " << to_string( X_horz ) << endl;
                }
            assert( (X_horz - (double)c_horz * MESH_SIZE) >= 0 );

            assert( ( X_horz - (double)c_horz * MESH_SIZE ) / MESH_SIZE >= 0 && ( X_horz - (double)c_horz * MESH_SIZE ) / MESH_SIZE <= 1 );
            */

            TEMP[r][c].u = GetWeightedCellAverage( (size_t)r_horz, (size_t)c_horz, ( X_horz - (double)c_horz * MESH_SIZE ) / MESH_SIZE, ( Y_horz - (double)r_horz*MESH_SIZE ) / MESH_SIZE ).u;//((c-c_horz)*(double)MESH_SIZE - X_horz) / (double)MESH_SIZE, ((r-r_horz)*(double)MESH_SIZE - Y_horz) / (double)MESH_SIZE ).u;

            TEMP[r][c].v = GetWeightedCellAverage( (size_t)r_vert, (size_t)c_vert, ( X_vert - (double)c_vert * MESH_SIZE ) / MESH_SIZE, ( Y_vert - (double)r_vert*MESH_SIZE ) / MESH_SIZE ).v;//((c-c_vert)*(double)MESH_SIZE - X_vert) / (double)MESH_SIZE, ((r-r_vert)*(double)MESH_SIZE - Y_vert) / (double)MESH_SIZE  ).v;

            TEMP[r][c].density = GetWeightedCellAverage( (size_t)r_den, (size_t)c_den, ( X_horz - (double)c_den * MESH_SIZE ) / MESH_SIZE, ( Y_vert - (double)r_den*MESH_SIZE ) / MESH_SIZE).density;//((c-c_den)*(double)MESH_SIZE - X_horz) / (double)MESH_SIZE, ((r-r_den)*(double)MESH_SIZE - Y_vert) / (double)MESH_SIZE ).density;

        }
    }
    MESH = TEMP;
}

cell GetWeightedCellAverage( size_t r, size_t c, double x, double y )// x and y give the local position within the "cell"
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
    BmpHeader header;
    IBmpInfoHeader infoHeader;

    infoHeader.Height = Img_Data.size();
    infoHeader.Width = Img_Data[0].size(); 
    infoHeader.horizontalResolution = 2 * abs(infoHeader.Width); // fits the image to 50 cm
    infoHeader.VerticalResolution = 2 * abs(infoHeader.Height);

    header.bfSize = 54 + abs(infoHeader.Width * infoHeader.Height * 3);

    // extra "padding" cuz rows must be a multiple of 4 (WHYYY?)

    ofstream file( filename, ios::binary);

    header.HeaderWrite( file );
    infoHeader.InfoHeaderWrite( file );

    Pixel pixel;

    cell MAX = MAX_CELL();
    cell MIN = MIN_CELL();
    
    //assert( MAX.p != MIN.p );
    
    MAX.u -= MIN.u * (MIN.u < 0);
    MAX.v -= MIN.v * (MIN.v < 0);
    MAX.p -= MIN.p * (MIN.p < 0);
    MAX.density -= MIN.density * (MIN.density < 0) ;
    MAX.divergence -= MIN.divergence * (MIN.divergence < 0);
    

    double BIGGEST_VEL = sqrt(MAX.u*MAX.u + MAX.v*MAX.v);
    double BIGGEST = ( MAX.p * PRINT_P + MAX.u * PRINT_U + MAX.v * PRINT_V + MAX.divergence * PRINT_DIV + BIGGEST_VEL * PRINT_VEL + PRINT_DEN );
    double VEL = 0;

    double PRINT_VAL;
    
    for( auto const &ROW : Img_Data )
    {
        for( cell const &VALUE : ROW )
        {

            // Assign color based on tile value
            VEL = sqrt( VALUE.u*VALUE.u + VALUE.v*VALUE.v );

            PRINT_VAL = ( (VALUE.p-MIN.p*(MIN.p<0)) * PRINT_P + (VALUE.u-MIN.u*(MIN.p<0)) * PRINT_U + (VALUE.v-MIN.v*(MIN.v<0)) * PRINT_V + (VALUE.divergence-MIN.divergence*(MIN.divergence<0)) * PRINT_DIV + VEL * PRINT_VEL )/BIGGEST + PRINT_DEN * VALUE.density;

            if ((PRINT_VAL * 255) > 255 || (PRINT_VAL * 255) < 0)
            {
                cout << "PRINT_VAL: " + to_string(PRINT_VAL) << endl;
            }

            pixel.Red = (uint8_t)( ( PRINT_VAL * 255 ) );
            pixel.Green = (uint8_t)( ( PRINT_VAL ) );
            pixel.Blue = (uint8_t)( ( 255 - PRINT_VAL * 255 ) );

            if( !VALUE.Fluid )
            {
                pixel.Red = 0; pixel.Blue = 0; pixel.Green = 0;
            }

            pixel.PixelWrite( file );
        }
        for( uint8_t pad = 0; pad < abs(infoHeader.Height) % 4; ++pad )
            file.write( (char *)0, 1 ); // definitly works??!!!
    }

    cout << "MAXIMUM divergence is: " << to_string( MAX.divergence ) << endl;
    
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
