#ifndef INCOMPRESSIBLE_INVISCID
#define INCOMPRESSIBLE_INVISCID

#include <stdint.h>
#include <cmath>

#include <cstdint>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>

struct cell
{  
    double u = 0; // x velocity in m/s
    double v = 0; // y velocity in m/s
    double density = 0;
    double divergence = 0;
    double p = 0;
    bool Fluid = 1;
};

int main(); // function declerations should be put before static variables
void ExternalForce();
void Divergence( size_t iterations, double Relaxation );
void Advection();
cell MAX_CELL();
cell MIN_CELL();

cell GetWeightedCellAverage( size_t r, size_t c, double x, double y );
void Bitmap( std::vector< std::vector< cell > > &Img_Data, std::string &filename );





struct BmpHeader//14 bytes
{
    char     bfType[2] = { 'B','M' };
    uint32_t bfSize;
    uint32_t ReservedBytes = 0;
    uint32_t Headeroffset = 54;

    
    void HeaderWrite( std::ofstream &file )
    {
        file.write((char *)&this->bfType, 2);
        file.write((char *)&this->bfSize, sizeof(uint32_t));
        file.write((char *)&this->ReservedBytes, sizeof(uint32_t));
        file.write((char *)&this->Headeroffset, sizeof(uint32_t));
    }
};

struct IBmpInfoHeader//40 bytes
{
    uint32_t Size = 40;
    int32_t  Width; //pixels
    int32_t  Height; //pixels
    uint16_t Planes = 1;
    uint16_t colorDepth = 24;
    uint32_t Compression = 0;
    uint32_t SizeImage = 0;
    int32_t  horizontalResolution;//pixels per meter
    int32_t  VerticalResolution;//pixels per meter
    uint32_t ClrTable = 0;
    uint32_t ClrImportant = 0;

    void InfoHeaderWrite( std::ofstream &file )
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
};

struct Pixel
{
    uint8_t Red;
    uint8_t Blue;
    uint8_t Green;

    void PixelWrite( std::ofstream &file )
    {
        file.write((char *)&this->Blue, 1);
        file.write((char *)&this->Green, 1);
        file.write((char *)&this->Red, 1);
    }
};


#endif