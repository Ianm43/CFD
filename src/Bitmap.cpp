#include "Bitmap.h"

using namespace std;


const int PIXELS_PER_TILE = 1; //must be a multiple of FUCKING FOUR
// add padding to each row = Height % 4

void Bitmap( vector< std::vector< uint8_t > > & Img_Data, std::string &filename )
{
    BmpHeader header;
    IBmpInfoHeader infoHeader;

    infoHeader.Width = Img_Data.size() * PIXELS_PER_TILE;
    infoHeader.horizontalResolution = 2 * infoHeader.Width; // fits the image to 50 cm
    infoHeader.Height = Img_Data[0].size() * PIXELS_PER_TILE; 
    infoHeader.VerticalResolution = 2 * infoHeader.Height;

    header.bfSize = 54 + (infoHeader.Width * infoHeader.Height * 3);

    // extra "padding" cuz rows must be a multiple of 4 (WHYYY?)

    ofstream file( filename, ios::binary);

    header.HeaderWrite( file );
    infoHeader.InfoHeaderWrite( file );

    Pixel pixel;
    
    for( std::vector< uint8_t > &ROW : Img_Data )
    {
        for( uint8_t &VALUE : ROW )
        {
            // Assign color based on tile value
            pixel.Red = (VALUE * 50) % 256;
            pixel.Green = (VALUE * 80) % 256;
            pixel.Blue = (VALUE * 110) % 256;

            pixel.PixelWrite( file );
        }
        for( uint8_t pad = 0; pad < infoHeader.Height % 4; ++pad )
            file.write( (char *)0, 1 ); // definitly works??!!!
    }
    
    file.close();    
}



