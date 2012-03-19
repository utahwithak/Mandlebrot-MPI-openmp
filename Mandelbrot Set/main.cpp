#include <iostream>
#include "mpi.h"
#include "omp.h"


#define CON_FCTR 255/30
using namespace std;

static double ImageWidth = 1020;
static double ImageHeight = 1020;

int main(int argc, char** argv) {
    double values[(int) ImageHeight][(int)ImageWidth];
	double MinRe		=	-2.0;
	double MaxRe		=	1.0;
	double MinIm		=	-1.2;
	double MaxIm		=	MinIm + (MaxRe - MinRe) * ImageHeight/ImageWidth;
	double Re_factor	=	(MaxRe - MinRe) / (ImageWidth - 1);
	double Im_factor	=	(MaxIm - MinIm) / (ImageHeight - 1);
	unsigned MaxIterations = 30;
	
    
    //OMP FOR Here
	for (unsigned y = 0; y < ImageHeight; y++){
		
        double c_im = MaxIm - y * Im_factor;
		
        for (unsigned x = 0; x < ImageWidth; ++x) {
            
            double c_re = MinRe + x * Re_factor;
            
            double Z_re = c_re, Z_im = c_im;
            int n=0;
            
            double Z_re2 = Z_re * Z_re, Z_im2 = Z_im * Z_im;
            while( Z_re2 + Z_im2 < 4 && n < MaxIterations) {
                Z_re2 = Z_re * Z_re, Z_im2 = Z_im * Z_im;
                Z_im = 2 * Z_re * Z_im + c_im;
                Z_re = Z_re2 - Z_im2 + c_re;
                n++;
            }
			
            values[y][x] = n ;
            if(1){
                
                cout<<"inside!: x:"<<x<<" Y: "<<y<<" n="<<n<<endl;
				//Send back the information
                
                //use drawing methods to draw x and y
				//values as an individual pixel
				//Z_re is x
				//Z_im is y
			}
		}
	}
    //print out the image
    FILE *f = fopen("/Users/cwieland/Desktop/out.ppm", "wb");
    fprintf(f, "P6\n%i %i 255\n",(int) ImageWidth , (int)ImageHeight);
    for (int y=0; y<ImageHeight; y++){
        for (int x=0; x<ImageWidth; x++)
        {
            if(values[y][x] >= MaxIterations/2){
                fputc((int)(values[y][x] * CON_FCTR) , f);
                fputc(255, f);
                fputc(255, f);
            }
            else{
                fputc((int)(values[y][x] * CON_FCTR) , f);
                fputc(0, f);
                fputc(0, f);
            }

        }
    }
    fclose(f);
    cout<<"DONE"<<endl;
    
	return 0;
}