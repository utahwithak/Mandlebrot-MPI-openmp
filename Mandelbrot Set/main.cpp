#include <iostream>
#include "mpi.h"
#include <sys/time.h>

#define CON_FCTR 255/30
using namespace std;

double When();


int main(int argc, char** argv) {
    
    int nproc, iproc;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    
    /*All for the sending / recieving */
    int* topLeft = (int*) malloc(2*sizeof(int));
    int* bottomRight = (int*) malloc(2*sizeof(int));
    
    
    
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    
    int Size = 1020;

	double startTime = When();
	int MaxIterations = 1000;
    int blocksize=Size/(nproc-1);
    
    if(iproc==0){
        fprintf(stderr,"\n%d:iproc\n", iproc);

        int* values = (int*)malloc(Size*Size*sizeof(int));
        //MASTER
        //COLUMNS

        for (int i = 1; i<nproc; i++) {
            //give them all a chunk
            topLeft[0] = (Size /((nproc-1) / 2) * (i-1));
            topLeft[1] = (Size /((nproc-1) / 2) * (i-1));
            cout<<topLeft[0]<<" tl sending "<<topLeft[1]<<"  "<<i<<endl;
            
            bottomRight[0] = (Size/(nproc/2)*i);
            bottomRight[1] = (Size/(nproc/2)*i);
            cout<<bottomRight[0]<<" BR sending "<<bottomRight[1]<<"  "<<i<<endl;

            MPI_Send(topLeft, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(bottomRight, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
            
        }
        fprintf(stderr,"\n%d:iproc Done Sending\n", iproc);

        //Now we get them all back
        int* dataBlock=(int*) malloc(blocksize*blocksize * sizeof(int));

        for (int i = 1; i<nproc; i++) {
            
            MPI_Recv(dataBlock, blocksize, MPI_INT,i , 0, MPI_COMM_WORLD, &status);
            fprintf(stderr,"%d: recieved:%i\n", iproc,i);
            #pragma omp parallel for num_threads(8)
            for (int y = 0; y < blocksize; y++) {
                for (int x = 0; x < blocksize; x++){
                    //put the values back inside the values array
                    //fprintf(stderr,"\n%d:x: %i Y: %i %i\n", iproc,x,y,dataBlock[blocksize*y+x]);

                    values[ (Size * ( y + (blocksize * (i-1)))) + (x + (blocksize * (i-1))) ] = dataBlock[blocksize*y+x];
                }
            }
          
            
        }
        fprintf(stderr,"\n%d:iproc Done recieving\n", iproc);

        //OUTPUT
        FILE *f = fopen("/Users/cwieland/Desktop/out.ppm", "wb");
        fprintf(f, "P6\n%i %i 255\n",(int) Size , (int)Size);
        for (int y=0; y<Size; y++){
            for (int x=0; x<Size; x++)
            {
                if(values[Size * y + x] >= MaxIterations/2){
                    fputc((int)(values[Size * y + x] * CON_FCTR) , f);
                    fputc(255, f);
                    fputc(255, f);
                }
                else{
                    fputc((int)(values[Size * y + x] * CON_FCTR) , f);
                    fputc(0, f);
                    fputc(0, f);
                }
            }
        }
        fclose(f);
        
    }
    else{
        fprintf(stderr,"\n%d:iproc\n", iproc);

        int* dataBlock=(int*) malloc(blocksize*blocksize * sizeof(int));
        //WORKER
        //Get our Chunk
        MPI_Recv( topLeft, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv( bottomRight, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        fprintf( stderr,"\n%d:iproc Done Receiving: TLX: %i TLY:%i BRX:%i BRY:%i\n", iproc,topLeft[0],topLeft[1],bottomRight[0],bottomRight[1]);

        int offy=0;
        int offx=0;
        #pragma omp parallel for num_threads(8)
        for (int y = topLeft[1]; y < bottomRight[1]; y++){
            double MinRe		=	-2.0;
            double MaxRe		=	1.0;
            double MinIm		=	-1.2;
            double MaxIm		=	MinIm + (MaxRe - MinRe) * Size/Size;
            double Re_factor	=	(MaxRe - MinRe) / (Size - 1);
            double Im_factor	=	(MaxIm - MinIm) / (Size - 1);
            
            double c_im = MaxIm - y * Im_factor;
            
            for (int x = topLeft[0]; x < bottomRight[0]; ++x, offx++) {
                
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
                dataBlock[ (blocksize * offy) + offx] = n;
                //cout<<"inside!: x:"<<x<<" Y: "<<y<<" n="<<n<<endl;
                
            }
            offy++;
        }
        //send our chunk back
        MPI_Send(dataBlock, Size/(nproc-1) , MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    cout<<"DONE "<<iproc<<"  "<<(When()-startTime)<<endl;

    MPI_Finalize();

	return 0;
}


/* Return the correct time in seconds, using a double precision number.       */
double When()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

