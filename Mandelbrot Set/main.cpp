#include <iostream>
#include "mpi.h"
#include <sys/time.h>
#include "math.h"
#define CON_FCTR 255/30
using namespace std;

double When();


int main(int argc, char** argv) {
    
    int nproc, iproc;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    
    /*All for the sending / recieving */
    int* topLeft = (int*) malloc(2*sizeof(int));

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    
    int Size = 8192;

	double startTime = When();
	int MaxIterations = 30;
    
    int blocksize= Size / ((nproc-1) / ((int)sqrt((nproc-1))))  ;
    
    //BLOCK DATA    
    int* dataBlock = (int*) malloc(blocksize*blocksize * sizeof(int));

    if(iproc == 0){

        int* values = (int*)malloc(Size*Size*sizeof(int));
        //MASTER
        //COLUMNS
        /*  1 2 
            3 4
         
         */

        int chunkSize=((nproc-1) / 2);
        for (int i = 1; i<nproc; i++) {
            //give them all a chunk
            //for 4 it should be
            /*
                x= 0, 510, 0, 510
                y= 0, 0, 510, 510
             */
            
            
            topLeft[0] = ((i-1) % chunkSize) * Size/chunkSize;
            topLeft[1] = ((i-1)/chunkSize) * Size/chunkSize ;

            
            cout<<topLeft[0]<<" tl sending "<<topLeft[1]<<"  >>>"<<i<<endl;
            MPI_Send(topLeft, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
            
        }
        fprintf(stderr,"\n%d:iproc Done Sending\n", iproc);


        for (int i = 1; i<nproc; i++) {
            MPI_Recv(dataBlock, blocksize*blocksize, MPI_INT, i , 0, MPI_COMM_WORLD, &status);
            int count =0;
            for (int y = 0; y < blocksize; y++) {
                for (int x = 0; x < blocksize; x++){
                    if(dataBlock[(blocksize*y) + x ] != 0){
                        //fprintf(stderr,"%d: Setting %i %i to:%i\n", iproc,x,y,dataBlock[blocksize*y+x]);
                        count++;
                    }

                    values[ (Size * ( y + (((i-1)/chunkSize) * Size/chunkSize))) +  ((((i-1) % chunkSize) * Size/chunkSize)+x) ] = dataBlock[blocksize*y+x];
                }
            }
            fprintf(stderr,"%d: %i->Count= %i\n", iproc,i,count);

          
            
        }
        fprintf(stderr,"\n%d:iproc Done recieving\n", iproc);
        cout<<"DONE "<<iproc<<"  "<<(When()-startTime)<<endl;

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
        cout<<"DONE Creating File"<<endl;

        
    }
    else{

        //WORKER
        //Get our Chunk
        MPI_Recv( topLeft, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        //fprintf( stderr,"\n%d:iproc Done Receiving: TLX: %i TLY:%i BRX:%i BRY:%i\n", iproc,topLeft[0],topLeft[1],topLeft[0]+blocksize,topLeft[1]+blocksize);

        int offy=0;
        int offx = 0;
        int count =0;
        for (int y = topLeft[1]; y < topLeft[1]+blocksize; y++){
           double MinRe		=	-2.0;
           double MaxRe		=	1.0;
           double MinIm		=	-1.2;
           double MaxIm		=	MinIm + (MaxRe - MinRe) * Size/Size;
           double Re_factor	=	(MaxRe - MinRe) / (Size - 1);
           double Im_factor	=	(MaxIm - MinIm) / (Size - 1);
            
            double c_im = MaxIm - y * Im_factor;
            offx=0;
            for (int x = topLeft[0]; x < topLeft[0] + blocksize ; x++) {
                
                //if(iproc==1)
                    // fprintf(stderr,"\t%d:%i\n", iproc,x);

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
            
                if(n!=0){
                    //fprintf(stderr,"%d:SETTING DATA BLOCK TO %i at: %i %i, offx %i, offy %i \n", iproc,n,x,y,offx,offy);
                    count++;
                }
                dataBlock[ (blocksize * offy) + offx ] = n;

                offx++;
            
               
            }
             
            offy++;
            //fprintf(stderr,"%d:%i and %i\n", iproc,y,topLeft[1]+blocksize);

        }
        //send our chunk back
        fprintf(stderr,"************************************%d:Sending Back count: %i\n", iproc,count);

        MPI_Send(dataBlock, blocksize*blocksize , MPI_INT, 0, 0, MPI_COMM_WORLD);
        cout<<"DONE "<<iproc<<"  "<<(When()-startTime)<<endl;
    }
    //free(dataBlock);

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

