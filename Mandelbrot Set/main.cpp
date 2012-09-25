#include <iostream>
#include "mpi.h"
#include <sys/time.h>
#include "math.h"
#include <omp.h>
#define CON_FCTR 255/30
using namespace std;

double When();
//$PBS_NODEFILE
//uniq($PBS_NODEFILE)>f
//mpirun -np nodes -machinefile f 

int main(int argc, char** argv) {
    
    int nproc, iproc;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    
    /*All for the sending / recieving */
    
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    
    int Size = 4096;
    
	double startTime = When();
	int MaxIterations = 1000;
    
    //int blocksize= Size / ((nproc-1) / ((int)sqrt((nproc-1))))  ;
    int chunkSize=32;

    //BLOCK DATA    
    
    if(iproc == 0){
        int* values = (int*)malloc(2*Size*Size*sizeof(int));
        int curSection=0;
        int tid=0;
        /* Fork a team of threads with each thread having a private tid variable */
        #pragma omp parallel private(tid) shared(values) num_threads(nproc-1)
        {
            int* dataBlock = (int*) malloc(2*chunkSize*chunkSize * sizeof(int));
            int* topLeft = (int*) malloc(2*sizeof(int));

            /* Obtain thread id */
            tid = omp_get_thread_num();
            //fprintf(stderr,"%d: tid=%i\n", iproc,tid);
            while (curSection < Size) {
                #pragma omp critical
                {
                    topLeft[0] = (curSection*chunkSize) % (Size)  ;
                    topLeft[1] = (curSection/chunkSize) * chunkSize;
                    curSection++;
                    MPI_Send(topLeft, 2, MPI_INT, tid+1, 0, MPI_COMM_WORLD);

                 }
                    //fprintf(stderr,"%d: sending:%i %i   >>>to iproc:%i  curSection:%i \n", iproc,topLeft[0],topLeft[1], tid+1,curSection);
                    //cout<<"sending: "<<topLeft[0]<<"  "<<topLeft[1]<<"  >>>tid: "<<tid<<" curSection: "<<curSection<<endl;
                #pragma omp critical
                {
                    MPI_Recv(dataBlock, chunkSize*chunkSize, MPI_INT, tid+1 , 0, MPI_COMM_WORLD, &status);
                }
                int count =0;
                int offy=0;
                int offx=0;
                for (int y = topLeft[1]; y < chunkSize + topLeft[1]; y++) {
                    offx=0;
                    for (int x = topLeft[0]; x < chunkSize + topLeft[0]; x++){ 
                        //fprintf(stderr,"x=%iy=%i offx:%ioffy%i\n",x,y,offx,offy);
                        values[Size*y + x] = dataBlock[(chunkSize*offy)+offx];
                        if(values[Size*y + x] != 0){

                            count++;
                        }
                    
                        offx++;
                    }
                    offy++;
                }
                //if(tid==0)
                //fprintf(stderr,"%d: %i->Count= %i\n", iproc,tid+1,count);

                //fprintf(stderr,"%d: %i end of while\n", iproc,tid);

            
            }
            topLeft[0] = -1;
            topLeft[1] = -1;
            #pragma omp critical
            {
                //cout<<"0:KILLING!:"<<topLeft[0]<<topLeft[1]<<"  >>>"<<tid<<endl;
                MPI_Send(topLeft, 2, MPI_INT, tid+1, 0, MPI_COMM_WORLD);
            }
        }//end pragma omp
        cout<<"FINISHED!! "<<iproc<<"  "<<(When()-startTime)<<endl;
        /*FILE *f = fopen("/~/out.ppm", "wb");
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
         */


    }
    else{
        int* topLeft = (int*) malloc(2*sizeof(int));
        int* dataBlock = (int*) malloc(2*chunkSize*chunkSize * sizeof(int));

        //WORKER
        //Get our Chunk

        MPI_Recv( topLeft, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        while (topLeft[0]>=0) {
            // fprintf( stderr,"\n%d:iproc Done Receiving: TLX: %i TLY:%i BRX:%i BRY:%i\n", iproc,topLeft[0],topLeft[1],topLeft[0]+chunkSize,topLeft[1]+chunkSize);
            
            int offy=0;
            int offx = 0;
            int count =0;
#pragma omp parallel for 
            for (int y = topLeft[1]; y < topLeft[1]+chunkSize; y++){
                double MinRe		=	-2.0;
                double MaxRe		=	1.0;
                double MinIm		=	-1.2;
                double MaxIm		=	MinIm + (MaxRe - MinRe) * Size/Size;
                double Re_factor	=	(MaxRe - MinRe) / (Size - 1);
                double Im_factor	=	(MaxIm - MinIm) / (Size - 1);
                
                double c_im = MaxIm - y * Im_factor;
                offx=0;
                for (int x = topLeft[0]; x < topLeft[0] + chunkSize ; x++) {

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
                    dataBlock[ (chunkSize * offy) + offx ] = n;
                    
                    offx++;
                    
                }
                
                offy++;
                //fprintf(stderr,"%d:%i and %i\n", iproc,y,topLeft[1]+blocksize);
                
            }
            //send our chunk back
            //if(iproc==1)
            //fprintf(stderr,"************************************%d:Sending Back count: %i\n", iproc,count);
            MPI_Send(dataBlock, chunkSize*chunkSize , MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(topLeft, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        }
              
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

