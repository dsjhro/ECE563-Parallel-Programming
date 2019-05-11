#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include "mkl.h"

#define dim 1600


void matInit(double* arr,long row,long col){

    int M = col;
    int i, j;
    for(i = 0; i < row; i++){
        for(j=0; j < col; j++)
	  arr[(M*i)+j] = (double)(rand()%2)+.2*i;
    }
}


void printMat(double *arr,long row,long col){

    int M = col;
    int i, j;
    for(i = 0; i < row; i++){
        for(j = 0; j < col; j++)
            printf("%f ", arr[(M*i)+j]);

        printf("\n");
	fflush(stdout);
    }
}


double* matMult(double* arr1,double* arr2,long row1,long col1,long row2,long col2){

    if (col1 == row2){
        ;
    } else {
        printf("Matrix Dimensions do not Match\n");
        exit(0);
    }

    int i, j, k;

    int M = col1;
    double* arrC = (double *)mkl_malloc(row1*col2*sizeof(double), 64);

    for(i = 0; i < row1; i++)
    {
        for(j = 0; j < col1; j++)
        {
            arrC[(M*i)+j] = 0;
            for(k = 0; k < col2; k++){
                arrC[(M*i)+j] += arr1[(M*i)+k] * arr2[(M*k)+j];
            }
        }
    }

    return(arrC);
}


double* matAdd(double* arr1,double* arr2,long xDim,long yDim){

    int M = xDim;
    double* arrC = (double*)mkl_malloc(xDim*yDim*sizeof(double), 64);
    int i = 0;
    int j = 0;

    for(i = 0; i < xDim; i++){
        arrC[(M*i)+j] = 0;
        for(j = 0; j < yDim; j++){
            arrC[(M*i)+j] = arr1[(M*i)+j] + arr2[(M*i)+j];
        }
    }

    return(arrC);
}

int main(int argc, char* argv[]){

    //Dimensions
    long r1 = dim, c1 = dim, r2 = dim, c2 = dim;

    //Timer
    double time;
    double globalTime;

    //Rank and Size of Processes
    int nump;
    int size;

    //Initialize Matrices
    double *arrA = NULL;
    double *arrB = NULL;

    //MPI Initilization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &nump);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Find Sqrt of processors for splititng up work
    int sqP = (int) sqrt(size);

    //Start Timer
    time = MPI_Wtime();

    //Create matrix A and B on master thread for distribution
    //to all processes
    if (nump == 0){
      arrA = (double*)mkl_malloc(r1*c1*sizeof(double),64);
      arrB = (double*)mkl_malloc(r2*c2*sizeof(double), 64);
        matInit(arrA, r1, c1);
        matInit(arrB, r2, c2);
        //printMat(arrA, r1, c1);
	// printf("\n\n");
	// printMat(arrB, r1, c1);
    }

    //Create blocks of array for distribution
    long xDim = c1/((int) sqrt(size));
    long yDim = r1/((int) sqrt(size));
    MPI_Datatype block;
    MPI_Datatype blocktype;

    //Creates sub array for chunks of arrays
    double * arrA_sub = (double*)mkl_malloc(xDim*yDim*sizeof(double),64);
    double* arrB_sub = (double*)mkl_malloc(xDim*yDim*sizeof(double), 64);

    //Create Block Vector to break up Array
    MPI_Type_vector(xDim, yDim, r1, MPI_DOUBLE, &block);
    MPI_Type_commit(&block);
    MPI_Type_create_resized(block, 0, xDim*sizeof(double), &blocktype);
    MPI_Type_commit(&blocktype);

    //Determine Scatter Stats and scatter array to all processes
    int disp[size];
    int scount[size];
    int rcount = ((r1*r1)/size);

    //Determine Display Stats, scount and rcount
    int rCt = 0;
    int eL = 0;
    int rowLim = (int) sqrt(size);
    int disect = (int) sqrt(size);
    int i;

    for(i = 0; i<size; i++){

      scount[i] = 1;
      disp[i] = rCt;
      rCt = rCt + 1;

      if(rCt == rowLim){
	eL = eL + (xDim*disect);
        rCt = eL;
	rowLim = rCt + disect;
      }
    }

    //Scatters blocks of A and B to all processes
    MPI_Scatterv(arrA,scount,disp,blocktype,arrA_sub,rcount,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatterv(arrB,scount,disp,blocktype,arrB_sub,rcount,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    //Create C  sub matrix and initialize to zero
    double *arrC_sub = (double*)mkl_malloc(xDim*yDim*sizeof(double),64);
    long N = xDim;
    int f, g;

    for(f = 0; f < xDim; f++){
	for(g = 0; g < yDim; g++){
	  arrC_sub[f*N + g] = 0.0;
	}
    }

    //Create Array C buffer
    double *arrC = NULL;
    if (nump == 0){
      arrC = (double*)mkl_malloc(r1*c1*sizeof(double), 64);
    }
   
    if (size == 1){
      //arrC_sub = matMult(arrA_sub, arrB_sub, xDim, yDim, xDim, yDim);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                xDim, xDim, xDim, 1.0, arrA_sub, xDim, arrB_sub, xDim, 0.0, arrC_sub, xDim);
	globalTime = MPI_Wtime() - time;
	//Print out Results
        //printf("\n");
        printf("%d x %d Multiplication, Using %d Processes in %f time\n", r1, r1, size, globalTime); 
        //printMat(arrC_sub,r1,c1);
        //printf("\n");
        fflush(stdout);
	return 0;

    } else {

      //Send and Recieve Coordinates
      int sendA=0, recA=0, sendB=0, recB=0;
      MPI_Status status;

      //Split Communicators into rows and columns for easier communication
      //initialize all variables needed to implement a COMM split command

      int rowNump = 0, colNump = 0;
      int rowSize = 0, colSize = 0;
      int colorRow = 0, colorCol = 0;
      MPI_Comm col, row;

      //Group rows into own unique communicator
      colorRow = (nump/sqP);
      MPI_Comm_split(MPI_COMM_WORLD, colorRow, nump, &row);

      //Determine New Row group rank and size
      MPI_Comm_size(row, &rowSize);
      MPI_Comm_rank(row, &rowNump);

      //Group columns into own unique communicator
      colorCol = nump%sqP;
      MPI_Comm_split(MPI_COMM_WORLD, colorCol, nump, &col);

      //Determine New Column group rank and size
      MPI_Comm_size(col, &colSize);
      MPI_Comm_rank(col, &colNump);

      //Determine were to send A sub matrices of row group
      if (rowNump == 0){
        sendA = (sqP-1);
      } else {
        sendA = rowNump - 1;
      }
      if (rowNump == (sqP-1)){
        recA = 0;
      } else {
        recA = rowNump + 1;
      }

      //Determine were to send B sub matrices of column group
      if (colNump == 0){
        sendB = (sqP-1);
      } else {
        sendB = colNump - 1;
      }
      if (colNump == (sqP-1)){
        recB = 0;
      } else {
        recB = colNump + 1;
      }
      
      //Initialize new sub matrix to hold sub resultant matrix
      double * arrA_sub_new = (double*)mkl_malloc(xDim*yDim*sizeof(double), 64);
      double * arrB_sub_new = (double*)mkl_malloc(xDim*yDim*sizeof(double), 64);
      double * arrC_sub_new = (double*)mkl_malloc(xDim*yDim*sizeof(double), 64);
      
      //Create inital skewed array
      //*******************************************************
      
      int left;
      
      //Original skew of row i; i places left
      for(left = 0; left < colorRow; left++){
	
	  MPI_Request reqA[2];
	  MPI_Status stA[2];

	  //Send Sub Matrix (A) one coordinate left (With Wrap Around) 
	  MPI_Irecv(arrA_sub_new, xDim*yDim, MPI_DOUBLE, recA, 1, row, &reqA[0]);
	  MPI_Isend(arrA_sub, xDim*yDim, MPI_DOUBLE, sendA, 1, row, &reqA[1]);
          
	  MPI_Waitall(2, reqA, stA);
	  memcpy(arrA_sub, arrA_sub_new, xDim*yDim*sizeof(double));

      }

      MPI_Barrier(MPI_COMM_WORLD);
      
      int up;

      //Original skew of col j; j places up
      for(up = 0; up < colorCol; up++){

	  MPI_Request reqB[2];
	  MPI_Status stB[2];

	  //Send Sub Matrix (B) one coordinate up (With Wrap Around) 
	  MPI_Irecv(arrB_sub_new, xDim*yDim, MPI_LONG, recB, 1, col, &reqB[0]);
	  MPI_Isend(arrB_sub, xDim*yDim, MPI_LONG, sendB, 1, col, &reqB[1]);
          
	  MPI_Waitall(2, reqB, stB);
	  memcpy(arrB_sub, arrB_sub_new, xDim*yDim*sizeof(double));	
      }

      MPI_Barrier(MPI_COMM_WORLD);
      
      //Multiply Initial Array before shift
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                xDim, xDim, xDim, 1.0, arrA_sub, xDim, arrB_sub, xDim, 0.0, arrC_sub, xDim); 

      MPI_Barrier(MPI_COMM_WORLD);
      
      int s;

      //Shift and calculate all sub matrices 
      for(s = 0; s < ((int)sqrt(size)-1); s++){
	  
	  MPI_Request reqAB[4];
	  MPI_Status stAB[4];
	  
	  //Send Sub Matrix (A) one coordinate left (With Wrap Around) 
	  MPI_Irecv(arrA_sub_new, xDim*yDim, MPI_LONG, recA, s, row, &reqAB[0]);
	  MPI_Isend(arrA_sub, xDim*yDim, MPI_LONG, sendA, s, row, &reqAB[1]);
	  
	  MPI_Barrier(MPI_COMM_WORLD);

	  //Send Sub Matrix (B) one coordinate up (With Wrap Around) 
	  MPI_Irecv(arrB_sub_new, xDim*yDim, MPI_LONG, recB, s+1, col, &reqAB[2]);
	  MPI_Isend(arrB_sub, xDim*yDim, MPI_LONG, sendB, s+1, col, &reqAB[3]);

	  MPI_Waitall(4, reqAB, stAB);
	  MPI_Barrier(MPI_COMM_WORLD);

	  //replace old array with newly recieved array
	  memcpy(arrA_sub, arrA_sub_new, xDim*yDim*sizeof(int));
	  memcpy(arrB_sub, arrB_sub_new, xDim*yDim*sizeof(int));

	  MPI_Barrier(MPI_COMM_WORLD);
	  
	  //Multiply new sub arrays and them to existing C array
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                xDim, xDim, xDim, 1.0, arrA_sub, xDim, arrB_sub, xDim, 0.0, arrC_sub_new, xDim); 
	  arrC_sub = matAdd(arrC_sub, arrC_sub_new, xDim, yDim);

	  MPI_Barrier(MPI_COMM_WORLD);

	}

      MPI_Barrier(MPI_COMM_WORLD);
    
        //Gather all sub arrays into resultant array C
        MPI_Gatherv(arrC_sub, xDim*yDim, MPI_DOUBLE, arrC,scount,disp, blocktype, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	
	//Stop Timer
	time = MPI_Wtime() - time;

	//Find Max time for global timer
	MPI_Reduce(&time, &globalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        //Print out Results
        if (nump == 0){
          printf("%d x %d Multiplication, Using %d Processes in %f time\n", r1, r1, size, globalTime);
	  //printf("\n");
          //printMat(arrC,r1,c1);
          //printf("\n");
          fflush(stdout);
        }



       }


    return 0;

}

