#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define dim 1600

void matInit(long* arr, int row, int col){

  int i;
  int j;

    int M = col;
    for(i = 0; i < row; i++){
        for(j=0; j < col; j++)
	  arr[(M*i)+j] = (rand()%2);
    }
}


void printMat(long *arr, int row, int col){

    int M = col;
    for(int i = 0; i < row; i++){
        for(int j=0; j < col; j++)
            printf("%d ", arr[(M*i)+j]);

        printf("\n");
    }
}

long* matMult(long* arr1, long* arr2, int row1, int col1, int row2, int col2){

    if (col1 == row2){
        ;
    } else {
        printf("Matrix Dimensions do not Match\n");
        exit(0);
    }

    int M = col1;
    long* arrC = (long*)malloc(row1*col2*sizeof(long));

    for (int i = 0; i < row1; i++)
    {
        for (int j = 0; j < col1; j++)
        {
            arrC[(M*i)+j] = 0;
            for (int k = 0; k < col2; k++){
	      //printf("%d x %d\n",arr1[(M*i)+k],arr2[(M*k)+j]);
                arrC[(M*i)+j] += arr1[(M*i)+k] * arr2[(M*k)+j];
            }
        }
    }

    return(arrC);
}

int main(int argc, char* argv[]){
  
    //Matrix Dimensions
    int r1 = dim;
    int c1 = dim;
    int r2 = dim;
    int c2 = dim;

    //Timer
    double time;
    double globalTime;
    
    //Rank and Size of Processes
    int nump;
    int size;

    long *arrA = NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &nump);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    time = MPI_Wtime();

    //Create matrix A on master thread for distribution
    //to all processes
    if (nump == 0){
        arrA = (long*)malloc(r1*c1*sizeof(long));
        matInit(arrA, r1, c1);
	// printMat(arrA, r1, c1);
    }

    // printf("\n");

    //Create Matrix B on all processes
    long *arrB = (long*)malloc(r2*c2*sizeof(long));
  
    if (nump == 0){
        matInit(arrB, r2, c2);
	// printMat(arrB, r2, c2);
	// printf("\n");
    }

    long *arrC = NULL;
    if (nump == 0){
        arrC = (long*)malloc(r1*c1*sizeof(long));
    }

    if(size ==1){
        arrC = matMult(arrA, arrB, r1, c1, r2, c2);
	
	//End Timer
	globalTime  = MPI_Wtime() - time;
    
        //print out results
        printf("%d x %d Multiplication, Using %d Processes in %f time\n", r1, r1, size, globalTime);

	return 0;

    }
    
    //Send array B to all processes
    MPI_Bcast(arrB, r1*c1, MPI_LONG,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Decide work of matrix A distributed on each process
    //create buffer for array
    int split = r1/size;
    long* arrSplit = (long*)malloc(split*c1*sizeof(long));

    //scatter array A to all process
    MPI_Scatter(arrA, split*c1, MPI_LONG, arrSplit, split*c1, MPI_LONG, 0, MPI_COMM_WORLD);
        
    //Compute Matrix multiplication on all sub arrays
    long *arrC_sub;
    arrC_sub = matMult(arrSplit, arrB, split, c1, r2, c2);


    //Gather all sub arrays into resultant array C
    MPI_Gather(arrC_sub, split*c1, MPI_LONG, arrC, split*c1, MPI_LONG, 0, MPI_COMM_WORLD);

    time = MPI_Wtime() - time;

    MPI_Reduce(&time, &globalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Finalize();
    
    //print out results
    if (nump == 0){
      printf("%d x %d Multiplication, Using %d Processes in %f time\n", r1, r1, size, globalTime);
	// printMat(arrC, r1, c2);
    }

    return 0;
}
