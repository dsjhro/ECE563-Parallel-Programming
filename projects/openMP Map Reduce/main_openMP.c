#include "standardHeaders.h"
#include "splitUpWork.h"
#include "queue.h"
#include "readerThreadOMP.h"
#include "mapOpenmp.h"
#include <time.h>

omp_lock_t messLock;
omp_lock_t queLock;
omp_lock_t mapLock;
omp_lock_t redLock;
omp_lock_t writeLock;
omp_lock_t lk0;
omp_lock_t lk1;

int readers; 
int mappers; 
int reducers; 
int writers;

int mapsRecieved = 0;
int reduceDone = 0; 

int main(int argc, char* argv[]){
  
  
   //MPI Environment Variables
    int nodeRankNum = 0; //change to 0 for openMP
    int clusterSize = 1; //change to 1 for openMP

    omp_init_lock(&queLock);
    omp_init_lock(&mapLock);
    omp_init_lock(&redLock);  
    omp_init_lock(&lk0);
    omp_init_lock(&lk1);
    omp_init_lock(&messLock);  

    int rCt;
    int mCt;
    int uCt;
    int wCt;
    int totalRead;
    int toalMapped;
    
    //All threads to change working directory
    //****************************************************************
  
    //Working Directory string
    char cwd_old[PATH_MAX];
    char cwd_new[PATH_MAX];

    if (getcwd(cwd_old, sizeof(cwd_old)) != NULL) {
        ;
    } else {
        perror("getcwd() error");
    }

    //Change current directory to data folder
    strcpy(cwd_new, cwd_old);
    strcat(cwd_new, "/RawText");
    chdir(cwd_new);
   
    //Initialize Count variable and get files count in Data Folder
    int fileCount = 0;

    //initialize file queue
    //struct que* fileQue = initQue();
    struct que* fileQue = getNodeFileQue(1, 0);

    //Set number of threads
    int threadNum = 16;//omp_get_max_threads(); //10
    omp_set_num_threads(threadNum);

    int mapDone[10];
    int readDone[10];
    int redDone[10];

    for(int d = 0; d<10; d++)
      {

	mapDone[d] = 0;
	readDone[d] = 0;
	redDone[d] = 0;
      }

    #pragma omp barrier

    //Start Timer
    double time = omp_get_wtime();

    #pragma omp parallel
    {
      #pragma omp barrier
      #pragma omp master
      {

        readers = 5; //3 
        mappers = 5; //3
        reducers = 5; //4
	struct node* tempFile[readers];
	struct que* q[readers];
	struct mapChain* mapTables[mappers]; 
	struct mapChain* masterTable; 
	struct node* wordString[mappers];
	int capacity = 163841;//224737;//1583539;
	struct que* readQue;
	struct que* mapQue;
	struct que* redQue;
	struct que* writerQue;

        //struct node* tempFile[readers];
        //struct que* q[readers];
        //struct mapChain* mapTables[mappers]; 
        //struct node* wordString[mappers];
	char* tempPtr;

	for(int mt = 0; mt< readers; mt++){
                mapTables[mt] = initMapTable(capacity);
		q[mt] = initQue();
        }

	masterTable = initMapTable(capacity);
	readQue = initQue();
	mapQue = initQue();
	redQue = initQue();
	writerQue = initQue();

	int rCt = 0;
	int mCt = 0;
	int uCt = 0;
	int wCt = 0;

	int totalRead = 0;

	//Start Timer
	double timeTest = omp_get_wtime();

	//Reader Threads 
	for(int rr = 0; rr<readers; rr++)
	{
            #pragma omp task
	    {
	      
	      int fqCt;
	      fqCt = rCt;
	      
	      //omp_set_lock(&queLock);
	      tempFile[fqCt] = dequ(fileQue);
	      //omp_unset_lock(&queLock);

	      while(!queEmpty(fileQue))
		{
	  
	          if(tempFile[fqCt] != NULL)
		  {
	  //omp_set_lock(&queLock);
		      reader(q[fqCt], tempFile[fqCt]->str, lk0);
		      //printf("Count: %d, Thread: %d, %s\n", fqCt, omp_get_thread_num(),tempFile[fqCt]->str);
		      //omp_unset_lock(&queLock);
              	  }

		  //omp_set_lock(&queLock);
		  tempFile[fqCt] = dequ(fileQue);
		  //omp_unset_lock(&queLock);

		}
	      
	      //omp_set_lock(&queLock);
	        readDone[rCt] = 1;
		printf("Read Done: %d, Time: %f\n", rCt, omp_get_wtime() - timeTest);
		char readCount[2];
		sprintf(readCount, "%d", rCt);
	        enqu(readQue, readCount);
		//omp_unset_lock(&queLock);

	    }
	    //#pragma omp taskwait
	    rCt++;


         }
     
	
	//Mapper Threads 
	for(int mm = 0; mm<mappers; mm++)
	{
          
            #pragma omp task
	    {
	      int mqCt;
	      mqCt = mCt;

	      //Start Timer
	      double mapTime = omp_get_wtime();
	      
	      int reDone = -1;
	      char rtDone[2];
	      sprintf(rtDone, "%d", mqCt);
	      while (mqCt != reDone)
	      {
	             if(!queEmpty(readQue))
	             {
		           omp_set_lock(&queLock);
	                   struct node* readCheck = readQue->head;
			   while(readCheck != NULL)
			   {
	                           
	                         if(strcmp(readCheck->str, rtDone) == 0)
	                         {
	                            reDone = mqCt;
	                         }
			         readCheck = readCheck->next;
				   //printf("MapCT: %d",mapCt);
	                   }
			   omp_unset_lock(&queLock);
	             }

	      }
	      
	      
	      //omp_set_lock(&mapLock);
	      wordString[mqCt] = dequ(q[mqCt]);
	      //omp_unset_lock(&mapLock);
	      //printf("Count: %d, Thread: %d\n", mqCt, omp_get_thread_num());
	    
	      //map results
	      while((!queEmpty(q[mqCt])))
	      {
		  //#pragma omp critical
		  if(wordString[mqCt] != NULL)
		  {
	  //omp_set_lock(&mapLock);
		     mapper(mapTables[mqCt], wordString[mqCt]->str, capacity);
		     //         omp_unset_lock(&mapLock);
	          }

		  //omp_set_lock(&mapLock);
		     wordString[mqCt] = dequ(q[mqCt]);
		     // omp_unset_lock(&mapLock);

	       }

	      
	      omp_set_lock(&mapLock);
	      mapDone[mCt] = 1;
	      enqu(mapQue, "DONE");
	      //printf("Map Done:%d\n", mapDone[mCt]);
	      printf("Map Done: %d, Time: %f\n", mqCt,omp_get_wtime() - mapTime);
	      omp_unset_lock(&mapLock);
	      
	  }
	  mCt++;

	}


	if(threadNum > 1)
	{
	
	//Reducer Threads 
	for(int rd = 0; rd<reducers; rd++)
	  {
            #pragma omp task
	    {
	      int rqCt;
	      rqCt = uCt;
	      
	      int mapCt = 0;

	      while (mapCt < mappers)
	      {
	             if(!queEmpty(mapQue))
	             {

	                   struct node* tempCheck = mapQue->head;
			   mapCt = 0;
	                   omp_set_lock(&mapLock);
			   while(tempCheck != NULL)
			   {
	                           mapCt++;
			           tempCheck = tempCheck->next;
	                   }
			   omp_unset_lock(&mapLock);
	             }

	      }

	      //Start Timer
	      double redTime = omp_get_wtime();


	      for(int red = rqCt; red<capacity; red+=reducers)
	      {
		  for(int redmap = 0; redmap<mappers; redmap++)
		    {
		      	        struct mapNode* curChain = (struct mapNode*) mapTables[redmap][red].head;
		    		if (curChain == NULL) {
				  ; //skip
				} else {
				    while (curChain != NULL){
				        combineWords(masterTable, curChain->mapStr,red,curChain->wordCount);  
				        curChain = curChain->next;
				    }
				}



		    }
		  

	      }
	      omp_set_lock(&redLock);
	      redDone[uCt] = uCt;
	      char redCount[2];
	      sprintf(redCount, "%d", uCt);
	      enqu(redQue, redCount);
	      printf("Reduce Done:%d, Time: %f\n", uCt, omp_get_wtime() - redTime);
	      omp_unset_lock(&redLock);
	      
	    }
	    uCt++;
	    
	  }

	//Writer Threads
	for(int ww = 0; ww<reducers; ww++)
	{
            #pragma omp task
	    {
	      int wqCt;
	      wqCt = wCt;
	      int rDone = -1;
	      char rtDone[2];
	      sprintf(rtDone, "%d", wqCt);
	      while (wqCt != rDone )
	      {
	             if(!queEmpty(redQue))
	             {

	                   struct node* tempCheck = redQue->head;
	                   omp_set_lock(&redLock);
			   while(tempCheck != NULL)
			   {
	                           
	                         if(strcmp(tempCheck->str, rtDone) == 0)
	                         {
	                            rDone = wqCt;
	                         }
			         tempCheck = tempCheck->next;
				   //printf("MapCT: %d",mapCt);
	                   }
			   omp_unset_lock(&redLock);
	             }

	      }
	      
	     

	      strcat(cwd_old, "/Output");
	      chdir(cwd_old);


	      saveMapToFile(masterTable, capacity, wqCt, reducers);
	      omp_set_lock(&writeLock);
	      enqu(writerQue, "DONE");
	      printf("Write Done: %d, Time: %f\n", wqCt, omp_get_wtime() - timeTest);
	      omp_unset_lock(&writeLock);
	      
	  }
	  wCt++;  
      }
      } else {
	    strcat(cwd_old, "/Output");
	    chdir(cwd_old);
	    saveMapToFile(mapTables[0], capacity, 0, 1);
      }
      
      if(threadNum > 1)
	{
	      int wriCt = 0;

	      while (wriCt < reducers)
	      {
	             if(!queEmpty(writerQue))
	             {

	                   struct node* wriCheck = writerQue->head;
			   wriCt = 0;
	                   omp_set_lock(&writeLock);
			   while(wriCheck != NULL)
			   {
	                           wriCt++;
			           wriCheck = wriCheck->next;
				   //printf("MapCT: %d",mapCt);
	                   }
			   omp_unset_lock(&writeLock);
	             }

	      }
	}



	printf("Time: %f\n", omp_get_wtime() - timeTest);
	/*
	for(int mpc = 0; mpc<10; mpc++){
	      free(wordString[mpc]);
	      free(mapTables[mpc]);
	    }
	
      
	*/
	
    
      }
    #pragma omp barrier
	
	
 
	   

    
    }



//free(fileQue);
    




    
    return 0;

}

