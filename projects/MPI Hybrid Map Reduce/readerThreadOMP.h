#include "standardHeaders.h"

void reader(struct que* q, char* filename, omp_lock_t writelock){
	
  //omp_lock_t writelock;
  //omp_init_lock(&writelock);
        
        FILE* filePtr;
	int fSize = fileSize(filename);
        
	char* fileAllStr = malloc(fSize);


        filePtr = fopen(filename, "rb");
        if (filePtr == NULL){
            printf("Could not open file: %s",filename);
        }
	
	//printf("FileName: %s\n",filename);
        fread(fileAllStr, 1, fSize, filePtr);
	enqu(q,fileAllStr);
	
	fclose(filePtr);
	
        /*
	int readSize = 200;
	int Lines = 5;
	char line[readSize+1];
	char* buffer = malloc(Lines*readSize+Lines);
        char* tempLine = NULL;
	
        filePtr = fopen(filename, "rb");
        if (filePtr == NULL){
            printf("Could not open file: %s",filename);
        }
	
        tempLine = fgets(line, readSize, filePtr);
	int lineCount = 0;

	//printf("Num: %d", omp_get_thread_num());

	while (tempLine != NULL){
	  //printf("%s", tempLine);
	          omp_set_lock(&writelock);            
                  enqu(q,tempLine);
                  omp_unset_lock(&writelock);
	  
	  /*
	   
	    strcat(buffer,tempLine);
	    lineCount++;
	    if(lineCount = Lines){
	          printf("%s\n", buffer);
	          //add to queue
	          omp_set_lock(&writelock);            
                  enqu(q,buffer);
                  omp_unset_lock(&writelock);
		  memset(&buffer[0], 0, sizeof(buffer));
		  lineCount = 0;
	    }
	  
	*/
	/*
            tempLine = fgets(line, readSize, filePtr);

        }
  
        fclose(filePtr);
*/
	return;
	
	
}
