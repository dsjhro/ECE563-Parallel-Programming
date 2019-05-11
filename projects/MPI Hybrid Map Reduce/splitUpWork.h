#ifndef SPLITUPWORK_H_INCLUDED
#define SPLITUPWORK_H_INCLUDED

#include "queue.h"

void getFileList(char fileList[][PATH_MAX]){
   
    //Working Directory string
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        ;
    } else {
        perror("getcwd() error");
    }

    //Change current directory to data folder
    strcat(cwd, "\\RawText");
    chdir(cwd);

    int ct = 0;

    //Directory Pointer
    DIR *directory;
    struct dirent *dir;

    //Open Current directory path
    directory = opendir(".");

    //Parse Directory and save to file array
    if (directory)
    {
        while ((dir = readdir(directory)) != NULL)
        {
            if (dir->d_name[0] != '.'){
                strcpy(fileList[ct],dir->d_name);
                //printf("File: %s\n", dir->d_name);
                ct++;
            }
        }
        closedir(directory);
    }

}

int fileSize(char file[PATH_MAX]){

    int siz = 0;

    FILE *fileP;
    fileP = fopen(file, "rb");

    fseek(fileP, 0L, SEEK_END);
    siz = ftell(fileP);
    rewind(fileP);

    fclose(fileP);

    return(siz);

}

int getWorkSize(char fileList[][PATH_MAX], int fileCount){

    int totalSize = 0;
    int i;

    for(i=0; i< fileCount; i++){
        totalSize += fileSize(fileList[i]);
    }

    return(totalSize);
}


int getFileCount(){

    //Working Directory string
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        ;
    } else {
        perror("getcwd() error");
        return 1;
    }

    //Change current directory to data folder
    strcat(cwd, "\\RawText");
    chdir(cwd);

    //File Count variable
    int fileCount = 0;

    //Count files in data folder directory
    DIR *dr;
    struct dirent *dir;

    //Open Working Directory
    dr = opendir(".");

    //Parse Directory
    while ((dir = readdir(dr)) != NULL)
    {
        if (dir->d_name[0] !='.'){
            //Count Files
            fileCount++;
        }
    }
    //Close Directory
    closedir(dr);

    return(fileCount);
}


struct que* getNodeFileQue(int clusterSize, int nodeRankNum)
{

    //Initialize Count variable and get files count in Data Folder
    int fileCount = 0;
    fileCount = getFileCount();

    //Init Array of File Path Names and fill with list of all files
    //in the data folder
    char fileList[fileCount][PATH_MAX];
    getFileList(fileList);

    int fileAmount = (int)ceil(fileCount/((float)clusterSize));
    int nodeFileList[fileAmount];   
   
    for(int in=0; in<fileAmount; in++)
    {
      nodeFileList[in] = -1;
    }

    int nodeFileCount = 0;

    for(int f=nodeRankNum; f<fileCount; f+=clusterSize)
    {
      nodeFileList[nodeFileCount] = f;
      nodeFileCount++;

    }

    struct que* fileQue = initQue();
  
    int j = 0;
    for(j=0; j<fileAmount; j++)
    {
      //printf("File: %s\n",fileList[j]);
      if(nodeFileList[j] != -1)
	{
       enqu(fileQue, fileList[nodeFileList[j]]);
	}
    } 
    
  
    return(fileQue);

}


#endif // SPLITUPWORK_H_INCLUDED
