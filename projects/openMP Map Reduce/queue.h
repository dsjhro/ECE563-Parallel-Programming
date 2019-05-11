#ifndef QUEUE_H_INCLUDED
#define QUEUE_H_INCLUDED

#include "standardHeaders.h"

struct node {
    char* str;
    struct node* next;
};

struct que{
    struct node* head;
    struct node* tail;
};

struct node* addWord(char* word){
    struct node* curNode = (struct node*)malloc(sizeof(struct node));
    curNode->str = malloc(strlen(word)+1);
    strcpy(curNode->str, word);
    curNode->next = NULL;
    return(curNode);
};

struct que* initQue(){
    struct que *q = (struct que*)malloc(sizeof(struct que));
    q->head = q->tail = NULL;
    return(q);
};

void enqu(struct que* q, char* word){

    //Add word to work queue
    struct node* curNode = addWord(word);

    //If it is the first word in queue
    //Point back to itself
    if (q->tail == NULL){
       q->head = q->tail = curNode;
       return;
    }

    //Else adjust pointer to point to next node
    //And make new tail equal to currently added node
    q->tail->next = curNode;
    q->tail = curNode;
}

// Function to remove a key from given queue q
struct node *dequ(struct que *q)
{

    omp_lock_t queLock;
    omp_init_lock(&queLock);
    omp_set_lock(&queLock);

    if (q->head == NULL){
       return NULL;
    }

    struct node *curHeadNode = q->head;
    q->head = q->head->next;

    if (q->head == NULL){
       q->tail = NULL;
    }
    
    omp_unset_lock(&queLock);

    return(curHeadNode);
}


int queEmpty(struct que *q)
{

  if(q->head == NULL && q->tail == NULL)
  {
      return(1);
  } else {
      return(0);
  }
  
}

#endif // QUEUE_H_INCLUDED
