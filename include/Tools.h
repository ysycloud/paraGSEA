/* File:     Tools.h
 * Purpose:  Header file for Tools.c.
 */
#ifndef _MY_Tools_H_
#define _MY_Tools_H_

#include "IO.h"


#include "sys/sysinfo.h"

void split_data(int size, int n, int rank, int* begin, int* end, int* local_n); 
void getTriples(int local_P, int genelen, int siglen, int profilenum, int linelen, int begin, int end,  char *file, struct Profile_triple * triples);
void getPartTriples(int genelen, int siglen, int profilenum, int linelen, int begin, int end,  char *file, struct Profile_triple * triples);
//file:begin1->begin1+len => triples: begin2->begin2+len  
void getFreeTriples(int genelen, int siglen, int profilenum, int linelen, int begin1, int begin2, int len, char *file, struct Profile_triple * triples);

int cmpset(int *set1,int *set2,int n);
int isInSet(int **set1,int *set2,int n,int iter);
unsigned long memoryAvailable(int base);

#endif
