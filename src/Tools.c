#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include "Tools.h"

void split_data(int size, int n, int my_rank, int* begin, int* end, int* local_n)
{
	*local_n = size / n;  
	int leave = size % n;
	if(my_rank < leave){
		(*local_n)++;   
		*begin = my_rank*(*local_n);
	}else{	
		*begin = my_rank*(*local_n) + leave; 
	}	 
	*end = *begin + *local_n;
}


//read the begin to end line part of profile file and get their triples
//file:begin->end => triples: 0->local_P  
void getTriples(int local_P, int genelen, int siglen, int profilenum, int linelen, int begin, int end,  char *file, struct Profile_triple * triples)
{
	int i;
	//allocate the temp memory
	short **profileSet = (short **)malloc(local_P*sizeof(short *));
	for(i=0;i<local_P;i++)
		profileSet[i] = (short *)malloc(genelen*sizeof(short));	

	//read file and get the proper data
	ReadFile_new(file, linelen, begin, end, profilenum, genelen, profileSet);
	
	//get the triple for every profile	
	for(i=0;i<local_P;i++)
		triples[i] = getTriple(profileSet[i], genelen, siglen);
	
	//free the temp memory
	for(i=0;i<local_P;i++)
		free(profileSet[i]);
	free(profileSet);	
}

//read the begin to end line part of profile file and get their triples in proper part of the  triples vector
//file:begin->end => triples: begin->end  
void getPartTriples(int genelen, int siglen, int profilenum, int linelen, int begin, int end,  char *file, struct Profile_triple * triples)
{
	int i;
	int len = end - begin;
	//allocate the temp memory
	short **profileSet = (short **)malloc(len*sizeof(short *));
	for(i=0;i<len;i++)
		profileSet[i] = (short *)malloc(genelen*sizeof(short));	
	
	//read file and get the proper data
	ReadFile_new(file, linelen, begin, end, profilenum, genelen, profileSet);
	
	//get the triple for every profile	
	for(i=begin;i<end;i++)
		triples[i] = getTriple(profileSet[i-begin], genelen, siglen);
	
	//free the temp memory
	for(i=0;i<len;i++)
		free(profileSet[i]);
	free(profileSet);		
}

//file:begin1->begin1+len => triples: begin2->begin2+len  
void getFreeTriples(int genelen, int siglen, int profilenum, int linelen, int begin1, int begin2, int len, char *file, struct Profile_triple * triples)
{
	int i;
	//allocate the temp memory
	short **profileSet = (short **)malloc(len*sizeof(short *));
	for(i=0;i<len;i++)
		profileSet[i] = (short *)malloc(genelen*sizeof(short));	
	
	//read file and get the proper data
	ReadFile_new(file, linelen, begin1, begin1 + len, profilenum, genelen, profileSet);
	
	//get the triple for every profile	
	for(i=0;i<len;i++)
		triples[begin2+i] = getTriple(profileSet[i], genelen, siglen);
	
	//free the temp memory
	for(i=0;i<len;i++)
		free(profileSet[i]);
	free(profileSet);		
}

/*
//read the begin to end line part of profile file and get their triples  
void getTriples(int local_P, int genelen, int siglen, int profilenum, int linelen, int begin, int end,  char *file, struct Profile_triple * triples)
{
	int i;
	//allocate the temp memory
	short **profileSet = (short **)malloc(local_P*sizeof(short *));
	for(i=0;i<local_P;i++)
		profileSet[i] = (short *)malloc(genelen*sizeof(short));	
	short **tmp_profiles = (short **)malloc(profilenum*sizeof(short *));
	for(i=0;i<profilenum;i++)
		tmp_profiles[i] = (short *)malloc(genelen*sizeof(short));
	
	//read file and get the proper data
	ReadFile(file, linelen, begin, end, profilenum, genelen, tmp_profiles);
	for(i=0; i<local_P; i++)	
		memcpy(profileSet[i],tmp_profiles[i+begin],genelen*sizeof(short));
	
	//get the triple for every profile	
	for(i=0;i<local_P;i++)
		triples[i] = getTriple(profileSet[i], genelen, siglen);
	
	//free the temp memory
	for(i=0;i<local_P;i++)
		free(profileSet[i]);
	free(profileSet);
	for(i=0;i<profilenum;i++)
		free(tmp_profiles[i]);
	free(tmp_profiles);		
}
//read the begin to end line part of profile file and get their triples in proper part of the  triples vector
void getPartTriples(int genelen, int siglen, int profilenum, int linelen, int begin, int end,  char *file, struct Profile_triple * triples)
{
	int i;
	//allocate the temp memory
	short **profileSet = (short **)malloc(profilenum*sizeof(short *));
	for(i=0;i<profilenum;i++)
		profileSet[i] = (short *)malloc(genelen*sizeof(short));	
	
	//read file and get the proper data
	ReadFile(file, linelen, begin, end, profilenum, genelen, profileSet);
	
	//get the triple for every profile	
	for(i=begin;i<end;i++)
		triples[i] = getTriple(profileSet[i], genelen, siglen);
	
	//free the temp memory
	for(i=0;i<profilenum;i++)
		free(profileSet[i]);
	free(profileSet);		
}
*/

//compare two set
int cmpset(int *set1,int *set2,int n)
{
	int i;
	quiksortINT(set1,0,n-1);
	quiksortINT(set2,0,n-1);
	for(i=0;i<n;i++)
		if(set1[i]!=set2[i])
			return 0;
	return 1;
}

//is set2 in set1
int isInSet(int **set1,int *set2,int n,int iter)
{
	int i;
	for(i=iter-1;i>=0;i--)
		if(cmpset(set1[i],set2,n)==1)
			return 1;
	return 0;
}

//test Available RAM
/******
0:byte
1:KB
2:MB
3:GB
******/
unsigned long memoryAvailable(int base)
{
	struct sysinfo si;
    sysinfo(&si);
	unsigned long ram = si.freeram;
	switch (base){
		case 0:
			return ram;
		case 1:
			return ram/1024;
		case 2:
			return ram/1024/1024;
		case 3:
			return ram/1024/1024/1024;
	}
    //printf("Totalram:       %ld\n", si.totalram/1024/1024/1024);
    //printf("Available:      %ld\n", si.freeram/1024/1024/1024);
}



