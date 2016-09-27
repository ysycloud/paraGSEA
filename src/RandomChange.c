#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include "RandomChange.h"

//generate 'size' not repeat sequence number randomly from 1-'total'
void GetRandomSequence(int total,int size , float* Seq)
{
	int i;
	int *sequence = (int *)malloc(total*sizeof(int));

    for (i = 0; i < total; i++)
	{
        sequence[i] = i+1;
    }
	
	struct timeval seed;
	gettimeofday(&seed, NULL);
	srand((unsigned)seed.tv_usec);  //must set microseconds as unit
    int end = total-1;

    for (i = 0; i < size; i++) 
	{
        int num = rand()%(end+1);
        Seq[i] = (float)sequence[num];
        sequence[num] = sequence[end];
        end--;
    }
	free(sequence);
}

//permutate the positions of array randomly
void changePosition(float *positions , int size) 
{   
	struct timeval seed;
	gettimeofday(&seed, NULL);
	srand((unsigned)seed.tv_usec);
	
	int index,rd;
	float temp;
    for( index=size-1; index>=size/2; index--) {
		rd = rand()%(index+1);
		temp = positions[rd];
		positions[rd] = positions[index];   
        positions[index] = temp;   
	}     
}   