#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include "mpi.h"
#include <omp.h>
#include "RandomChange.h"
#include "GSEA.h"
#include "IO.h"
 
void Usage(char prog_name[]);
void split_data(int size, int n, int rank, int* begin, int* end, int* local_n); 
void getTriples(int local_P, int genelen, int siglen, int profilenum, int linelen, int begin, int end,  char *file, struct Profile_triple * triples);
void getPartTriples(int genelen, int siglen, int profilenum, int linelen, int begin, int end,  char *file, struct Profile_triple * triples);

int main(int argc,char *argv[])
{	
	int i,j;
	int genelen;
	int profilenum1,profilenum2;
	int linelen1,linelen2;
	struct Profile_triple *triples1,*triples2;
	float **local_ES_Matrix;		//part of the ES_Matrix in this process
	int	my_rank;   /* My process rank           */
    int	p;         /* The number of processes   */
    int source,dest;  
    int tag = 0;
    MPI_Status  status;
	int local_P;	//the data number of each processes must hand
	int begin,end;
	int parameternum;
	int corenum;
	int siglen;

	double start,finish,duration;
	
	/* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	/* check parameter*/
	if(my_rank == 0)
	{
		parameternum = argc;
		if(parameternum!=6)
			Usage(argv[0]);
	}
	MPI_Bcast(&parameternum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(parameternum!=6)
	{
		MPI_Finalize();
		exit(0);
	}
	
	corenum = atoi(argv[1]);
	siglen = atoi(argv[2]);
	//barrier all processes to compute time
	MPI_Barrier(MPI_COMM_WORLD); 
	if(my_rank == 0){
		printf("Profile Set is Loading...!\n");
		GET_TIME(start);
	}
	
	//read file parameters in all processes
	ReadFilePara(argv[3], &profilenum1, &genelen, &linelen1);
	ReadFilePara(argv[4], &profilenum2, &genelen, &linelen2);

	// compute the local size 、up boundary and down boundary for every process in dataset1
	split_data(profilenum1, p, my_rank, &begin, &end, &local_P);
	
	/*****read the local part file of dataset1 in every process and get their triples****************/
	triples1 = (struct Profile_triple *)malloc(sizeof(struct Profile_triple)*local_P);	
	getTriples(local_P, genelen, siglen, profilenum1, linelen1, begin, end, argv[3], triples1);
	
	/********************para load profile dataset2 by openmp******************************/
	//allocate the triples memory for dataset2
	triples2 = (struct Profile_triple *)malloc(sizeof(struct Profile_triple)*profilenum2);
	#pragma omp parallel num_threads(corenum)
	{
		int local_t;	//the data number of each thread must hand
		int begin_t,end_t;
		int threadID = omp_get_thread_num();
		
		// compute the local size 、up boundary and down boundary for every thread in dataset2
		split_data(profilenum2, corenum, threadID, &begin_t, &end_t, &local_t);
		
		// compute the begin_t to end_t triples
		getPartTriples(genelen, siglen, profilenum2, linelen2, begin_t, end_t, argv[4], triples2);		
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == 0){
		GET_TIME(finish);
		//compute the IO time
		duration = finish-start;     
		printf("loading IO and prework time : %.4f s\n",duration);

		printf("Paral compute the ES_Matrix is Starting...!\n");
		GET_TIME(start);
	}
	
	/*
	if(my_rank == 0){
		int k;
		for(k=0;k<siglen;k++)
			printf("%d ",triples1[0].gsUp[k]);
		printf("\n");
		for(k=0;k<genelen;k++)
			printf("%d ",triples2[profilenum2-1].index[k]);
		printf("\n");
	}
	*/
	
	/********************para compute the part of ES_Matrix******************************/
	//allocate the local_ES_Matrix memory
	local_ES_Matrix = (float **)malloc(local_P*sizeof(float *));
	for(i=0;i<local_P;i++)
		local_ES_Matrix[i] = (float *)malloc(profilenum2*sizeof(float));
	#pragma omp parallel num_threads(corenum)
	{
		int k,t;
		int local_t;	//the data number of each thread must hand
		int begin_t,end_t;
		int threadID = omp_get_thread_num();
		
		// compute the local size 、up boundary and down boundary for every thread in dataset2
		split_data(profilenum2, corenum, threadID, &begin_t, &end_t, &local_t);
		
		// compute the part of the ES matrix
		for(k=0;k<local_P;k++)
			for(t=begin_t;t<end_t;t++)
				local_ES_Matrix[k][t] = ES_Profile_triple(triples1[k],triples2[t],genelen,siglen);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == 0){
		GET_TIME(finish);
		//compute the compute time
		duration = finish-start;     
		printf("Paral compute the ES_Matrix time : %.4f s\n",duration);
		
		printf("Writing file is Starting...!\n");
		GET_TIME(start);
	}
	
	/*
	if(my_rank == 0){
		int k;
		for(k=0;k<profilenum2;k++)
			printf("%f ",local_ES_Matrix[0][k]);
		printf("\n");
	}
	*/
	
	char Res[128];
	sprintf(Res,"%s_%d.txt",argv[5],my_rank);
	WritetxtResult(0, local_P, profilenum2, Res, local_ES_Matrix);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == 0){
		GET_TIME(finish);
		//compute the write time
		duration = finish-start;     
		printf("Write Result spent: %.4f s\n",duration);
	}
	
	//free the memory
	for(i=0;i<local_P;i++)
		free(local_ES_Matrix[i]);
	free(local_ES_Matrix);
	free(triples1);
	free(triples2);
	
	MPI_Finalize();
	return 0;
}

void Usage(char prog_name[]) {
	fprintf(stderr, "usage: <total_process_num> <per_num_in_each_process> <hostfile> %s\n",prog_name);
	fprintf(stderr, " <thread_num>  <Expression Signature Length>\n");
	fprintf(stderr, " <inputfile1>  <inputfile2>\n");
	fprintf(stderr, " <outputfile(ES_Matrix)>\n");
}  /* Usage */

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
