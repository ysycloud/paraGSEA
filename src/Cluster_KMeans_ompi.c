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

#define MAXITER 3000
 
void split_data(int size, int n, int rank, int* begin, int* end, int* local_n); 
int cmpset(int *set1,int *set2,int n);
int isInSet(int **set1,int *set2,int n,int iter);

int main(int argc,char *argv[])
{	
	int i,j;
	int genelen;
	int local_profilenum,profilenum;
	int linelen;
	float **local_ES_Matrix;		//part of the ES_Matrix in this process
	int *local_classflag;
	int *global_classflag;
	float *loc_avesimilar;
	float *global_avesimilar;
	int	my_rank;   /* My process rank           */
    int	p;         /* The number of processes   */
    int source,dest;  
    int tag = 0;
    MPI_Status  status;
	int begin,end;
	int *cluster_center;
	int isbreak=0;
	int iternum=0;
	
	int leave,global_begin;
	int **cluster_centers_history;
		
	int corenum = atoi(argv[1]);
	int cluster_center_num = atoi(argv[2]);

	double start,finish,duration;
	
	/* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	MPI_Barrier(MPI_COMM_WORLD); 
	if(my_rank == 0){
		
		//max iteration correspond to a cluster centers set
		cluster_centers_history = (int **)malloc(MAXITER*sizeof(int*));
		for(i=0;i<MAXITER;i++)
			cluster_centers_history[i] = (int *)malloc(cluster_center_num*sizeof(int));
		
		printf("Matrix is Loading...!\n");
		GET_TIME(start);
	}
	

	
	char myfile[128];
	sprintf(myfile,"%s_%d.txt",argv[3],my_rank);
	
	//read file parameters in all processes
	ReadMatrixFilePara(myfile, &local_profilenum, &profilenum, &linelen);
	//printf("%d\t%d\t%d\n",local_profilenum,profilenum,linelen);
	
	//calculate the global begin line of profiles in each process
	leave = profilenum % p;
	if(my_rank < leave)
		global_begin = my_rank * local_profilenum;
	else
		global_begin = my_rank * local_profilenum + leave;
	
	//allocate memory for ES Matrix
	local_ES_Matrix = (float**)malloc(local_profilenum*sizeof(float *));
	for(i=0;i<local_profilenum;i++)
		local_ES_Matrix[i] = (float*)malloc(profilenum*sizeof(float));
		
	//read file and get the proper data
	ReadMatrixFile(myfile, linelen, 0, local_profilenum, local_profilenum, profilenum, local_ES_Matrix);
	
	//prework for memory malloc
	local_classflag = (int *)malloc(local_profilenum*sizeof(int));
	global_classflag = (int *)malloc(profilenum*sizeof(int));
	loc_avesimilar = (float *)malloc(local_profilenum*sizeof(float));
	global_avesimilar = (float *)malloc(profilenum*sizeof(float));

	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == 0){
		GET_TIME(finish);
		//compute the IO time
		duration = finish-start;     
		printf("loading IO and prework time : %.4f s\n",duration);

		printf("Paral KMeans compute the Cluster Centers is Starting...!\n");
		GET_TIME(start);
	}

	//generate the init cluster center
	cluster_center = (int *)malloc(cluster_center_num*sizeof(int));
	if(my_rank == 0){		
		float* tmp = (float *)malloc(cluster_center_num*sizeof(float));
		GetRandomSequence(profilenum, cluster_center_num ,tmp);
		//adjust the center num,because GetRandomSequence function get the seq from 1-total,but we need 0-total-1
		for(i=0;i<cluster_center_num;i++)
			cluster_center[i] = (int)tmp[i]-1;
		free(tmp);
	}
	MPI_Bcast(cluster_center, cluster_center_num, MPI_INT, 0 ,MPI_COMM_WORLD);
	
	if(my_rank == 0)
	{
		printf("Init cluster centers is:\n");
		for(i=0;i<cluster_center_num;i++)
			printf("%d ",cluster_center[i]);
		printf("\n");
	}
	
	while(!isbreak)
	{
		if(my_rank==0)
		{	
			memcpy(cluster_centers_history[iternum],cluster_center,cluster_center_num*sizeof(int));	
		}
		
		iternum++;
		/***************************paral split the class***********************************/
		#pragma omp parallel num_threads(corenum)
		{
			int k,t;
			int local_t;	
			int begin_t,end_t;
			int flag;
			int threadID = omp_get_thread_num();
			
			// compute the local size、up boundary and down boundary for every thread in dataset2
			split_data(local_profilenum, corenum, threadID, &begin_t, &end_t, &local_t);
			
			// compute the part of the local_classflag vector
			for(k=begin_t;k<end_t;k++)
			{
				flag = 0;		
				for(t=1;t<cluster_center_num;t++)
				{
					if((global_begin+k)==cluster_center[t]) //if the center is current profile ,break directly.
					{						
						flag = t;  
						break;
					}
					if(local_ES_Matrix[k][cluster_center[t]] > local_ES_Matrix[k][cluster_center[flag]])
						//biger ES，more similar，shorter distance
						flag = t;
				}									
				local_classflag[k] = cluster_center[flag];
			}				
		}
			
		/*****************Allgather local_classflag in all process**************/
		if(profilenum % p == 0){ //can split blanced
			MPI_Allgather(local_classflag,local_profilenum,MPI_INT,global_classflag,local_profilenum,MPI_INT,MPI_COMM_WORLD);
		}else{	//can not split blanced 
			if(my_rank == 0){
				/*****************gather local_classflag in process0**************/
				int leave = profilenum % p;
				//tmp memory for processes before leave
				int *tmp1 = (int *)malloc(local_profilenum*sizeof(int));
				//tmp memory for processes after leave
				int *tmp2 = (int *)malloc((local_profilenum-1)*sizeof(int));
				//copy local_classflag vector in process0 to global_classflag vector 
				memcpy(global_classflag,local_classflag,local_profilenum*sizeof(int));
				/********receive and copy local_classflag vector in other processes to global_classflag vector******/
				for(i=1; i<p; i++){
					if(i < leave){
						MPI_Recv(tmp1 , local_profilenum, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&global_classflag[i*local_profilenum],tmp1,local_profilenum*sizeof(int));
					}else{
						MPI_Recv(tmp2 , local_profilenum-1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&global_classflag[i*(local_profilenum-1)+leave],tmp2,(local_profilenum-1)*sizeof(int));
					}	
				}
				free(tmp1);
				free(tmp2);							
			}else{
				MPI_Send(local_classflag, local_profilenum, MPI_INT, 0, tag, MPI_COMM_WORLD);
			}
			//Bcast the global_classflag vector to all process
			MPI_Bcast( global_classflag, profilenum, MPI_INT, 0 ,MPI_COMM_WORLD);
		}
			
		/******************paral compute the average distance for every profile in his class**************************/
		#pragma omp parallel num_threads(corenum)
		{
			int k,t;
			int local_t;	
			int begin_t,end_t;
			float ave;
			int count;
			int myclass;
			int threadID = omp_get_thread_num();
			
			// compute the local size 、up boundary and down boundary for every thread in dataset2
			split_data(local_profilenum, corenum, threadID, &begin_t, &end_t, &local_t);
			
			// compute the average distance
			for(k=begin_t;k<end_t;k++)
			{
				myclass = local_classflag[k];
				ave = 0;
				count = 0;
				for(t=0;t<profilenum;t++)
					if(global_classflag[t] == myclass)
					{
						ave += local_ES_Matrix[k][t];
						count++;
					}
				loc_avesimilar[k] = ave/count;					
			}				
		}
		
		/*
		printf("rank_%d:\n",my_rank);
		for(i=0;i<local_profilenum;i++)
			printf("%f ",loc_avesimilar[i]);
		printf("\n");
		*/
		
		/*****************Allgather loc_avesimilar to process0**************/
		if(profilenum % p == 0){ //can split blanced
			MPI_Gather(loc_avesimilar,local_profilenum,MPI_FLOAT,global_avesimilar,local_profilenum,MPI_FLOAT,0, MPI_COMM_WORLD);
		}else{	//can not split blanced 
			if(my_rank == 0){
				/*****************gather loc_avesimilar in process0**************/
				int leave = profilenum % p;
				//tmp memory for processes before leave
				float *tmp1 = (float *)malloc(local_profilenum*sizeof(float));
				//tmp memory for processes after leave
				float *tmp2 = (float *)malloc((local_profilenum-1)*sizeof(float));
				//copy loc_avesimilar vector in process0 to global_avesimilar vector 
				memcpy(global_avesimilar,loc_avesimilar,local_profilenum*sizeof(float));
				/********receive and copy loc_avesimilar vector in other processes to global_avesimilar vector******/
				for(i=1; i<p; i++){
					if(i < leave){
						MPI_Recv(tmp1 , local_profilenum, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&global_avesimilar[i*local_profilenum],tmp1,local_profilenum*sizeof(float));
					}else{
						MPI_Recv(tmp2 , local_profilenum-1, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&global_avesimilar[i*(local_profilenum-1)+leave],tmp2,(local_profilenum-1)*sizeof(float));
					}	
				}
				free(tmp1);
				free(tmp2);							
			}else{
				MPI_Send(loc_avesimilar, local_profilenum, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			}
		}
			
		/*****************find the new cluster centers********************/
		if(my_rank==0)
		{	
			/*
			printf("average distance is:\n");
			for(i=0;i<profilenum;i++)
				printf("%d->%d--%f\t",i,global_classflag[i],global_avesimilar[i]);
			printf("\n");
			*/
						
			int *indexClass = (int *)malloc(profilenum*sizeof(int));//the num of class x
			float *maxsimilar = (float *)malloc(cluster_center_num*sizeof(float));// ith cluster's max average similarity
			int *cluster_center_new = (int *)malloc(cluster_center_num*sizeof(int));
			
			//get the index of every cluster center
			for(i=0;i<cluster_center_num;i++)
				indexClass[cluster_center[i]] = i;
			
			//init the tmp memory
			//memset(maxsimilar,-1,cluster_center_num*sizeof(float));
			for(i=0;i<cluster_center_num;i++)
				maxsimilar[i] = -2;
		//	memcpy(cluster_center_new,cluster_center,cluster_center_num*sizeof(int));
		
			int cluster_num;
			for(i=0;i<profilenum;i++)
			{
				cluster_num = indexClass[global_classflag[i]];
				if(global_avesimilar[i]>maxsimilar[cluster_num])
				{
					maxsimilar[cluster_num] = global_avesimilar[i];
					cluster_center_new[cluster_num] = i;
				}
			}
						
			printf("%dth iteraction cluster_center_new is:\n",iternum);
			for(i=0;i<cluster_center_num;i++)
				printf("%d ",cluster_center_new[i]);
			printf("\n");
			
			if(isInSet(cluster_centers_history,cluster_center_new,cluster_center_num,iternum)||iternum>=MAXITER){
				isbreak = 1;
			}
			memcpy(cluster_center,cluster_center_new,cluster_center_num*sizeof(int));	
		}
		MPI_Bcast( &isbreak, 1, MPI_INT, 0 ,MPI_COMM_WORLD);
		MPI_Bcast(cluster_center, cluster_center_num, MPI_INT, 0 ,MPI_COMM_WORLD);
	}
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == 0){
		WritetxtClusterResult(global_classflag ,profilenum, argv[4]);
		GET_TIME(finish);
		//compute the Write time
		duration = finish-start;     
		printf("Paral KMeans compute the Cluster Centers Spent: %.4f s\n",duration);
	}
	
	//free the memory
	if(my_rank == 0)
	{
		for(i=0;i<MAXITER;i++)
			free(cluster_centers_history[i]);
		free(cluster_centers_history);
	}	
	free(loc_avesimilar);	
	free(global_avesimilar);
	for(i=0;i<local_profilenum;i++)
		free(local_ES_Matrix[i]);
	free(local_ES_Matrix);
	free(local_classflag);
	free(global_classflag);
	free(cluster_center);
}

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