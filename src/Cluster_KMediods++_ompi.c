#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include "mpi.h"
#include <omp.h>
#include <unistd.h> 
#include <getopt.h> 
#include "Tools.h"

#define MAXHIST 10000

#define ERRM "Cluster error:"

char *USAGE =
"\n"
"Usage:"
"  Cluster_KMediods++_ompi [options]\n"
"\n"
"  general options before command by MPI:\n"
"	 -n process_num : Total number of processes. [ default 1 ]\n"
"	 -ppn pernum: the number of processes in each node. [ default 1 ]\n"
"	 -hostfile hostfile:  list the IP or Hostname of nodes. [ default localhost ]"
"\n"
"  general options:\n"
"    -t --thread: the number of threads in per process_num. [ default 1 ]\n"
"	 -c	--cluster: the number of clusters we want to get. [ default 5 ]\n"
"	 -w	--write: whether output the results . [ default 1]\n"				
"\n"
"  input/output options: \n"
"    -i --input1: distributed ES_Matrix file we get from stage 2(Compare Profiles)\n"
"	 -o --output: output class flags file of every profiles in root node\n"
"    -s --sample: input file/a parsed sample sequence number file from pretreatment stage. \n"
"    -r --reference: input a directory includes referenced files about genesymbols and cids. \n";

void Usage();

int main(int argc,char *argv[])
{	
	int i,j;
	int genelen;
	int local_profilenum,profilenum;
	int linelen;
	float **local_ES_Matrix;		//part of the ES_Matrix in this process
	int *local_classflag;
	int *global_classflag;
	float *local_leastcenter;
	float *global_leastcenter;
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
	int histnum=0;
	int ismaxhist=0;
	
	int leave,global_begin;
	int **cluster_centers_history;
	int parameternum;
	int corenum;
	int cluster_center_num;
	int ifwrite;

	double start,finish,duration;
	
	FILE *fp;
	
	/* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	/* check parameter */
	if(my_rank == 0)
	{
		parameternum = argc;
		if(parameternum==1)
			Usage();
	}
	MPI_Bcast(&parameternum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(parameternum==1)
	{
		MPI_Finalize();
		exit(0);
	}
	
	// Unset flags (value -1).
	corenum = -1;
	cluster_center_num = -1;
	ifwrite = -1;	
    // Unset options (value 'UNSET').
	char * const UNSET = "unset";
    char * input   = UNSET;
	char * output   = UNSET;
	char * sample   = UNSET;
	char * reference   = UNSET;
	
	int c;
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"thread",             required_argument,        0, 't'},
			{"cluster",             required_argument,       0, 'c'},
			{"write",              required_argument,        0, 'w'},
			{"input",             required_argument,         0, 'i'},
			{"sample",             required_argument,        0, 's'},
			{"reference",			required_argument,        0, 'r'},
			{"output",             required_argument,        0, 'o'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "t:c:w:i:s:r:o:",
            long_options, &option_index);
	
		if(c==-1)	break;
		
		switch (c) {
		
		case 0:
			// A flag was set. //
			break;

		case 'i':
			if (input == UNSET) 
			{
				input = optarg;
			}
			else 
			{
				if(my_rank==0)
				{
					fprintf(stderr, "%s --input set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();
				exit(0);
			}
			break;
			
		case 's':
			if (sample == UNSET) 
			{
				sample = optarg;
			}
			else 
			{
				if(my_rank==0)
				{
					fprintf(stderr, "%s --sample set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();
				exit(0);
			}
			break;
			
		case 'r':
			if (reference == UNSET) 
			{
				reference = optarg;
			}
			else 
			{
				if(my_rank==0)
				{
					fprintf(stderr, "%s --reference set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();
				exit(0);
			}
			break;
		
		case 'o':
			if (output == UNSET) 
			{
				output = optarg;
			}
			else 
			{
				if(my_rank==0)
				{
					fprintf(stderr, "%s --output set more than once\n", ERRM);
					Usage();
				}
				MPI_Finalize();
				exit(0);
			}
			break;
		
		case 't':
			if (corenum < 0) {
				corenum = atoi(optarg);
				if (corenum < 1) {
					if(my_rank==0)
					{
						fprintf(stderr, "%s --thread must be a positive integer\n", ERRM);
						Usage();
					}	
					MPI_Finalize();
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --thread set more " "than once\n", ERRM);
					Usage();
				}		
				MPI_Finalize();
				exit(0);
			}
			break;
			
		case 'c':
			if (cluster_center_num < 0) {
				cluster_center_num = atoi(optarg);
				if (cluster_center_num < 1) {
					if(my_rank==0)
					{
						fprintf(stderr, "%s --cluster must be a positive integer\n", ERRM);
						Usage();
					}	
					MPI_Finalize();
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --cluster set more " "than once\n", ERRM);
					Usage();
				}	
				MPI_Finalize();
				exit(0);
			}
			break;

		case 'w':
			if (ifwrite < 0) {
				ifwrite = atof(optarg);
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --write set more " "than once\n", ERRM);
					Usage();
				}		
				MPI_Finalize();
				exit(0);
			}
			break;
			
		default:
			// Cannot parse. //
			if(my_rank==0)
				Usage();
			MPI_Finalize();
			exit(0);
		}		
	}

	//check the parameters
	if(corenum == -1)
		corenum = 1;
	
	if(cluster_center_num == -1)
		cluster_center_num = 5;

	if(ifwrite == -1)
		ifwrite = 1;
	
	if((fp=fopen(sample,"r"))==NULL)
	{
		if(my_rank==0)
			fprintf(stderr, "can not open sample sequence number '%s' file\n",sample);
		MPI_Finalize();
		exit(0);
	}
	fclose(fp);
	
	char genelistfile[FileName_LEN];
	sprintf(genelistfile,"%s/Gene_List.txt",reference);
	
	if((fp=fopen(genelistfile,"r"))==NULL)
	{
		if(my_rank==0)
			fprintf(stderr, "the reference directory may be incorrect!\n");
		MPI_Finalize();
		exit(0);
	}
	fclose(fp);
	
	if(output == UNSET)
	{
		if(my_rank==0)
			fprintf(stderr,"Not Set output parameter!\n");
		MPI_Finalize();
		exit(0);
	}
	
	MPI_Barrier(MPI_COMM_WORLD); 
	if(my_rank == 0){
		
		//max iteration correspond to a cluster centers set
		cluster_centers_history = (int **)malloc(MAXHIST*sizeof(int*));
		for(i=0;i<MAXHIST;i++)
			cluster_centers_history[i] = (int *)malloc(cluster_center_num*sizeof(int));
		
		printf("Matrix is Loading...!\n");
		GET_TIME(start);
	}
	
	
	char myfile[128];
	sprintf(myfile,"%s_%d.txt",input,my_rank);
	
	//read file parameters in all processes
	ReadMatrixFilePara(myfile, &local_profilenum, &profilenum, &linelen);
	//printf("%d\t%d\t%d\n",local_profilenum,profilenum,linelen);
	
	//input file check
	if( local_profilenum <= 0 || profilenum <= 0 || linelen<=0)
	{
		if(my_rank==0)
			fprintf(stderr,"this file input is not exist!\n");
		MPI_Finalize();
		exit(0);
	}
	
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
	local_leastcenter = (float *)malloc(local_profilenum*sizeof(float));	
	global_leastcenter = (float *)malloc(profilenum*sizeof(float));
	loc_avesimilar = (float *)malloc(local_profilenum*sizeof(float));
	global_avesimilar = (float *)malloc(profilenum*sizeof(float));

	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == 0){
		GET_TIME(finish);
		//compute the IO time
		duration = finish-start;     
		printf("loading IO and prework time : %.4f s\n",duration);

		printf("Paral KMediods++ compute the Cluster Centers is Starting...!\n");
		GET_TIME(start);
	}

	/****************generate the init cluster center*****************************************/
	cluster_center = (int *)malloc(cluster_center_num*sizeof(int));
	int currentInitCenterNum = 0;
	int currentInitCenter;
	if(my_rank==0)
	{
		//generate the first center
		float tmp;
		GetRandomSequence(profilenum,1,&tmp);
		//adjust the center num,because GetRandomSequence function get the seq from 1-total,but we need 0-total-1
		currentInitCenter = (int)tmp-1;
	}
		
	for(; currentInitCenterNum<cluster_center_num; currentInitCenterNum++)
	{
		//bcast the current center to all processes
		MPI_Bcast(&currentInitCenter, 1 , MPI_INT, 0 ,MPI_COMM_WORLD);
		cluster_center[currentInitCenterNum] = currentInitCenter;
		
		if(currentInitCenterNum==(cluster_center_num-1))
			break;
		
		/******************paral split the each data to compute the least center***********************************/
		#pragma omp parallel num_threads(corenum)
		{
			int k,t;
			int local_t;	
			int begin_t,end_t;
			float least;
			int threadID = omp_get_thread_num();
			
			// compute the local size、up boundary and down boundary for every thread in dataset2
			split_data(local_profilenum, corenum, threadID, &begin_t, &end_t, &local_t);
			
			// compute the part of the least center vector
			for(k=begin_t;k<end_t;k++)
			{
				least = -1;	 //init  a most far distance	
				for(t=0;t<=currentInitCenterNum;t++)
				{
					if((global_begin+k)==cluster_center[t])
					//if the center is current profile ,break and set the least dis is 1,
					//stop it be choosed as center again						
					{
						least = 1;
						break;
					}						
					if(local_ES_Matrix[k][cluster_center[t]] > least)
						//biger ES，more similar，shorter distance
						least = local_ES_Matrix[k][cluster_center[t]];
				}									
				local_leastcenter[k] = least;
			}				
		}
		
		/*****************gather local_leastcenter to root process**************/
		if(profilenum % p == 0){ //can split blanced
			MPI_Gather(local_leastcenter,local_profilenum,MPI_FLOAT,global_leastcenter,local_profilenum,MPI_FLOAT,0,MPI_COMM_WORLD);
		}else{	//can not split blanced 
			if(my_rank == 0){
				/*****************gather local_leastcenter in process0**************/
				int leave = profilenum % p;
				//tmp memory for processes before leave
				float *tmp1 = (float *)malloc(local_profilenum*sizeof(float));
				//tmp memory for processes after leave
				float *tmp2 = (float *)malloc((local_profilenum-1)*sizeof(float));
				//copy local_leastcenter vector in process0 to global_leastcenter vector 
				memcpy(global_leastcenter,local_leastcenter,local_profilenum*sizeof(float));
				/********receive and copy local_leastcenter vector in other processes to global_leastcenter vector******/
				for(i=1; i<p; i++){
					if(i < leave){
						MPI_Recv(tmp1 , local_profilenum, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&global_leastcenter[i*local_profilenum],tmp1,local_profilenum*sizeof(float));
					}else{
						MPI_Recv(tmp2 , local_profilenum-1, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&global_leastcenter[i*(local_profilenum-1)+leave],tmp2,(local_profilenum-1)*sizeof(float));
					}	
				}
				free(tmp1);
				free(tmp2);							
			}else{
				MPI_Send(local_leastcenter, local_profilenum, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			}
		}
		
		/*****************generate the current init center**************/
		if(my_rank==0)
		{
			float far = 1;
			currentInitCenter = 0;
			//find the farest least center
			for(i=0;i<profilenum;i++)
			{
				if(global_leastcenter[i]<far)
				{
					far = global_leastcenter[i];
					currentInitCenter = i;
				}
			}			
		}
	}
		
	free(local_leastcenter);
	free(global_leastcenter);
	/**************************generate the init cluster center END************************************/
	MPI_Barrier(MPI_COMM_WORLD);
	
	
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
			memcpy(cluster_centers_history[histnum],cluster_center,cluster_center_num*sizeof(int));				
		}
		histnum++;
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
			
			if(ismaxhist == 1){
				if(isInSet(cluster_centers_history,cluster_center_new,cluster_center_num,MAXHIST)){
					isbreak = 1;
				}
			}else{
				if(isInSet(cluster_centers_history,cluster_center_new,cluster_center_num,histnum)){
					isbreak = 1;
				}
			}
			
			if(histnum >= MAXHIST)
			{
				histnum = 0;
				ismaxhist = 1;
			}
				
			memcpy(cluster_center,cluster_center_new,cluster_center_num*sizeof(int));	
		}
		MPI_Bcast( &isbreak, 1, MPI_INT, 0 ,MPI_COMM_WORLD);
		MPI_Bcast(cluster_center, cluster_center_num, MPI_INT, 0 ,MPI_COMM_WORLD);
	}
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == 0){
		
		if(ifwrite==1)
			WritetxtClusterResult(global_classflag ,profilenum,cluster_center_num , output, sample, reference);
		else
			printf("Just run for test, no results output\n");

		GET_TIME(finish);
		//compute the Write time
		duration = finish-start;     
		printf("Paral KMediods++ compute the Cluster Centers Spent: %.4f s\n",duration);
	}
	
	//free the memory
	if(my_rank == 0)
	{
		for(i=0;i<MAXHIST;i++)
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
	
	MPI_Finalize();
	return 0;
}

void Usage() {
	fprintf(stderr, "%s\n", USAGE);
}  /* Usage */