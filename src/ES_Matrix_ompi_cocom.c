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


#define ERRM "ES Matrix error:"

char *USAGE =
"\n"
"Usage:"
"  ES_Matrix_ompi_cocom [options]\n"
"\n"
"  general options before command by MPI:\n"
"	 -n process_num : Total number of processes. [ default 1 ]\n"
"	 -ppn pernum: the number of processes in each node. [ default 1 ]\n"
"	 -hostfile hostfile:  list the IP or Hostname of nodes. [ default localhost ]"
"\n"
"  general options:\n"
"    -t --thread: the number of threads in per process_num. [ default 1 ]\n"
"	 -l	--siglen: the length of Gene Expression Signature. [ default 50 ]\n"
"	 -a	--loadtime: the load time of dataset2. [ default 1 (>=1)]\n"
"	 -p	--proportion: the proportion of dataset be used . [ default 1 (0,1]]\n"
"	 -w	--write: whether output the results . [ default 1]\n"
"\n"
"  input/output options: \n"
"    -1 --input1: a parsed profiles's file from pretreatment stage.\n"
"    -2 --input2: another parsed profiles's file from pretreatment stage.\n"
"	 -o --output: output file ,distributed in every nodes ,with ES Matrix\n";

void Usage();
void Build_derived_type(
		struct Profile_triple* m_ptr, 			 /*  in  */
		MPI_Datatype* triple_mpi_t_ptr 		 /*  out  */);

int main(int argc,char *argv[])
{	
	int i,j;
	int genelen;
	int profilenum1,profilenum2;
	int linelen1,linelen2;
	struct Profile_triple *triples1,**triples2,**local_triples2;
	float **local_ES_Matrix;		//part of the ES_Matrix in this process
	float *ES_test;  //ES used for testing without writing
	int	my_rank;   /* My process rank           */
    int	p;         /* The number of processes   */
    int source,dest;  
    int tag = 0;
    MPI_Status  status;
	MPI_Datatype triple_mpi_t;
	int local_P;	//the data number of each processes must hand
	int begin,end;
	int begin_d2,end_d2,local_P_d2;
	int corenum;
	int siglen;
	int load_time;
	int parameternum;
	float proportion;
	int ifwrite;
	
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
		if(parameternum == 1)
			Usage();
	}
	MPI_Bcast(&parameternum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(parameternum == 1)
	{
		MPI_Finalize();
		exit(0);
	}
	
	// Unset flags (value -1).
	corenum = -1;
	siglen = -1;
	load_time = -1;	
	proportion = -1;
	ifwrite = -1;
    // Unset options (value 'UNSET').
	char * const UNSET = "unset";
    char * input1   = UNSET;
	char * input2   = UNSET;
	char * output   = UNSET;
	
	int c;
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"thread",             required_argument,        0, 't'},
			{"siglen",             required_argument,        0, 'l'},
			{"loadtime",           required_argument,        0, 'a'},
			{"proportion",         required_argument,        0, 'p'},
			{"write",              required_argument,        0, 'w'},
			{"input1",             required_argument,        0, '1'},
			{"input2",             required_argument,        0, '2'},
			{"output",             required_argument,        0, 'o'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "t:l:a:p:w:1:2:o:",
            long_options, &option_index);
	
		if(c==-1)	break;
		
		switch (c) {
		
		case 0:
			// A flag was set. //
			break;

		case '1':
			if (input1 == UNSET) 
			{
				input1 = optarg;
			}
			else 
			{
				if(my_rank==0)
				{
					fprintf(stderr, "%s --input1 set more than once\n", ERRM);
					Usage();
				}	
				MPI_Finalize();
				exit(0);
			}
			break;
		
		case '2':
			if (input2 == UNSET) 
			{
				input2 = optarg;
			}
			else 
			{
				if(my_rank==0)
				{
					fprintf(stderr, "%s --input2 set more than once\n", ERRM);
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
			
		case 'l':
			if (siglen < 0) {
				siglen = atoi(optarg);
				if (siglen < 1) {
					if(my_rank==0)
					{
						fprintf(stderr, "%s --siglen must be a positive integer\n", ERRM);
						Usage();
					}		
					MPI_Finalize();
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --siglen set more " "than once\n", ERRM);
					Usage();
				}	
				MPI_Finalize();
				exit(0);
			}
			break;
		
		case 'a':
			if (load_time < 0) {
				load_time = atoi(optarg);
				if (load_time < 1) {
					if(my_rank==0)
					{
						fprintf(stderr, "%s --load time must be a positive integer\n", ERRM);
						Usage();
					}		
					MPI_Finalize();
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --load time set more " "than once\n", ERRM);
					Usage();
				}		
				MPI_Finalize();
				exit(0);
			}
			break;

		case 'p':
			if (proportion < 0) {
				proportion = atof(optarg);
				if (proportion > 1 || proportion <= 0) {
					if(my_rank==0)
					{
						fprintf(stderr, "%s -- proportion must be kept in (0,1]\n", ERRM);
						Usage();
					}		
					MPI_Finalize();
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --proportion set more " "than once\n", ERRM);
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
	
	if(siglen == -1)
		siglen = 50;
	
	if(load_time == -1)
		load_time = 1;

	if(proportion == -1)
		proportion = 1;

	if(ifwrite == -1)
		ifwrite = 1;
	
	if(output == UNSET)
	{
		if(my_rank==0)
			fprintf(stderr," [ param error : -o ] Not Set output parameter!\n");
		MPI_Finalize();
		exit(0);
	}
	
	triples2 = (struct Profile_triple **)malloc(load_time*sizeof(struct Profile_triple *));
	local_triples2 = (struct Profile_triple **)malloc(load_time*sizeof(struct Profile_triple *));
	
	//barrier all processes to compute time
	MPI_Barrier(MPI_COMM_WORLD); 
	if(my_rank == 0){
		printf("Profile Set is Loading...!\n");
		GET_TIME(start);
	}
	
	//read file parameters in all processes
	ReadFilePara(input1, &profilenum1, &genelen, &linelen1);
	ReadFilePara(input2, &profilenum2, &genelen, &linelen2);
	
	profilenum1 *= proportion;
	profilenum2 *= proportion;
	
	//input file check
	if( profilenum1 <= 0 || genelen <= 0)
	{
		if(my_rank==0)
			fprintf(stderr," [ param error : -1 ] this file input1 is not exist!\n");
		MPI_Finalize();
		exit(0);
	}
	
	if( profilenum2 <= 0 || genelen <= 0)
	{
		if(my_rank==0)
			fprintf(stderr," [ param error : -2 ] this file input2 is not exist!\n");
		MPI_Finalize();
		exit(0);
	}

	if(my_rank==0)
	{
		printf("Genelen:	%d\n", genelen);
		printf("Profiles1 length:	%d\n", profilenum1);
		printf("Profiles2 length:	%d\n", profilenum2);
	}
	
	// compute the local size 、up boundary and down boundary for every process in dataset1
	split_data(profilenum1, p, my_rank, &begin, &end, &local_P);
	
	if(my_rank==0)
		printf("Memory check......\n");
	
	unsigned long memavail = memoryAvailable(1);
	
	unsigned long memneed = sizeof(struct Profile_triple)/1024*(local_P+profilenum2/load_time) + local_P/1024*profilenum2*sizeof(float);
	
	unsigned long memallneed = sizeof(struct Profile_triple)/1024*(profilenum1+profilenum2) + profilenum1/1024*profilenum2*sizeof(float);
	
	if(my_rank==0)
	{
		printf("Available Memory:      %ld KB\n", memavail);
		printf("Needed Memory:      %ld KB\n", memneed);
		printf("All Needed Memory:      %ld KB\n", memallneed);
	}

	unsigned long mem1 = sizeof(struct Profile_triple)/1024*(local_P+profilenum2/load_time);
	unsigned long mem2 = profilenum1/1024*profilenum2*sizeof(float);
	
	int nodenum = (int)(mem2/(memavail-mem1)+1);
	if( memneed > memavail )
	{
		if( my_rank==0 )
		{
			//printf("mem1:      %ld KB\n", mem1);
			//printf("mem2:      %ld KB\n", mem2);
			printf("available memory is not enough to store all results, recommend to use more than %d nodes!!!\n", nodenum);
		}
		if(ifwrite==1){
			MPI_Finalize();
			exit(0);
		}else{
			if( my_rank==0 )
			{
				printf("because we are just testing without writing, we will continue!!!\n");
			}
		}
	}
	
	/*****read the local part file of dataset1 in every process and get their triples****************/
	triples1 = (struct Profile_triple *)malloc(sizeof(struct Profile_triple)*local_P);	
	getTriples(local_P, genelen, siglen, profilenum1, linelen1, begin, end, input1, triples1);
	Build_derived_type(&triples1[0],&triple_mpi_t);  //derive the new MPI Type
	
	
	if(ifwrite==1)
	{
		//not test, then allocate the local_ES_Matrix memory
		local_ES_Matrix = (float **)malloc(local_P*sizeof(float *));
		for(i=0;i<local_P;i++)
			local_ES_Matrix[i] = (float *)malloc(profilenum2*sizeof(float));
	}else{
		ES_test = (float *)malloc(corenum*sizeof(float));
	}
	
	int current_time = 0;
	int begin_localfile2, end_localfile2, len_localfile2;
	
	while( current_time < load_time )
	{
		/********************para load profile dataset2 by openmp******************************/
		split_data(profilenum2, load_time, current_time, &begin_localfile2, &end_localfile2, &len_localfile2);
		//allocate the triples memory for dataset2
		triples2[current_time] = (struct Profile_triple *)malloc(sizeof(struct Profile_triple)*len_localfile2);
	
		// compute the local size 、up boundary and down boundary for every process in dataset2
		split_data(len_localfile2, p, my_rank, &begin_d2, &end_d2, &local_P_d2);
	
		//allocate the triples memory for tmp local dataset2
		local_triples2[current_time] = (struct Profile_triple *)malloc(sizeof(struct Profile_triple)*local_P_d2);
		
		/********************every process para load profile each-time dataset2 by openmp******************************/		
		#pragma omp parallel num_threads(corenum)
		{
			int local_t;	//the data number of each thread must hand
			int begin_t,end_t;
			int threadID = omp_get_thread_num();
		
			// compute the local size 、up boundary and down boundary for every thread in dataset2
			split_data(local_P_d2, corenum, threadID, &begin_t, &end_t, &local_t);
		
			// compute the begin_t to end_t triples
			getFreeTriples(genelen, siglen, profilenum2, linelen2, begin_localfile2 + begin_d2 + begin_t, begin_t, local_t, input2, local_triples2[current_time]);		
		}
	
	
		/*****************Allgather part of data2's triple in all process**************/
		if(len_localfile2 % p == 0){ //can split blanced
			MPI_Allgather(local_triples2[current_time],local_P_d2,triple_mpi_t,triples2[current_time],local_P_d2,triple_mpi_t,MPI_COMM_WORLD);
			free(local_triples2[current_time]);
		}else{	//can not split blanced 
			if(my_rank == 0){
				/*****************gather data2's triple in process0**************/
				int leave = len_localfile2 % p;
				//tmp memory for processes before leave
				struct Profile_triple *tmp1 = (struct Profile_triple*)malloc(local_P_d2*sizeof(struct Profile_triple));
				//tmp memory for processes after leave
				struct Profile_triple *tmp2 = (struct Profile_triple*)malloc((local_P_d2-1)*sizeof(struct Profile_triple));
				//copy local triples vector in process0 to global triples vector 
				memcpy(triples2[current_time],local_triples2[current_time],local_P_d2*sizeof(struct Profile_triple));
				/********receive and copy data2's local triple vector in other processes to global triple vector******/
				for(i=1; i<p; i++){
					if(i < leave){
						MPI_Recv(tmp1 , local_P_d2, triple_mpi_t, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&triples2[current_time][i*local_P_d2],tmp1,local_P_d2*sizeof(struct Profile_triple));
					}else{
						MPI_Recv(tmp2 , local_P_d2-1, triple_mpi_t, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&triples2[current_time][i*(local_P_d2-1)+leave],tmp2,(local_P_d2-1)*sizeof(struct Profile_triple));
					}	
				}
				free(tmp1);
				free(tmp2);							
			}else{
				MPI_Send(local_triples2[current_time], local_P_d2, triple_mpi_t, 0, tag, MPI_COMM_WORLD);
			}
			free(local_triples2[current_time]);
			//Bcast the dataset2‘s triples to all process
			MPI_Bcast( triples2[current_time], len_localfile2, triple_mpi_t, 0 ,MPI_COMM_WORLD);
		}
	 	
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0){
			GET_TIME(finish);
			//compute the IO time
			duration = finish-start;     
			printf("phase %d -->  loading IO and prework time in collective comunication way: %.4f s\n", current_time+1, duration);

			printf("phase %d -->  Paral compute the ES_Matrix is Starting...!\n", current_time+1);
			GET_TIME(start);
		}
	
		/*
		if(my_rank == 0){
			int k;
			for(k=0;k<siglen;k++)
				printf("%d ",triples1[0].gsUp[k]);
			printf("\n");
			for(k=0;k<genelen;k++)
				printf("%d ",triples2[current_time][len_localfile2-1].index[k]);
			printf("\n");
		}
		*/
	
		/********************para compute the part of ES_Matrix******************************/
		#pragma omp parallel num_threads(corenum)
		{
			int k,t;
			int local_t;	//the data number of each thread must hand
			int begin_t,end_t;
			int threadID = omp_get_thread_num();
		
			// compute the local size 、up boundary and down boundary for every thread in dataset2
			split_data(len_localfile2, corenum, threadID, &begin_t, &end_t, &local_t);
		
			// compute the part of the ES matrix
			if(ifwrite==1){
				for(k=0;k<local_P;k++)
					for(t=begin_t;t<end_t;t++)
						local_ES_Matrix[k][begin_localfile2 + t] = ES_Profile_triple(triples1[k],triples2[current_time][t],genelen,siglen);
			}else{ //just calculate for testing
				for(k=0;k<local_P;k++)
					for(t=begin_t;t<end_t;t++)
						ES_test[threadID] = ES_Profile_triple(triples1[k],triples2[current_time][t],genelen,siglen);
			}
		}
	
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0){
			GET_TIME(finish);
			//compute the compute time
			duration = finish-start;     
			printf("phase %d --> Paral compute the ES_Matrix time: %.4f s\n", current_time+1, duration);
		
			if(current_time==load_time-1)
				if(ifwrite == 1)
					printf("Writing file is Starting...!\n");
			
			GET_TIME(start);		
		}
			
		free(triples2[current_time]);
		current_time++;	
	}
	
/*
	if(my_rank == 0){
		int k;
		for(k=0;k<profilenum2;k++)
			printf("%f ",local_ES_Matrix[0][k]);
		printf("\n");
	}
*/	
	
	if(ifwrite==1)
	{
		char Res[128];
		sprintf(Res,"%s_%d.txt",output,my_rank);
		WritetxtResult(0, local_P, profilenum2, Res, local_ES_Matrix);
	
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0){
			GET_TIME(finish);
			//compute the write time
			duration = finish-start;     
			printf("Write Result spent: %.4f s\n",duration);
		}
	}else{
		if(my_rank == 0){   
			printf("Just run for test, no results output\n");
		}
	}
	
	//free the memory
	if(ifwrite==1)
	{
		for(i=0;i<local_P;i++)
			free(local_ES_Matrix[i]);
		free(local_ES_Matrix);
	}else{
		free(ES_test);
	}
	free(triples1);
	free(triples2);
	
	MPI_Finalize();
	return 0;
}

void Build_derived_type(
		struct Profile_triple* m_ptr, 	 /*  in  */
		MPI_Datatype* triple_mpi_t_ptr /*  out  */)
{
		
	int block_lengths[4];
	MPI_Aint displacements[4];
	MPI_Datatype typelist[4];
		
	MPI_Aint start_address;
	MPI_Aint address;
		
	block_lengths[0] = MAX_GENESET;
	block_lengths[1] = MAX_GENESET;
	block_lengths[2] = MAX_GENE;
	block_lengths[3] = 1;
						
	typelist[0] = MPI_SHORT;
	typelist[1] = MPI_SHORT;
	typelist[2] = MPI_SHORT;
	typelist[3] = MPI_INT;
		
	displacements[0] = 0;
		
	MPI_Address(m_ptr, &start_address);
	MPI_Address(m_ptr->gsDown, &address);
	displacements[1] = address - start_address;
	MPI_Address(m_ptr->index, &address);
	displacements[2] = address - start_address;
	MPI_Address(&m_ptr->cid, &address);
	displacements[3] = address - start_address;
		
	MPI_Type_struct(4,block_lengths,displacements,typelist,triple_mpi_t_ptr);
	MPI_Type_commit(triple_mpi_t_ptr);
}

void Usage() {
	fprintf(stderr, "%s\n", USAGE);
}  /* Usage */