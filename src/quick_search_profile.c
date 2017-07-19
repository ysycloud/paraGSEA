#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include <unistd.h> 
#include <getopt.h> 
#include "mpi.h"
#include <omp.h>
#include "Tools.h"

#define ERRM "quick search profile error:"

char *USAGE =
"\n"
"Usage:"
"  quick_search_profile [options]\n"
"\n"
"  general options before command by MPI:\n"
"	 -n process_num : Total number of processes. [ default 1 ]\n"
"	 -ppn pernum: the number of processes in each node. [ default 1 ]\n"
"	 -hostfile hostfile:  list the IP or Hostname of nodes. [ default localhost ]"
"\n"
"  general options:\n"
"    -n --topn: The first and last N GSEA records ordered by ES. [ default 10 ]\n"
"    -t --thread: the number of threads in per process_num. [ default 1 ]\n"
"	 -l	--siglen: the length of Gene Expression Signature. [ default 50 ]\n"
"\n"
"  input/output options: \n"
"    -i --input: input file/a parsed profiles's file from pretreatment stage. \n"
"    -s --sample: input file/a parsed sample sequence number file from pretreatment stage. \n"
"    -r --reference: input a directory includes referenced files about genesymbols and cids. \n";

void Usage();
void Build_derived_type(
		struct ES_RESULT* m_ptr, 			 /*  in  */
		MPI_Datatype* es_mpi_t_ptr 		 /*  out  */);
 
int main(int argc,char *argv[])
{	
	int i,profilenum,genelen,linelen;
	struct Profile_triple *profileSet;
	struct Profile_triple searchProfile;
	char gsStr[1024];
	
	struct ES_RESULT *es_result;
	struct ES_RESULT *local_es;
	
	short target_profile[L1000_LEN];
	int target_profile_len;
	struct original_Profile original_target_profile[L1000_LEN];	
	
	MPI_Datatype es_mpi_t;
	int	my_rank;   /* My process rank           */
    int	p;         /* The number of processes   */
    int source,dest ;  
    int tag = 0;
    MPI_Status  status;
	int local_n;	//the data number of each processes must hand
	int parameternum;
	double start,finish,duration;
	
	
	FILE *fp;
	char conditionsfile[FileName_LEN];
	char offsetfile[FileName_LEN];
	char genelistfile[FileName_LEN];
	int input_way;
	
	char conditions[L1000_CONDITION_LEN];
	long cidnum;
	long offset;
	
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
	int TopN = -1;
	int corenum = -1;
	int siglen = -1;
    // Unset options (value 'UNSET').
	char * const UNSET = "unset";
    char * input   = UNSET;
	char * sample   = UNSET;
	char * reference   = UNSET;
	
	int c;
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"topn",              required_argument,        0, 'n'},
			{"thread",             required_argument,        0, 't'},
			{"siglen",             required_argument,        0, 'l'},
			{"input",             required_argument,        0, 'i'},
			{"sample",             required_argument,        0, 's'},
			{"reference",             required_argument,        0, 'r'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "n:t:l:i:s:r:",
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
		
		case 'n':
			if (TopN < 0) {
				TopN = atoi(optarg);
				if (TopN < 1) {
					if(my_rank==0)
					{
						fprintf(stderr, "%s --topn must be a positive integer\n", ERRM);
						Usage();
					}			
					MPI_Finalize();					
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --topn set more " "than once\n", ERRM);
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
			
		default:
			// Cannot parse. //
			if(my_rank==0)
				Usage();
			MPI_Finalize();
			exit(0);
		}		
	}

	//check the parameters
	if(TopN==-1)
		TopN = 10;
	
	if(corenum == -1)
		corenum = 1;
	
	if(siglen == -1)
		siglen = 50;
	
	if((fp=fopen(sample,"r"))==NULL)
	{
		if(my_rank==0)
			fprintf(stderr, " [ param error : -s ] can not open sample sequence number '%s' file\n",sample);
		MPI_Finalize();
		exit(0);
	}
	fclose(fp);
	
	sprintf(genelistfile,"%s/Gene_List.txt",reference);
	
	if((fp=fopen(genelistfile,"r"))==NULL)
	{
		if(my_rank==0)
			fprintf(stderr, " [ param error : -r ] the reference directory may be incorrect!\n");
		MPI_Finalize();
		exit(0);
	}
	fclose(fp);
	
	sprintf(conditionsfile,"%s/Samples_Condition.txt",reference);
	sprintf(offsetfile,"%s/Samples_RowByteOffset.txt",reference);

	//barrier all processes to compute time
	MPI_Barrier(MPI_COMM_WORLD); 
	if(my_rank == 0){
		printf("Profile Set is Loading...!\n");
		GET_TIME(start);
	}
	
	//read file parameters in all processes
	ReadFilePara(input, &profilenum, &genelen, &linelen);	
	
	if( profilenum <= 0 || genelen <= 0)
	{
		if(my_rank==0)
			fprintf(stderr," [ param error : -i ] this file is not exist!\n");
		MPI_Finalize();
		exit(0);
	}
	
	if(my_rank == 0){
		printf("profilenum:%d\t genelen:%d\n",profilenum,genelen);
		
		printf("Memory check......\n");	
	}
	
	unsigned long memavail = memoryAvailable(1);
	unsigned long memneed = (2*sizeof(short)*(profilenum/p+1)*genelen + profilenum*sizeof(struct ES_RESULT))/1024;
	
	if(my_rank == 0 )
	{
		printf("Available Memory:      %ld KB\n", memavail);		
		printf("Needed Memory:      %ld KB\n", memneed);
		
		if(memavail < memneed)
			printf("available memory is not enough!!! Please use more nodes!!!\n");
		else
			//malloc global GSEA para Vector for process0
			es_result = (struct ES_RESULT*)malloc(profilenum*sizeof(struct ES_RESULT));
	}
	
	if(memavail < memneed)	
		return;	
			
	// compute the local size ¡¢up boundary and down boundary for every process
	int begin,end;	
	split_data(profilenum, p, my_rank, &begin, &end, &local_n);
	int leave = profilenum % p;
	
	local_es = (struct ES_RESULT*)malloc(local_n*sizeof(struct ES_RESULT));
	
	/********************para load profile dataset******************************/
	//para read the file to local profileSet triple memory
	profileSet = (struct Profile_triple *)malloc(sizeof(struct Profile_triple)*local_n);
	
	/********************para load profile dataset by openmp in every process******************************/
	#pragma omp parallel num_threads(corenum)
	{
		int local_t;	//the data number of each thread must hand
		int begin_t,end_t;
		int threadID = omp_get_thread_num();
		
		// compute the local size ¡¢up boundary and down boundary for every thread
		split_data(local_n, corenum, threadID, &begin_t, &end_t, &local_t);
		
		// compute the begin_t to end_t triples(file:begin1->begin1+len => triples: begin2->begin2+len)
		getFreeTriples(genelen, siglen, profilenum, linelen, begin+begin_t, begin_t, local_t, input, profileSet);
	}
	
	//file:begin->end => triples: 0->local_P
	//getTriples(local_n, genelen, siglen, profilenum, linelen, begin, end, input, profileSet);	
	/************************************************************************************/
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == 0){
		GET_TIME(finish);
		//compute the IO time
		duration = finish-start;     
		printf("loading IO and prework time: %.4f s\n",duration);		
	}
	
	if(my_rank == 0){		
		
		printf("which type of profile do you want to input ( 0 -> gene symbol list, others -> gene symbols and their expression levels ):");
		scanf("%d", &input_way);
		printf("input the path of file that has target profile until 'exit'(each line has a Gene Symbol/name):\n");
		scanf("%s",gsStr);
	}
	
	MPI_Bcast( gsStr, 1024, MPI_CHAR, 0 ,MPI_COMM_WORLD);  //Bcast the gsStr to all process
	MPI_Bcast( &input_way, 1, MPI_INT, 0 ,MPI_COMM_WORLD);  //Bcast the input_way to all process
	
	
	while(strcmp(gsStr,"exit")!=0)
	{
		//get the target profile
		if(input_way==0)
		{
			getProfile(target_profile,&target_profile_len,gsStr,genelistfile);
		}			
		else
		{
			getProfilewithExpression(original_target_profile,&target_profile_len, gsStr,genelistfile);
		}	
		
		if(target_profile_len != genelen)
		{
			if(my_rank==0)
			{
				getchar();    //remove the Enter from stdin
				printf("This profile is too short, please make sure it has same %d genes with L1000 library profiles!\n",genelen);
				printf("input the path of file that has target profile until 'exit'(each line has a Gene Symbol/name):\n");
				scanf("%s",gsStr);
			}
			MPI_Bcast( gsStr, 1024, MPI_CHAR, 0 ,MPI_COMM_WORLD);
			continue;
		}
			
		MPI_Barrier(MPI_COMM_WORLD);
		/********************run the GSEA algorithm by mpi****************************/
		GET_TIME(start);
		
		if(input_way != 0) //get sorted id profile
		{
			quiksort_profile(original_target_profile, 0, genelen-1);
			for(i=0; i< genelen; i++)
				target_profile[i] = original_target_profile[i].id;
		}
		searchProfile = getTriple(target_profile, genelen, siglen);
		
		//para using the GSEA algorithm for all process by multi-thread way	
		#pragma omp parallel num_threads(corenum)
		{
			int t_c;
			int local_t;	//the data number of each thread must hand
			int begin_t,end_t;
			int threadID = omp_get_thread_num();
		
			// compute the local size ¡¢up boundary and down boundary for every thread
			split_data(local_n, corenum, threadID, &begin_t, &end_t, &local_t);
				
			for(t_c=begin_t; t_c<end_t; t_c++)
			{		
				local_es[t_c].ES = ES_Profile_triple(searchProfile, profileSet[t_c], genelen, siglen);
				local_es[t_c].cid = begin+t_c+1;
			}
		}
		
		/*********************gather gsea to process0 in es_result*******************************/
		Build_derived_type(&local_es[0],&es_mpi_t); //derive the new MPI Type
		
		if(my_rank == 0){
			//tmp memory for processes before leave
			struct ES_RESULT *tmp_es1 = (struct ES_RESULT*)malloc(local_n*sizeof(struct ES_RESULT));
			//tmp memory for processes after leave
			struct ES_RESULT *tmp_es2 = (struct ES_RESULT*)malloc((local_n-1)*sizeof(struct ES_RESULT));
			//copy local es vector in process0 to global es vector 
			memcpy(es_result,local_es,local_n*sizeof(struct ES_RESULT));
			/*****************receive and copy local es vector in other processes to global es vector******/
			if(leave == 0 ){  //even process0 not increase the local_n
				for(i=1; i<p; i++){
					MPI_Recv(tmp_es1 , local_n, es_mpi_t, i, tag, MPI_COMM_WORLD, &status);
					memcpy(&es_result[i*local_n],tmp_es1,local_n*sizeof(struct ES_RESULT));	
				}
			}else{ //at last ,the process0 had a increased local_n
				for(i=1; i<p; i++){
					if(i < leave){
						MPI_Recv(tmp_es1 , local_n, es_mpi_t, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&es_result[i*local_n],tmp_es1,local_n*sizeof(struct ES_RESULT));
					}else{
						MPI_Recv(tmp_es2 , local_n-1, es_mpi_t, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&es_result[i*(local_n-1)+leave],tmp_es2,(local_n-1)*sizeof(struct ES_RESULT));
					}	
				}			
			}			
			free(tmp_es1);
			free(tmp_es2);
		}else{
			MPI_Send(local_es, local_n, es_mpi_t, 0, tag, MPI_COMM_WORLD);
		}	
		
		if(my_rank == 0){
			//sort the gsea result
			quiksort_es(es_result,0,profilenum-1);
		
			/********************print the TopN GSEA result*************************/
			printf("\nprintf the high level of TopN similar profile result:\n");
			for(i = profilenum-1; i > profilenum-1-TopN; i--)
			{
				cidnum = readByteOffsetFile(sample,es_result[i].cid);
				offset = readByteOffsetFile(offsetfile,cidnum);
				getSampleConditions(conditionsfile, offset, conditions);
				printf("\nNO.%d -> SampleConditions: %s  ES:%f\n", profilenum-i, conditions, es_result[i].ES );
			}
			printf("\nprintf the low level of TopN similar profile result:\n");
			for(i=0; i<TopN; i++)
			{
				cidnum = readByteOffsetFile(sample,es_result[i].cid);
				offset = readByteOffsetFile(offsetfile,cidnum);
				getSampleConditions(conditionsfile, offset, conditions);
				printf("\nNO.%d -> SampleConditions: %s  ES:%f\n", i+1, conditions, es_result[i].ES ); 
			}
		}
	
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0){
			GET_TIME(finish);
			//compute the ES time 
			duration = finish-start;     
			printf("finish ES time: %.4f s\n",duration);
			
			getchar();    //remove the Enter from stdin
			//get the profile path
			printf("input the path of file that has target profile until 'exit'(each line has a Gene Symbol/name):\n");
			scanf("%s",gsStr);
		}
		MPI_Bcast( gsStr, 1024, MPI_CHAR, 0 ,MPI_COMM_WORLD);  //Bcast the gsStr to all process
	}
	
	//free the memory allocate dyn.
	if(my_rank == 0)
		free(es_result);
	free(local_es);
	free(profileSet);

	MPI_Finalize();
	return 0;

}

void Build_derived_type(
		struct ES_RESULT* m_ptr, 	 /*  in  */
		MPI_Datatype* es_mpi_t_ptr /*  out  */)
{
		
	int block_lengths[2];
	MPI_Aint displacements[2];
	MPI_Datatype typelist[2];
		
	MPI_Aint start_address;
	MPI_Aint address;
		
	block_lengths[0] = 1;
	block_lengths[1] = 1;
						
	typelist[0] = MPI_FLOAT;
	typelist[1] = MPI_INT;
		
	displacements[0] = 0;
		
	MPI_Address(m_ptr, &start_address);
	MPI_Address(&(m_ptr->cid), &address);
	displacements[1] = address - start_address;
		
	MPI_Type_struct(2,block_lengths,displacements,typelist,es_mpi_t_ptr);
	MPI_Type_commit(es_mpi_t_ptr);
}

void Usage() {
	fprintf(stderr, "%s\n", USAGE);
}  /* Usage */