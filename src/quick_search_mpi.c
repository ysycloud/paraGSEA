#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include <unistd.h> 
#include <getopt.h> 
#include "mpi.h"
#include "RandomChange.h"
#include "GSEA.h"
#include "IO.h"

#define ERRM "quick search error:"

char *USAGE =
"\n"
"Usage:"
"  quick_search_mpi [options]\n"
"\n"
"  general options before command by MPI:\n"
"	 -n process_num : Total number of processes\n"
"	 -ppn pernum: the number of processes in each node\n"
"	 -hostfile hostfile:  list the IP or Hostname of nodes"
"\n"
"  general options:\n"
"    -n --topn: The first and last N GSEA records ordered by ES\n"
"\n"
"  input/output options \n"
"    -i --input: input file/a parsed profiles's file from pretreatment stage. \n";

void Usage();
void Build_derived_type(
		struct GSEA_RESULT* m_ptr, 			 /*  in  */
		MPI_Datatype* gsea_mpi_t_ptr 		 /*  out  */);
 
int main(int argc,char *argv[])
{	
	int i,profilenum,genelen,linelen,siglen;
	short **profileSet,**tmp_global_profile;
	short **indexSet;
	char gsStr[1024];
	short gs[MAX_GENESET];	
	struct GSEA_RESULT *gsea_result;
	struct GSEA_RESULT *local_gsea;
	float global_ES[Global_ES_SIZE];
	MPI_Datatype gsea_mpi_t;
	int	my_rank;   /* My process rank           */
    int	p;         /* The number of processes   */
    int source,dest ;  
    int tag = 0;
    MPI_Status  status;
	int local_n;	//the data number of each processes must hand
	int leave;		//the leave data number come from the data can not be divided totally by process number
	int parameternum;
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
	int TopN = -1;
    // Unset options (value 'UNSET').
	char * const UNSET = "unset";
    char * input   = UNSET;
	
	int c;
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"topn",              required_argument,        0, 'n'},
			{"input",             required_argument,        0, 'i'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "n:i:",
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
					exit(0);
				}
			}
			else {
				if(my_rank==0)
				{
					fprintf(stderr,"%s --topn set more " "than once\n", ERRM);
					Usage();
				}		
				exit(0);
			}
			break;
			
		default:
			// Cannot parse. //
			if(my_rank==0)
				Usage();
			exit(0);
		}		
	}

	//check the parameters
	if(TopN==-1)
	{
		if(my_rank==0)
			fprintf(stderr,"Not Set TopN parameter!\n");
		exit(0);
	}

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
			fprintf(stderr,"this file is not exist!\n");
		exit(0);
	}
	
	if(my_rank == 0){
		printf("profilenum:%d\t genelen:%d\n",profilenum,genelen);
		//malloc global GSEA para Vector for process0
		gsea_result = (struct GSEA_RESULT*)malloc(profilenum*sizeof(struct GSEA_RESULT));
	}
			
	// compute the local size ¡¢up boundary and down boundary for every process	
	local_n = profilenum / p;  
	leave = profilenum % p;
	int begin,end;
	if(my_rank < leave){
		local_n++;   
		begin = my_rank*local_n;
	}else{	
		begin = my_rank*local_n + leave; 
	}	 
	end = begin + local_n;
	
	//malloc profile dataset memory for all processes using local size
	profileSet = (short **)malloc(local_n*sizeof(short *));
	for(i=0;i<local_n;i++)
		profileSet[i] = (short *)malloc(genelen*sizeof(short));
	//malloc index set for profile dataset 
	indexSet = (short **)malloc(local_n*sizeof(short *));
	for(i=0;i<local_n;i++)
		indexSet[i] = (short *)malloc(genelen*sizeof(short));
	local_gsea = (struct GSEA_RESULT*)malloc(local_n*sizeof(struct GSEA_RESULT));
	//malloc tmp all profiles dataset memory for all processes in order to read
	tmp_global_profile = (short **)malloc(profilenum*sizeof(short *));
	for(i=0;i<profilenum;i++)
		tmp_global_profile[i] = (short *)malloc(genelen*sizeof(short));
	
	/********************para load profile dataset by openmp******************************/
	//para read the file to local profileSet memory
	ReadFile(input, linelen, begin, end, profilenum, genelen, tmp_global_profile);
	for(i=0; i<local_n; i++)
		memcpy(profileSet[i],tmp_global_profile[i+begin],genelen*sizeof(short));
	
	//para compute the index for profile sets
	for(i=0; i<local_n; i++)
		getIndex(profileSet[i],indexSet[i],genelen);
	/************************************************************************************/
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == 0){
		GET_TIME(finish);
		//compute the IO time
		duration = finish-start;     
		printf("loading IO and prework time by MPI: %.4f s\n",duration);		
	}
	
	if(my_rank == 0){		
		//get the geneset , split by space
		printf("input the GeneSet( a integer[1-genelen] string split by space ):\n");
		scanf("%[^\n]",gsStr);
		getchar();
	}
	
	MPI_Bcast( gsStr, 1024, MPI_CHAR, 0 ,MPI_COMM_WORLD);  //Bcast the gsStr to all process
		
	while(strcmp(gsStr,"exit")!=0)
	{
		//get the geneset
		getGeneSet(gs,&siglen,gsStr);
	
		MPI_Barrier(MPI_COMM_WORLD);
		/********************run the GSEA algorithm by mpi****************************/
		GET_TIME(start);
		
		//compute the global ES	in all process
		getGlobalES( genelen, siglen , global_ES);

		//para using the GSEA algorithm for all process
		for(i=0; i<local_n; i++){
			GSEA( gs, indexSet[i], genelen, siglen, &(local_gsea[i].ES), &(local_gsea[i].NES), &(local_gsea[i].pv), global_ES );
			local_gsea[i].cid = begin+i+1;
		}

		/*********************gather gsea to process0 in gsea_result*******************************/
		Build_derived_type(&local_gsea[0],&gsea_mpi_t); //derive the new MPI Type
		if(my_rank == 0){
			//tmp memory for processes before leave
			struct GSEA_RESULT *tmp_gsea1 = (struct GSEA_RESULT*)malloc(local_n*sizeof(struct GSEA_RESULT));
			//tmp memory for processes after leave
			struct GSEA_RESULT *tmp_gsea2 = (struct GSEA_RESULT*)malloc((local_n-1)*sizeof(struct GSEA_RESULT));
			//copy local gsea vector in process0 to global gsea vector 
			memcpy(gsea_result,local_gsea,local_n*sizeof(struct GSEA_RESULT));
			/*****************receive and copy local gsea vector in other processes to global gsea vector******/
			if(leave == 0 ){  //even process0 not increase the local_n
				for(i=1; i<p; i++){
					MPI_Recv(tmp_gsea1 , local_n, gsea_mpi_t, i, tag, MPI_COMM_WORLD, &status);
					memcpy(&gsea_result[i*local_n],tmp_gsea1,local_n*sizeof(struct GSEA_RESULT));	
				}
			}else{ //at last ,the process0 had a increased local_n
				for(i=1; i<p; i++){
					if(i < leave){
						MPI_Recv(tmp_gsea1 , local_n, gsea_mpi_t, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&gsea_result[i*local_n],tmp_gsea1,local_n*sizeof(struct GSEA_RESULT));
					}else{
						MPI_Recv(tmp_gsea2 , local_n-1, gsea_mpi_t, i, tag, MPI_COMM_WORLD, &status);
						memcpy(&gsea_result[i*(local_n-1)+leave],tmp_gsea2,(local_n-1)*sizeof(struct GSEA_RESULT));
					}	
				}			
			}			
			free(tmp_gsea1);
			free(tmp_gsea2);
		}else{
			MPI_Send(local_gsea, local_n, gsea_mpi_t, 0, tag, MPI_COMM_WORLD);
		}	
		
		if(my_rank == 0){
			//sort the gsea result
			quiksort_gsea(gsea_result,0,profilenum-1);
		
			/********************print the TopN GSEA result*************************/
			printf("printf the high level of TopN GSEA result:\n");
			for(i = profilenum-1; i > profilenum-1-TopN; i--)
				printf("NO.%d -> cid:%d  ES:%f  NES:%f  pv:%.10lf\n", profilenum-i, gsea_result[i].cid, gsea_result[i].ES, gsea_result[i].NES, gsea_result[i].pv);
			printf("printf the low level of TopN GSEA result:\n");
			for(i=0; i<TopN; i++)
				printf("NO.%d -> cid:%d  ES:%f  NES:%f  pv:%.10lf\n", i+1, gsea_result[i].cid, gsea_result[i].ES, gsea_result[i].NES, gsea_result[i].pv);  
		}
	
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0){
			GET_TIME(finish);
			//compute the GSEA time 
			duration = finish-start;     
			printf("finish GSEA time by MPI: %.4f s\n",duration);
			
			//get the geneset , split by space
			printf("input the GeneSet(split by space):\n");
			scanf("%[^\n]",gsStr);
			getchar();
		}
		MPI_Bcast( gsStr, 1024, MPI_CHAR, 0 ,MPI_COMM_WORLD);  //Bcast the gsStr to all process
	}
	
		
	//free the memory allocate dyn.
	if(my_rank == 0)
		free(gsea_result);
	free(local_gsea);
	for(i=0; i<local_n; i++){
		free(profileSet[i]);
		free(indexSet[i]);
	}
	free(profileSet);
	free(indexSet);
	for(i=0; i<profilenum; i++){
		free(tmp_global_profile[i]);
	}
	free(tmp_global_profile);

	MPI_Finalize();
	return 0;

}

void Build_derived_type(
		struct GSEA_RESULT* m_ptr, 	 /*  in  */
		MPI_Datatype* gsea_mpi_t_ptr /*  out  */)
{
		
	int block_lengths[4];
	MPI_Aint displacements[4];
	MPI_Datatype typelist[4];
		
	MPI_Aint start_address;
	MPI_Aint address;
		
	block_lengths[0] = 1;
	block_lengths[1] = 1;
	block_lengths[2] = 1;
	block_lengths[3] = 1;
						
	typelist[0] = MPI_FLOAT;
	typelist[1] = MPI_FLOAT;
	typelist[2] = MPI_DOUBLE;
	typelist[3] = MPI_INT;
		
	displacements[0] = 0;
		
	MPI_Address(m_ptr, &start_address);
	MPI_Address(&(m_ptr->NES), &address);
	displacements[1] = address - start_address;
	MPI_Address(&(m_ptr->pv), &address);
	displacements[2] = address - start_address;
	MPI_Address(&(m_ptr->cid), &address);
	displacements[3] = address - start_address;
		
	MPI_Type_struct(4,block_lengths,displacements,typelist,gsea_mpi_t_ptr);
	MPI_Type_commit(gsea_mpi_t_ptr);
}

void Usage() {
	fprintf(stderr, "%s\n", USAGE);
}  /* Usage */