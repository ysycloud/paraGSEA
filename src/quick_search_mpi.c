#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include <unistd.h> 
#include <getopt.h> 
#include "mpi.h"
#include "Tools.h"

#define ERRM "quick search error:"

char *USAGE =
"\n"
"Usage:"
"  quick_search_mpi [options]\n"
"\n"
"  general options before command by MPI:\n"
"	 -n process_num : Total number of processes. [ default 1 ]\n"
"	 -ppn pernum: the number of processes in each node. [ default 1 ]\n"
"	 -hostfile hostfile:  list the IP or Hostname of nodes. [ default localhost ]"
"\n"
"  general options:\n"
"    -n --topn: The first and last N GSEA records ordered by ES. [ default 10 ]\n"
"\n"
"  input/output options: \n"
"    -i --input: input file/a parsed profiles's file from pretreatment stage. \n"
"    -s --sample: input file/a parsed sample sequence number file from pretreatment stage. \n"
"    -r --reference: input a directory includes referenced files about genesymbols and cids. \n";

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
			{"input",             required_argument,        0, 'i'},
			{"sample",             required_argument,        0, 's'},
			{"reference",             required_argument,        0, 'r'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "n:i:s:r:",
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
	unsigned long memneed = (2*sizeof(short)*(profilenum/p+1)*genelen + profilenum*sizeof(struct GSEA_RESULT))/1024;
	
	if(my_rank == 0 )
	{
		printf("Available Memory:      %ld KB\n", memavail);		
		printf("Needed Memory:      %ld KB\n", memneed);
		
		if(memavail < memneed)
			printf("available memory is not enough!!! Please use more nodes!!!\n");
		else
			//malloc global GSEA para Vector for process0
			gsea_result = (struct GSEA_RESULT*)malloc(profilenum*sizeof(struct GSEA_RESULT));
	}
	
	if(memavail < memneed)	
		return;
	
			
	// compute the local size ��up boundary and down boundary for every process	
	int begin,end;
	int leave = profilenum % p;
	split_data(profilenum, p, my_rank, &begin, &end, &local_n);
	
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
		
		printf("which way do you want to input the GeneSet( 0 -> standard input , others -> file input ):");
		scanf("%d", &input_way);
	
		if(input_way==0)
		{
			//get the geneset , split by space
			getchar();
			printf("input the GeneSet until 'exit'( a string of each Gene Symbol split by space ):\n");
			scanf("%[^\n]",gsStr);	
		}else
		{
			printf("input the path of file that has GeneSet until 'exit'(each line has a Gene Symbol/name):\n");
			scanf("%s",gsStr);
		}
	}
	
	MPI_Bcast( gsStr, 1024, MPI_CHAR, 0 ,MPI_COMM_WORLD);  //Bcast the gsStr to all process
	MPI_Bcast( &input_way, 1, MPI_INT, 0 ,MPI_COMM_WORLD);  //Bcast the input_way to all process
		
	while(strcmp(gsStr,"exit")!=0)
	{
		//get the geneset
		if(input_way==0)
		{
			getGeneSet(gs,&siglen,gsStr,genelistfile);
			if(siglen==0)
			{
				if(my_rank==0)
				{
					getchar();    //remove the Enter from stdin
					printf("There is no gene be hitted, please make sure the GeneSet have at least one Gene in Profile!\n");
					printf("input the GeneSet until 'exit'( a string of each Gene Symbol split by space ):\n");
					scanf("%[^\n]",gsStr);
				}
				MPI_Bcast( gsStr, 1024, MPI_CHAR, 0 ,MPI_COMM_WORLD);
				continue;
			}
		}else
		{
			getGeneSetbyFile(gs,&siglen,gsStr,genelistfile);
			if(siglen==0)
			{
				if(my_rank==0)
				{
					getchar();    //remove the Enter from stdin
					printf("There is no gene be hitted, please make sure the GeneSet have at least one Gene in Profile!\n");
					printf("input the path of file that has GeneSet until 'exit'(each line has a Gene Symbol/name):\n");
					scanf("%s",gsStr);
				}
				MPI_Bcast( gsStr, 1024, MPI_CHAR, 0 ,MPI_COMM_WORLD);
				continue;
			}
		}
	
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
			printf("\nprintf the high level of TopN GSEA result:\n");
			for(i = profilenum-1; i > profilenum-1-TopN; i--)
			{
				cidnum = readByteOffsetFile(sample,gsea_result[i].cid);
				offset = readByteOffsetFile(offsetfile,cidnum);
				getSampleConditions(conditionsfile, offset, conditions);
				printf("\nNO.%d -> SampleConditions: %s  ES:%f  NES:%f  pv:%.10lf\n", profilenum-i, conditions, gsea_result[i].ES, gsea_result[i].NES, gsea_result[i].pv);
			}
			printf("\nprintf the low level of TopN GSEA result:\n");
			for(i=0; i<TopN; i++)
			{
				cidnum = readByteOffsetFile(sample,gsea_result[i].cid);
				offset = readByteOffsetFile(offsetfile,cidnum);
				getSampleConditions(conditionsfile, offset, conditions);
				printf("\nNO.%d -> SampleConditions: %s  ES:%f  NES:%f  pv:%.10lf\n", i+1, conditions, gsea_result[i].ES, gsea_result[i].NES, gsea_result[i].pv); 
			}
		}
	
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0){
			GET_TIME(finish);
			//compute the GSEA time 
			duration = finish-start;     
			printf("finish GSEA time by MPI: %.4f s\n",duration);
			
			getchar();    //remove the Enter from stdin
			//get the geneset
			if(input_way==0)
			{
				//get the geneset , split by space
				printf("input the GeneSet until 'exit'( a string of each Gene Symbol split by space ):\n");
				scanf("%[^\n]",gsStr);	
			}else
			{
				printf("input the path of file that has GeneSet until 'exit'(each line has a Gene Symbol/name):\n");
				scanf("%s",gsStr);
			}	
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