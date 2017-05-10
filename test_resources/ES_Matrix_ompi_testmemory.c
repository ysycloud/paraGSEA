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

#define ERRM "ES Matrix error:"

char *USAGE =
"\n";
 
void Usage();

int main(int argc,char *argv[])
{	
	int i,j;
	int genelen;
	int profilenum1;
	int linelen1;
	struct Profile_triple **triples1;
	int	my_rank;   /* My process rank           */
    int	p;         /* The number of processes   */
	int begin,end;
	int parameternum;
	int corenum;
	int siglen;	
	int load_time;
	float proportion;
	

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
			{"input1",             required_argument,        0, '1'},
			{"input2",             required_argument,        0, '2'},
			{"output",             required_argument,        0, 'o'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "t:l:a:p:1:2:o:",
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
	
	if(output == UNSET)
	{
		if(my_rank==0)
			fprintf(stderr," [ param error : -o ] Not Set output parameter!\n");
		MPI_Finalize();
		exit(0);
	}
	
	//barrier all processes to compute time
	MPI_Barrier(MPI_COMM_WORLD); 
	if(my_rank == 0){
		printf("Profile Set is Loading...!\n");
		GET_TIME(start);
	}
	
	//read file parameters in all processes
	ReadFilePara(input1, &profilenum1, &genelen, &linelen1);
	
	profilenum1 *= proportion;
	
	//input file check
	if( profilenum1 <= 0 || genelen <= 0)
	{
		if(my_rank==0)
			fprintf(stderr," [ param error : -1 ] this file input1 is not exist!\n");
		MPI_Finalize();
		exit(0);
	}
	
		
	
	if(my_rank == 0){
		
		triples1 = (struct Profile_triple **)malloc(sizeof(struct Profile_triple*)*30);
		
		for(i=0;i<100;i++)
		{
			/*****read the local part file of dataset1 in every process and get their triples****************/
			triples1[i] = (struct Profile_triple *)malloc(sizeof(struct Profile_triple)*profilenum1);	
			getTriples(profilenum1, genelen, siglen, profilenum1, linelen1, 0, profilenum1, input1, triples1[i]);
			
			printf("term %d\n",i );
		}

	}
	
	for(i=0;i<30;i++)
		free(triples1[i]);
	//free the memory
	free(triples1);
	
	MPI_Finalize();
	return 0;
}

void Usage() {
	fprintf(stderr, "%s\n", USAGE);
}  /* Usage */
