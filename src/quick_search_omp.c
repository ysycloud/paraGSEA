#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include <omp.h>
#include <unistd.h> 
#include <getopt.h> 
#include "RandomChange.h"
#include "GSEA.h"
#include "IO.h"

#define ERRM "quick search error:"

char *USAGE =
"\n"
"Usage:"
"  quick_search_omp [options]\n"
"\n"
"  general options:\n"
"	 -t --thread: the number of threads\n"
"    -n --topn: The first and last N GSEA records ordered by ES\n"
"\n"
"  input/output options \n"
"    -i --input: input file/a parsed profiles's file from pretreatment stage. \n";

float global_ES[Global_ES_SIZE];

void Usage();

int main(int argc,char *argv[])
{	
	int i,j,profilenum,genelen,linelen,siglen;
	short **profileSet;
	short **indexSet;
	short gs[MAX_GENESET];
	char gsStr[1024];
	struct GSEA_RESULT *gsea_result;
	
	double start,finish,duration;
	
	int input_way;
	int	paral_way = 1;
	
	// Unset flags (value -1).
	int TopN = -1;
	int thread_count = -1;
    // Unset options (value 'UNSET').
	char * const UNSET = "unset";
    char * input   = UNSET;

	
	if (argc == 1) 
	{
		Usage();
		exit(0);
    }
	
	int c;
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"thread",            required_argument,        0, 't'},
			{"topn",              required_argument,        0, 'n'},
			{"input",             required_argument,        0, 'i'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "n:i:t:",
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
				fprintf(stderr, "%s --input set more than once\n", ERRM);
				Usage();
				exit(0);
			}
			break;
		
		case 'n':
			if (TopN < 0) {
				TopN = atoi(optarg);
				if (TopN < 1) {
					fprintf(stderr, "%s --topn must be a positive integer\n", ERRM);
					Usage();
					exit(0);
				}
			}
			else {
				fprintf(stderr,"%s --topn set more " "than once\n", ERRM);
				Usage();
				exit(0);
			}
			break;
			
		case 't':
			if (thread_count < 0) {
				thread_count = atoi(optarg);
				if (thread_count < 1) {
					fprintf(stderr, "%s --thread must be a positive integer\n", ERRM);
					Usage();
					exit(0);
				}
			}
			else {
				fprintf(stderr,"%s --thread set more " "than once\n", ERRM);
				Usage();
				exit(0);
			}
			break;
			
		default:
			// Cannot parse. //
			Usage();
			exit(0);
		}		
	}

	//check the parameters
	if(TopN==-1)
	{
		fprintf(stderr,"Not Set TopN parameter!\n");
		exit(0);
	}
	
	if(thread_count==-1)
	{
		fprintf(stderr,"Not Set Thread parameter!\n");
		exit(0);
	}
		
	printf("Profile Set is Loading...!\n");
	
	GET_TIME(start);
	//read file parameters
	ReadFilePara(input, &profilenum, &genelen, &linelen);	
	
	if( profilenum <= 0 || genelen <= 0)
	{
		fprintf(stderr,"this file is not exist!\n");
		exit(0);
	}
	
	printf("profilenum:%d\t genelen:%d\n",profilenum,genelen);
	
	//malloc profile dataset memory
	profileSet = (short **)malloc(profilenum*sizeof(short *));
	for(i=0;i<profilenum;i++)
		profileSet[i] = (short *)malloc(genelen*sizeof(short));
	//malloc index set for profile dataset 
	indexSet = (short **)malloc(profilenum*sizeof(short *));
	for(i=0;i<profilenum;i++)
		indexSet[i] = (short *)malloc(genelen*sizeof(short));
	
	//malloc GSEA para Vector
	gsea_result = (struct GSEA_RESULT*)malloc(profilenum*sizeof(struct GSEA_RESULT));
	
	/********************para load profile dataset by openmp******************************/
	#pragma omp parallel num_threads(thread_count)
	{
		int k;
		
		if(paral_way == 1){
			#pragma omp for
			for( k=0; k<profilenum; k++){
				ReadFile(input, linelen, k, k+1, profilenum, genelen, profileSet);
				getIndex(profileSet[k],indexSet[k],genelen);
			}
		}else{
			int local_n;	//the data number of each thread must hand
			int leave;		//the leave data number come from the data can not be divided totally by thread number
			local_n = profilenum / thread_count;  
			leave = profilenum % thread_count;
			int threadID = omp_get_thread_num();
			int begin,end;
			// compute the up boundary and down boundary for every thread	
			if(threadID < leave){
				local_n++;   
				begin = threadID*local_n;
			}else{	
				begin = threadID*local_n + leave; 
			}	 
			end = begin + local_n;
		
			//para read the file to global profileSet memory
			ReadFile(input, linelen, begin, end, profilenum, genelen, profileSet);
		
			//para compute the index for profile sets
			for(k=begin; k<end; k++)
				getIndex(profileSet[k],indexSet[k],genelen);
		}		
	}
	
	GET_TIME(finish);
	//compute the IO time
	duration = finish-start;     
	printf("loading IO and prework time by openmp: %.4f s\n",duration); 
	
		
	//get the geneset , split by space
	printf("input the GeneSet( a integer[1-genelen] string split by space ):\n");
	scanf("%[^\n]",gsStr);	
	//gets(gsStr);
	while(strcmp(gsStr,"exit")!=0)
	{
		//get the geneset
		getGeneSet(gs,&siglen,gsStr);
		
		GET_TIME(start);
		/********************run the GSEA algorithm by openmp*****************************/
		//compute the global ES	
		getGlobalES( genelen, siglen , global_ES);
		
		#pragma omp parallel num_threads(thread_count)
		{
			int k;
			
			if(paral_way == 1){
				#pragma omp for
				for(k=0; k<profilenum; k++){
					GSEA( gs, indexSet[k], genelen, siglen, &(gsea_result[k].ES), &(gsea_result[k].NES), &(gsea_result[k].pv), global_ES );
					gsea_result[k].cid = k+1;
				}
			}else{
				int local_n;	//the data number of each thread must hand
				int leave;		//the leave data number come from the data can not be divided totally by thread number
				local_n = profilenum / thread_count;  
				leave = profilenum % thread_count;
				int threadID=omp_get_thread_num();
				int begin,end;
				// compute the up boundary and down boundary for every thread
				if(threadID < leave){
					local_n++;   
					begin = threadID*local_n;
				}else{	
					begin = threadID*local_n + leave; 
				}	 
				end = begin + local_n;
				//para using the GSEA algorithm
				for(k=begin; k<end; k++){
					GSEA( gs, indexSet[k], genelen, siglen, &(gsea_result[k].ES), &(gsea_result[k].NES), &(gsea_result[k].pv), global_ES );
					gsea_result[k].cid = k+1;
				}
			}
		}
		
		//sort the gsea result
		quiksort_gsea(gsea_result,0,profilenum-1);
		
		/********************print the TopN GSEA result*************************/
		printf("printf the high level of TopN GSEA result:\n");
		for(i = profilenum-1; i > profilenum-1-TopN; i--)
			printf("NO.%d -> cid:%d  ES:%f  NES:%f  pv:%.10lf\n", profilenum-i, gsea_result[i].cid, gsea_result[i].ES, gsea_result[i].NES, gsea_result[i].pv);
		printf("printf the low level of TopN GSEA result:\n");
		for(i=0; i<TopN; i++)
			printf("NO.%d -> cid:%d  ES:%f  NES:%f  pv:%.10lf\n", i+1, gsea_result[i].cid, gsea_result[i].ES, gsea_result[i].NES, gsea_result[i].pv);  
				
		GET_TIME(finish);
		duration = finish-start;    //compute the GSEA time 
		printf("finish GSEA time by openmp: %.4f s\n",duration); 
		
		getchar();    //remove the Enter from stdin
		printf("input the GeneSet(split by space):\n");
		scanf("%[^\n]",gsStr);
		//gets(gsStr);		
	}
	
	//free the memory allocate dyn.
	free(gsea_result);
	for(i=0; i<profilenum; i++){
		free(profileSet[i]);
		free(indexSet[i]);
	}
	free(profileSet);
	free(indexSet);
	
	return 0;
}

void Usage() {
	fprintf(stderr, "%s\n", USAGE);
}  /* Usage */
