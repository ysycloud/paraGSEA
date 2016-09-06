#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include <unistd.h> 
#include <getopt.h> 
#include "RandomChange.h"
#include "GSEA.h"
#include "IO.h"

#define ERRM "quick search error:"

char *USAGE =
"\n"
"Usage:"
"  quick_search_serial [options]\n"
"\n"
"  general options:\n"
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
	

	// Unset flags (value -1).
	int TopN = -1;
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
	
	printf("Profile Set is Loading...!\n");
	
	GET_TIME(start);
	//read file parameters
	ReadFilePara(input, &profilenum, &genelen, &linelen);
	
	if( profilenum <= 0 || genelen <= 0 )
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
	
	//load profile dataset
	ReadFile(input, linelen, 0 , profilenum , profilenum, genelen, profileSet); 
		
	//compute the index for profile sets
	for(i=0; i<profilenum; i++)
		getIndex(profileSet[i],indexSet[i],genelen);
	
	GET_TIME(finish);
	//compute the IO time and prework time
	duration = finish-start;     
	printf("loading IO and prework time: %.4f s\n",duration); 
	
	
	//get the geneset , split by space
	printf("input the GeneSet( a integer[1-genelen] string split by space ):\n");
	scanf("%[^\n]",gsStr);	
	//gets(gsStr);
	while(strcmp(gsStr,"exit")!=0)
	{
		//get the geneset
		getGeneSet(gs,&siglen,gsStr);
		
		GET_TIME(start);
		/********************run the GSEA algorithm*****************************/
		//compute the global ES	
		getGlobalES( genelen, siglen , global_ES);
		
		for(i=0; i<profilenum; i++){
			GSEA( gs, indexSet[i], genelen, siglen, &(gsea_result[i].ES), &(gsea_result[i].NES), &(gsea_result[i].pv), global_ES );
			gsea_result[i].cid = i+1;
		}
		
		//printf("cid:%d  ES:%f  NES:%f  pv:%.10lf\n",gsea_result[19999].cid, gsea_result[19999].ES, gsea_result[19999].NES, gsea_result[19999].pv);
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
		printf("finish GSEA time: %.4f s\n",duration); 
		
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