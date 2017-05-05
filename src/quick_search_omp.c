#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include <omp.h>
#include <unistd.h> 
#include <getopt.h> 
#include "Tools.h"

#define ERRM "quick search error:"

char *USAGE =
"\n"
"Usage:"
"  quick_search_omp [options]\n"
"\n"
"  general options:\n"
"	 -t --thread: the number of threads. [ default 1 ] \n"
"    -n --topn: The first and last N GSEA records ordered by ES. [ default 10 ]\n"
"\n"
"  input/output options: \n"
"    -i --input: input file/a parsed profiles's file from pretreatment stage. \n"
"    -s --sample: input file/a parsed sample sequence number file from pretreatment stage. \n"
"    -r --reference: input a directory includes referenced files about genesymbols and cids. \n";

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
	
	FILE *fp;
	char conditions[L1000_CONDITION_LEN];
	char conditionsfile[FileName_LEN];
	char offsetfile[FileName_LEN];
	char genelistfile[FileName_LEN];
	long cidnum;
	long offset;
	int input_way;

	int	paral_way = 1;
	
	// Unset flags (value -1).
	int TopN = -1;
	int thread_count = -1;
    // Unset options (value 'UNSET').
	char * const UNSET = "unset";
    char * input   = UNSET;
	char * sample   = UNSET;
	char * reference   = UNSET;	
	
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
			{"sample",             required_argument,        0, 's'},
			{"reference",             required_argument,        0, 'r'},
			{"input",             required_argument,        0, 'i'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "n:i:t:s:r:",
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
			
		
		case 's':
			if (sample == UNSET) 
			{
				sample = optarg;
			}
			else 
			{
				fprintf(stderr, "%s --sample set more than once\n", ERRM);
				Usage();
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
				fprintf(stderr, "%s --reference set more than once\n", ERRM);
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
		TopN = 10;
	
	if(thread_count==-1)
		thread_count = 1;
	
	if((fp=fopen(sample,"r"))==NULL)
	{
		fprintf(stderr, " [ param error : -s ] can not open sample sequence number '%s' file\n",sample);
		exit(0);
	}
	fclose(fp);
	
	sprintf(genelistfile,"%s/Gene_List.txt",reference);
	
	if((fp=fopen(genelistfile,"r"))==NULL)
	{
		fprintf(stderr, "[ param error : -r ] the reference directory may be incorrect!\n");
		exit(0);
	}
	fclose(fp);
	
	sprintf(conditionsfile,"%s/Samples_Condition.txt",reference);
	sprintf(offsetfile,"%s/Samples_RowByteOffset.txt",reference);
		
	printf("Profile Set is Loading...!\n");
	
	GET_TIME(start);
	//read file parameters
	ReadFilePara(input, &profilenum, &genelen, &linelen);	
	
	if( profilenum <= 0 || genelen <= 0)
	{
		fprintf(stderr," [ param error : -i ] this file is not exist!\n");
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
		
	while(strcmp(gsStr,"exit")!=0)
	{
		//get the geneset
		if(input_way==0)
		{
			getGeneSet(gs,&siglen,gsStr,genelistfile);
			if(siglen==0)
			{
				getchar();    //remove the Enter from stdin
				printf("There is no gene be hitted, please make sure the GeneSet have at least one Gene in Profile!\n");
				printf("input the GeneSet until 'exit'( a string of each Gene Symbol split by space ):\n");
				scanf("%[^\n]",gsStr);
				continue;
			}
		}else
		{
			getGeneSetbyFile(gs,&siglen,gsStr,genelistfile);
			if(siglen==0)
			{
				getchar();    //remove the Enter from stdin
				printf("There is no gene be hitted, please make sure the GeneSet have at least one Gene in Profile!\n");
				printf("input the path of file that has GeneSet until 'exit'(each line has a Gene Symbol/name):\n");
				scanf("%s",gsStr);
				continue;
			}
		}
		
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
				
		GET_TIME(finish);
		duration = finish-start;    //compute the GSEA time 
		printf("finish GSEA time by openmp: %.4f s\n",duration); 
		
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
