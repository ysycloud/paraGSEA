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
	
	// Unset flags (value -1).
	int TopN = -1;
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
	
	if((fp=fopen(sample,"r"))==NULL)
	{
		fprintf(stderr, "can not open %s file\n",sample);
		exit(0);
	}
	fclose(fp);
	
	sprintf(genelistfile,"%s/Gene_List.txt",reference);
	
	if((fp=fopen(genelistfile,"r"))==NULL)
	{
		fprintf(stderr, "can not open %s file\n",genelistfile);
		exit(0);
	}
	fclose(fp);
	
	sprintf(conditionsfile,"%s/Samples_Condition.txt",reference);
	sprintf(offsetfile,"%s/Samples_RowByteOffset.txt",reference);
	
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
		}else
		{
			getGeneSetbyFile(gs,&siglen,gsStr,genelistfile);
		}
		
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
		printf("\nprintf the high level of TopN GSEA result:\n");
		for(i = profilenum-1; i > profilenum-1-TopN; i--)
		{
			cidnum = readByteOffsetFile(sample,gsea_result[i].cid);
			offset = readByteOffsetFile(offsetfile,cidnum);
			getSampleConditions(conditionsfile, offset, conditions);
			printf("NO.%d -> SampleConditions: %s  ES:%f  NES:%f  pv:%.10lf\n", profilenum-i, conditions, gsea_result[i].ES, gsea_result[i].NES, gsea_result[i].pv);
		}
			
		printf("\nprintf the low level of TopN GSEA result:\n");
		for(i=0; i<TopN; i++)
		{
			cidnum = readByteOffsetFile(sample,gsea_result[i].cid);
			offset = readByteOffsetFile(offsetfile,cidnum);
			getSampleConditions(conditionsfile, offset, conditions);
			printf("NO.%d -> SampleConditions: %s  ES:%f  NES:%f  pv:%.10lf\n", i+1, conditions, gsea_result[i].ES, gsea_result[i].NES, gsea_result[i].pv); 
		}
			 				
		GET_TIME(finish);
		duration = finish-start;    //compute the GSEA time 
		printf("finish GSEA time: %.4f s\n",duration); 
		
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