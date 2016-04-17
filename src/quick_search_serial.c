#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include "RandomChange.h"
#include "GSEA.h"
#include "IO.h"

float global_ES[Global_ES_SIZE];

void Usage(char prog_name[]);

int main(int argc,char *argv[])
{	
	int i,j,profilenum,genelen,linelen,siglen;
	short **profileSet;
	short **indexSet;
	short gs[MAX_GENESET];
	char gsStr[1024];
	struct GSEA_RESULT *gsea_result;	
	double start,finish,duration;
	
	if(argc!=3)
	{
		Usage(argv[0]);
		exit(0);
	}
	
	int	TopN = atoi(argv[2]);
				
	printf("Profile Set is Loading...!\n");
	
	GET_TIME(start);
	//read file parameters
	ReadFilePara(argv[1], &profilenum, &genelen, &linelen);	
	
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
	ReadFile(argv[1], linelen, 0 , profilenum , profilenum, genelen, profileSet); 
		
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
		for(i=0;i<siglen;i++)
			printf("%d\t",gs[i]);
		printf("\n");
		
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

void Usage(char prog_name[]) {
	fprintf(stderr, "usage:  %s <inputfile> <TopN>\n", prog_name);
}  /* Usage */