#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include <omp.h>
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
	
	if(argc!=5)
	{
		Usage(argv[0]);
		exit(0);
	}

	int thread_count = atoi(argv[2]);
	int	TopN = atoi(argv[3]);
	int	paral_way = atoi(argv[4]);
		
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
	
	/********************para load profile dataset by openmp******************************/
	#pragma omp parallel num_threads(thread_count)
	{
		int k;
		
		if(paral_way == 1){
			#pragma omp for
			for( k=0; k<profilenum; k++){
				ReadFile(argv[1], linelen, k, k+1, profilenum, genelen, profileSet);
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
			ReadFile(argv[1], linelen, begin, end, profilenum, genelen, profileSet);
		
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
		for(i=0;i<siglen;i++)
			printf("%d\t",gs[i]);
		printf("\n");
		
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

void Usage(char prog_name[]) {
	fprintf(stderr, "usage:  %s <inputfile> <TopN> <Thread_num> <ParaWay>\n", prog_name);
	fprintf(stderr, " <ParaWay>: 0->split data in balance load way by ourselves; 1->using #pragma omp for\n");
}  /* Usage */
