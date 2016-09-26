#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "IO.h"
#include "RandomChange.h"
#include "GSEA.h"

//read profile dataset and get the relative parameters
int ReadFilePara(char path[], int *profilenum, int *genelen, int *LineLength)
{
	FILE *fp; 
	char StrLine[135168]; 
	int strlength;           
	char c[] = " ";
	int line,col;
	
	int i,j;
	
	//read File
	if((fp = fopen(path,"r")) == NULL) 
	{ 
		printf("file %s error!\n",path); 
		return -1; 
	} 
	fgets(StrLine,256,fp);
	*profilenum = atoi(strtok(StrLine,c));
	*genelen = atoi(strtok(NULL,c));
	
	if(*genelen>1000)
		strlength = 135168; //level3
	else
		strlength = 6144; //L1000

	fgets(StrLine,strlength,fp);
	*LineLength = strlen(StrLine);
	fclose(fp);
	return 0;
}

//read part of profile file and load to profileSet array
int ReadFile(char path[],int LineLength,int BeginLine,int EndLine,int profilenum, int genelen, short **profileSet)
{ 
	FILE *fp; 
	char *StrLine; 
	char *saveptr;
	char c[] = "\t";
	int line,col;
	int i,j;
	
	StrLine = (char *)malloc((LineLength+1)*sizeof(char));
	
	//read File
	if((fp = fopen(path,"r")) == NULL) 
	{ 
		printf("file %s error!\n", path); 
		return -1; 
	} 
	fgets(StrLine,256,fp);
	line = BeginLine;
	fseek(fp,BeginLine*(LineLength),SEEK_CUR); 
	while (!feof(fp)) 
    { 
		fgets(StrLine,LineLength+1,fp);  //read one line
		col = 0;
		if( line<profilenum && col< genelen )  
			profileSet[line][col++] = atoi(strtok_r(StrLine,c,&saveptr));
			//profileSet[line][col++] = atoi(strsep(&StrLine,c));
		char *p = strtok_r(NULL,c,&saveptr);
		//char *p = strsep(&StrLine,c);
		while(p)
		{
			if( line< profilenum &&col< genelen )
				profileSet[line][col++] = atoi(p);
			p = strtok_r(NULL,c,&saveptr);
			//p = strsep(&StrLine,c);			
		}
		line++;
		if(line>=EndLine) break;
	}
	fclose(fp);
	free(StrLine);
	return 0;
}

//write part of ES_Matrix to txt file£¬can be open by excel...
void WritetxtResult(int sourceBegin ,int sourceEnd, int matlen, char writepath[], float **ES_Matrix)
{
	int i,j;
	int local_sourcelen = sourceEnd - sourceBegin;
	
	FILE *fp;
	if((fp=fopen(writepath,"w"))==NULL)
	{
		printf("can not open the %s file to write!\n",writepath);
		exit(1);
	}
	
	fprintf(fp,"%10d\t%10d\n",local_sourcelen,matlen);
	
	for( i=sourceBegin; i < sourceEnd; i++){
		for( j=0; j < matlen-1; j++)
			fprintf(fp,"%5.3f\t",ES_Matrix[i][j]);
		fprintf(fp,"%5.3f\n",ES_Matrix[i][j]);
	}	
	fclose(fp);
}

//get GeneSet from standard input
void getGeneSet(short gs[],int *count, char gsStr[])
{
	short existflag[MAX_GENE];
	short gstmp[MAX_GENESET];
	char genelist[L1000_LEN][12];
	char c[] = " ", *p;
	char gene[12];
	int line;
	int tmp,i,j;
	*count = 0;
	
	//initial flag vector
	memset(existflag, 0, MAX_GENE * sizeof(short));
	readGeneListFile(genelist ,&line,"data/prepareForNewDataSet/Gene_List.txt");

	//get gstmp and gs count,remove the repeat elements
	strcpy(gene,strtok(gsStr,c));
	j=0;
	while(gene[j++]!='\0');
	gene[j-1]='\n';
	gene[j]='\0';
	
	for(i=0;i<line;i++)
	{
		if(strcmp(gene,genelist[i])==0)
		{
			gstmp[(*count)++] = i+1;
			existflag[gstmp[(*count)-1]] = 1;
			break;
		}
	}
	
	p = strtok(NULL,c);	
	while(p)
	{
		strcpy(gene,p);
		j=0;
		while(gene[j++]!='\0');
		gene[j-1]='\n';
		gene[j]='\0';
		for(i=0;i<line;i++)
		{
			if(strcmp(gene,genelist[i])==0)
			{
				if(existflag[i+1]==0)
				{	//this gene not input
					gstmp[(*count)++] = i+1;
					existflag[gstmp[(*count)-1]] = 1;
					break;
				}
			}
		}
		p = strtok(NULL,c);
	}
	
	//get gs 
	memcpy(gs,gstmp,(*count)*sizeof(short));
}

//get GeneSet from text file input
void getGeneSetbyFile(short gs[],int *count, char filename[])
{
	short existflag[MAX_GENE];
	short gstmp[MAX_GENESET];
	char genelist[L1000_LEN][12];
	char genesetlist[MAX_GENESET][12];
	char gene[12];
	int line1,line2;
	int tmp,i,j;
	*count = 0;
	
	//initial flag vector
	memset(existflag, 0, MAX_GENE * sizeof(short));
	readGeneListFile(genelist, &line1, "data/prepareForNewDataSet/Gene_List.txt");
	readGeneListFile(genesetlist, &line2, filename);

	//get gstmp and gs count,remove the repeat elements
	for(j=0;j<line2;j++)
	{
		for(i=0;i<line1;i++)
		{
			if(strcmp(genesetlist[j],genelist[i])==0)
			{
				if(existflag[i+1]==0)
				{	//this gene not input
					gstmp[(*count)++] = i+1;
					existflag[gstmp[(*count)-1]] = 1;
					break;
				}
			}
		}
	}
	//get gs 
	memcpy(gs,gstmp,(*count)*sizeof(short));
}

//read Matrix txt dataset and get the relative parameters
int ReadMatrixFilePara(char path[], int *profilenum1, int *profilenum2, int *LineLength)
{
	FILE *fp; 
	char Str[256];
	char *StrLine; 
	int strlength;           
	char c[] = " ";
	int line,col;
	
	int i,j;
	
	//read File
	if((fp = fopen(path,"r")) == NULL) 
	{ 
		printf("file %s error!\n",path); 
		return -1; 
	} 
	fgets(Str,256,fp);
	*profilenum1 = atoi(strtok(Str,c));
	*profilenum2 = atoi(strtok(NULL,c));
	
	strlength = ((*profilenum2)+1)*6; 
	StrLine = (char*)malloc(strlength*sizeof(char));

	fgets(StrLine,strlength,fp);
	*LineLength = strlen(StrLine);
	//printf("%d__%s\n",strlength,StrLine);
	fclose(fp);
	return 0;
}

//read part of Matrix txt file and load to matrix array
int ReadMatrixFile(char path[],int LineLength,int BeginLine,int EndLine,int profilenum1, int profilenum2, float **Matrix)
{ 
	FILE *fp; 
	char *StrLine; 
	char *saveptr;
	char c[] = "\t";
	int line,col;
	int i,j;
	
	StrLine = (char *)malloc((LineLength+1)*sizeof(char));
	
	//read File
	if((fp = fopen(path,"r")) == NULL) 
	{ 
		printf("file %s error!\n",path); 
		return -1; 
	} 
	fgets(StrLine,256,fp);
	line = BeginLine;
	fseek(fp,BeginLine*(LineLength),SEEK_CUR); 
	while (!feof(fp)) 
    { 
		fgets(StrLine,LineLength+1,fp);  //read one line
		col = 0;
		if( line<profilenum1 && col< profilenum2 )  
			Matrix[line][col++] = atof(strtok_r(StrLine,c,&saveptr));
		char *p = strtok_r(NULL,c,&saveptr);
		while(p)
		{
			if( line< profilenum1 &&col< profilenum2 )
				Matrix[line][col++] = atof(p);
			p = strtok_r(NULL,c,&saveptr);		
		}
		line++;
		if(line>=EndLine) break;
	}
	fclose(fp);
	free(StrLine);
	return 0;
}

void readGeneListFile(char genelist[][12] ,int *line, char path[])
{
	FILE *fp;
	if((fp=fopen(path,"r"))==NULL)
	{
		printf("can not open %s file\n",path);
		exit(1);
	}
	*line=0;

	while(fgets(genelist[(*line)++], 20, fp)!= NULL);
	(*line)--;

	fclose(fp);
}

void getByteOffsetFile(char path1[],char path2[])
{
	FILE *fp1,*fp2;
	char str[100];
	if((fp1=fopen(path1,"r"))==NULL)
	{
		printf("can not open %s file\n",path1);
		exit(1);
	}
	if((fp2=fopen(path2,"w"))==NULL)
	{
		printf("can not open %s file\n",path2);
		exit(1);
	}
	
	int line=0;

	long offset=0;

	fprintf( fp2, "%10d\t", offset );

	while(fgets(str, L1000_CONDITION_LEN , fp1)!= NULL)  //read every line
	{
		offset = ftell(fp1);		
		fprintf( fp2, "%10d\t", offset );   
	}

	fclose(fp1);
	fclose(fp2);
}


long readByteOffsetFile(char path[],int row_num)
{
	FILE *fp;
	char str[20];
	long offset;
	char c[] = "\t";

	if((fp=fopen(path,"r"))==NULL)
	{
		printf("can not open %s file\n",path);
		exit(1);
	}

	fseek(fp,(row_num-1)*11,SEEK_CUR); 
	fgets(str,20,fp);
	offset = atol(strtok(str,c));

	fclose(fp);

	return offset;
}

void getSampleConditions(char path[], long offset, char conditions[])
{

	FILE *fp;
	if((fp=fopen(path,"r"))==NULL)
	{
		printf("can not open %s file\n",path);
		exit(1);
	}

	fseek(fp,offset,SEEK_CUR); 
	fgets(conditions,L1000_CONDITION_LEN,fp);

	fclose(fp);
}

//write classflag vector
void WritetxtClusterResult(int classflag[] ,int len, int cluster, char writepath[])
{
	int i,j,cluster_num=0;
	char conditions[L1000_CONDITION_LEN];
	int *flag,*classes;	
	long cidnum;
	long offset;
	flag = (int *) malloc(len*sizeof(int));
	classes = (int *) malloc(len*sizeof(int));
	memset(flag,0,len*sizeof(int));	
	
	FILE *fp;
	if((fp=fopen(writepath,"w"))==NULL)
	{
		printf("can not open the %s file to write!\n",writepath);
		exit(1);
	}
	
	for( i=0; i < len; i++)
	{
		if(flag[classflag[i]]==0)
			flag[classflag[i]]=(++cluster_num);
		classes[i] = flag[classflag[i]];
	}
	
	for( i=0; i<cluster; i++)
	{
		
		fprintf(fp,"cluster %d :\n",i+1);
		for( j=0; j < len; j++)
		{
			if(classes[j]==(i+1))
			{
				cidnum = readByteOffsetFile("data/data_for_test_cidnum.txt",gsea_result[i].cid);
				offset = readByteOffsetFile("data/prepareForNewDataSet/Samples_RowByteOffset.txt",cidnum);
				getSampleConditions("data/prepareForNewDataSet/Samples_Condition.txt", offset, conditions);
				fprintf(fp,"%s\n",conditions);
			}				
		}
			
	}
	
	fclose(fp);
}