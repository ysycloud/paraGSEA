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
		printf("file error!\n"); 
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
		printf("file error!\n"); 
		return -1; 
	} 
	fgets(StrLine,256,fp);
	line = BeginLine;
	//不能因为前面取LineLength是strlen就以为有\0,这里+1，每行的后面本来也没有\0
	fseek(fp,BeginLine*(LineLength),SEEK_CUR); 
	while (!feof(fp)) 
    { 
		fgets(StrLine,LineLength+1,fp);  //read one line
		col = 0;
		if( line<profilenum && col< genelen )   //防止有些空行或多余后列导致超界
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

//write part of ES_Matrix to txt file，can be open by excel...
void WritetxtResult(int sourceBegin ,int sourceEnd, int matlen, char writepath[], float **ES_Matrix)
{
	int i,j;
	int local_sourcelen = sourceEnd - sourceBegin;
	
	FILE *fp;
	if((fp=fopen(writepath,"w"))==NULL)
	{
		printf("can not open the file to write!");
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
	char c[] = " ", *p;
	int tmp,i;
	*count = 0;
	
	//initial flag vector
	memset(existflag, 0, MAX_GENE * sizeof(short));
		
	//get gstmp and gs count,remove the repeat elements
	gstmp[(*count)++] = atoi(strtok(gsStr,c));
	existflag[gstmp[(*count)-1]] = 1;
	p = strtok(NULL,c);
	while(p)
	{
		tmp = atoi(p);
		if(existflag[tmp]==0){	//this gene not input
			gstmp[(*count)++] = tmp;
			existflag[gstmp[(*count)-1]] = 1;
		}
		p = strtok(NULL,c);
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
		printf("file error!\n"); 
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
		printf("file error!\n"); 
		return -1; 
	} 
	fgets(StrLine,256,fp);
	line = BeginLine;
	//不能因为前面取LineLength是strlen就以为有\0,这里+1，每行的后面本来也没有\0
	fseek(fp,BeginLine*(LineLength),SEEK_CUR); 
	while (!feof(fp)) 
    { 
		fgets(StrLine,LineLength+1,fp);  //read one line
		col = 0;
		if( line<profilenum1 && col< profilenum2 )   //防止有些空行或多余后列导致超界
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

//write classflag vector
void WritetxtClusterResult(int classflag[] ,int len, char writepath[])
{
	int i;
	
	FILE *fp;
	if((fp=fopen(writepath,"w"))==NULL)
	{
		printf("can not open the file to write!");
		exit(1);
	}
	
	for( i=0; i < len; i++)
		fprintf(fp,"%d\n",classflag[i]);

	fclose(fp);
}
