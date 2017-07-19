/* File:     GSEA.h
 * Purpose:  Header file for GSEA.c, which implements a quickGSEA function.
 */
#ifndef _MY_GSEA_H_
#define _MY_GSEA_H_

//#define MAX_GENE 22300
#define MAX_GENE 1000
#define MAX_GENESET 250 
#define Global_ES_SIZE 10000
#define permutation_SIZE 1000
#define L1000_LEN 1000
#define GENE_INFO_LEN 200
#define L1000_CONDITION_LEN 250
#define FileName_LEN 150

struct original_Profile{
	short id;
	double expression;
};

struct Profile_triple{
	short gsUp[MAX_GENESET];
	short gsDown[MAX_GENESET];
	short index[MAX_GENE];
	int cid;
};

struct ES_RESULT{
	float ES;
	int cid;
};

struct GSEA_RESULT{
	float ES;
	float NES;
	double pv;
	int cid;
};

void quiksortINT(int a[],int low,int high);
void quiksort(float a[],int low,int high);
void quiksort_profile(struct original_Profile *profile, int low,int high);
void quiksort_es(struct ES_RESULT *es,int low,int high);
void quiksort_gsea(struct GSEA_RESULT *gsea,int low,int high);
float quickGeneSet(float isgs[], int len,int sig);
void getIndex(short s[],short indexS[],int len);
float ES_GeneSet(short gs[], short indexS[],int len,int sig);

void GetRandomSequence(int total,int size , float* Seq);
void changePosition(float *positions , int size);
void getGlobalES(int len,int sig,float global_ES[]);
void GSEA(short gs[], short indexS[],int len,int sig, float *ES, float *NES, double *pv, float global_ES[]);

struct Profile_triple getTriple(short s[], int len, int sig);
float ES_Profile_original(short s1[], short s2[],int len,int sig);
float ES_Profile(short s1[], short s2[],int len,int sig);
float ES_Profile_triple(struct Profile_triple s1, struct Profile_triple s2 ,int len,int sig);


#endif
