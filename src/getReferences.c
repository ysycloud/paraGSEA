#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <sys/time.h> 
#include <unistd.h> 
#include <getopt.h> 
#include "Tools.h"

#define ERRM "getReferences error:"

char *USAGE =
"\n"
"Usage:"
"  getReferences [options]\n"
"\n"
"  input options: \n"
"    -1 --input1: input file/a separate txt file of gene info. \n"
"    -2 --input2: input file/a separate txt file of inst(profile treatment condition) info. \n"
"\n"
"  output options: \n"
"    -r --reference: input a directory used for outputing referenced files about genesymbols and cids. \n";


void Usage();

int main(int argc,char *argv[])
{	

	FILE *fp;
	char conditionsfile[FileName_LEN];
	char offsetfile[FileName_LEN];
	char genelistfile[FileName_LEN];
	
    // Unset options (value 'UNSET').
	char * const UNSET = "unset";
    char * input1   = UNSET;
	char * input2   = UNSET;
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
			{"input1",              required_argument,        0, '1'},
			{"input2",              required_argument,        0, '2'},
			{"reference",           required_argument,        0, 'r'},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "1:2:r:",
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
				fprintf(stderr, "%s --input1 set more than once\n", ERRM);
				Usage();
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
				fprintf(stderr, "%s --input2 set more than once\n", ERRM);
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
			
		default:
			// Cannot parse. //
			Usage();
			exit(0);
		}		
	}
	
	//check the parameters
	
	
	if((fp=fopen(input1,"r"))==NULL)
	{
		fprintf(stderr, "[ param error : -1 ] can not open gene info file: '%s' !\n", input1);
		exit(0);
	}
	fclose(fp);
	
	if((fp=fopen(input2,"r"))==NULL)
	{
		fprintf(stderr, "[ param error : -2 ] can not open inst info file: '%s' !\n", input2);
		exit(0);
	}
	fclose(fp);
		
	sprintf(genelistfile,"%s/Gene_List.txt",reference);
	sprintf(conditionsfile,"%s/Samples_Condition.txt",reference);
	sprintf(offsetfile,"%s/Samples_RowByteOffset.txt",reference);
	
	//getGeneListFile("GSE92742_Broad_LINCS_gene_info.txt",2,genelistfile);
	//getConditionReference("GSE92742_Broad_LINCS_inst_info.txt",conditionsfile,offsetfile);
	getGeneListFile(input1,2,genelistfile);
	getConditionReference(input2,conditionsfile,offsetfile);
	
	return 0;
}

void Usage() {
	fprintf(stderr, "%s\n", USAGE);
}  /* Usage */