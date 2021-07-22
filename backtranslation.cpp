#include<iostream.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
// Author: Lalitha Viswanathan
// Org: TCS (Tata Consultancy Services)
//reads in codon usage table and amino acid sequence
//computes codon as one which has maximum chances of coding for that amino acid
#define asp 0
#define cys 1  
#define gln 2 
#define glu 3 
#define gly 4  
#define his 5  
#define ile 6  
#define leu 7  
#define lys 8  
#define met 9  
#define phe 10 
#define pro 11 
#define ser 12 
#define thr 13 
#define trp 14 
#define tyr 15 
#define val 16 
#define asn 17
#define ala 18
#define arg 19
#define stop 20
#define MAX 10000

char salpha;
class codonusage
{
	char* DNAseq[MAX],*leastprobDNAseq[MAX];
	char* aminoacidseq[MAX];
	char firststr[MAX][100];
	int aasize[MAX]; int dnasize[MAX];
	int noaminoacids;
	typedef struct aminoacids{
		double values[30],percentage[30];
		char acid;
		int count;
		char codon[30][4];
	}aa;
       aa codonarray[25];
	public:
	codonusage(char* fname);
	void readcodontable(char* fname);
	void setcount(int i){codonarray[i].count++;}
	int getcount(int i){return codonarray[i].count;}
	void backtranslate();
	void printtofile(char*);
	void print();
	~codonusage();
};

int code(char *);
int map(char);
// Print codon usage to file
void codonusage::printtofile(char* fname)
{
	FILE* fp=fopen(fname,"w+");
	for(int i=0;i<noaminoacids;i++)
	{
		fputc('\n',fp);
		fputc('>',fp);
		fprintf(fp,"%s",firststr[i]);
		for(int counter=0;counter<aasize[i];counter++)
			fprintf(fp,"%c",aminoacidseq[i][counter]);
		fputs("\n\t MOST PROBABLE DNA SEQUENCE \t\n",fp);
		for(int j=0;j<dnasize[i];j++)
			fprintf(fp,"%c",DNAseq[i][j]);
		fputs("\n\t LEAST PROBABLE DNA SEQUENCE \t\n",fp);
		fprintf(fp,"%c",leastprobDNAseq[i][j]);
	}
	fclose(fp);
}
codonusage::codonusage(char * fname)
{
	FILE* filepointer;
	if(!(filepointer=fopen(fname,"r")))
	{
		printf("%s does not exist\n",fname);
		exit(1);
	}
	char ch,str[500];
	noaminoacids=-1;
	int counter=0;
	for(int i=0;i<MAX;i++)
	aminoacidseq[i]=(char *)malloc(sizeof(char)*1);
	while(filepointer) 
	{
		ch=fgetc(filepointer);
		if(feof(filepointer))
			break;

		if(ch=='>')
		{
				noaminoacids++;
				fgets(firststr[noaminoacids],100,filepointer);
				aasize[noaminoacids-1]=counter;
				counter=0;
		}
		else if(isalpha(ch))
		{
				aminoacidseq[noaminoacids][counter]=ch;
				counter++;
				aminoacidseq[noaminoacids]=(char *)realloc(aminoacidseq[noaminoacids],sizeof(char)*(strlen(aminoacidseq[noaminoacids])+1));
		}
	}
fclose(filepointer);
for(int aacounter=0;aacounter<noaminoacids;aacounter++)
{
		DNAseq[aacounter]=(char *)malloc(sizeof(char)*(aasize[aacounter]*3));
}
}
//	Constructor
	codonusage::~codonusage()
	{
	}

	// Read codon table 
void codonusage::readcodontable(char* fname)
{
	FILE* fp=fopen(fname,"r");
	char* s, codn[4],str[100];
	double num, percent;
	for(int codonarrcounter=0;codonarrcounter<23;codonarrcounter++)
			codonarray[codonarrcounter].count=0;
	while(fp)
	{
		char* s=fgets(str,100,fp);
		if(feof(fp))
				break;
		sscanf(s,"\n%s %lf (%lf)",codn,&num,&percent);
		//num is number of times it(codn) occurs in organism
		int index=code(codn);
		codonarray[index].values[codonarray[index].count]=num;
		//values is number of time it(codn) occurs in organism
		//percentage is the percentage occurance of the codon in organism
		codonarray[index].percentage[codonarray[index].count]=percent;
		//salpha is the alphabet that repsents aminon acid
		codonarray[index].acid=salpha;
		//copy the codn sequence into the codon array of amino acid
		strcpy(codonarray[index].codon[codonarray[index].count], codn);
		//increment the count of number of codons coding for amino acid
		codonarray[index].count++;
	}
	fclose(fp);
	printf("before sorting the values array and percentage array\n");
	for(int i=0;i<23;i++)
	{
		double percentage,numb;
		char tempcodon[4],salpatemp;
		//after this code is executed, for each amino acids, first set of values will be the one corresponding to the most frequently used codon
		//and last entry wud be the one corresponding to the least used codon.
		for(int j=0;j<codonarray[i].count;j++)
			for(int k=j+1;k<codonarray[i].count;k++)
				if(codonarray[i].values[k]>=codonarray[i].values[j])
				{
					numb=codonarray[i].values[k];
					percentage=codonarray[i].percentage[k];
					strcpy(tempcodon,codonarray[i].codon[k]);

					codonarray[i].values[k]=codonarray[i].values[j];
					codonarray[i].percentage[k]=codonarray[i].values[j];
					strcpy(codonarray[i].codon[k],codonarray[i].codon[j]);
					
					codonarray[i].values[j]=numb;
					codonarray[i].percentage[j]=percentage;
					strcpy(codonarray[i].codon[j],tempcodon);

				}
	}
}
// print codon usage
void codonusage::print()
{
		for(int counter=0;counter<20;counter++)
		{
				printf("Amino acid= %c codons =\n",codonarray[counter].acid);
				for(int innercounter=0;innercounter<codonarray[counter].count;innercounter++)
						printf("\t codon %s percentage %lf values %lf\n",codonarray[counter].codon[innercounter],codonarray[counter].percentage[innercounter],codonarray[counter].values[innercounter]);
		}
}
// calculate codon usage
void codonusage::backtranslate()
{
	char cod[4],leastprobcod[4];
	for(int counter=0;counter<=noaminoacids;counter++)
	{
		int n=0;
		DNAseq[counter]=(char*)malloc(sizeof(char)*(aasize[counter]*3));
		leastprobDNAseq[counter]=(char *)malloc(sizeof(char)*(aasize[counter]*3));
		dnasize[counter]=aasize[counter]*3;
		for(int innercounter=0;innercounter<aasize[counter];innercounter++)
		{
			for(int aminoacidcounter=0;aminoacidcounter<22;aminoacidcounter++)
				if(codonarray[aminoacidcounter].acid==aminoacidseq[counter][innercounter])
				{
					strcpy(cod,codonarray[aminoacidcounter].codon[0]);
					strcpy(leastprobcod,codonarray[aminoacidcounter].codon[codonarray[aminoacidcounter].count-1]);
					DNAseq[counter][n]=cod[0];
					leastprobDNAseq[counter][n]=leastprobcod[0];
				    n++;	
					DNAseq[counter][n]=cod[1];	
					leastprobDNAseq[counter][n]=leastprobcod[1];
					n++;
					DNAseq[counter][n]=cod[2];
					leastprobDNAseq[counter][n]=leastprobcod[2];
					n++;
					break;	
				}	
		}
		dnasize[counter]=n;
	}
}
int map(char aminoacid)
{
	if(aminoacid=='A')
		return 0;
	else if(aminoacid=='T')
		return 1;
	else if(aminoacid=='U')
		return 1;
	else if(aminoacid=='G')
		return 2;
	else
		return 3;
}
// Amino Acid code
int code(char* codon)
{
	int number=(map(codon[2])*1) + (map(codon[1])*4)  + (map(codon[0])*16);
	switch(number)
	{
		case 0:
			salpha='K';
			return lys;
		case 1:
			salpha='N';
			return asn;
		case 2:
			salpha='K';
			return lys;
		case 3:
			salpha='N';
			return asn;
		case 4:
			salpha='I';
			return ile;
		case 5:
			salpha='I';
			return ile;
		case 6:
			salpha='M';
			return met;
		case 7:
			salpha='I';
			return ile;
		case 8:
			salpha='R';
			return arg;
		case 9:
			salpha='S';
			return ser;
		case 10:
			salpha='R';
			return arg;
		case 11:
			salpha='S';
			return ser;
		case 12:
			salpha='T';
			return thr;
		case 13:
			salpha='T';
			return thr;
		case 14: 
			salpha='T';
			return thr;
		case 15:
			salpha='T';
			return thr;
		case 16:
			salpha='X';
			return stop;
		case 17:
			salpha='Y';
			return tyr;
		case 18:
			salpha='X';
			return stop;
		case 19:
			salpha='Y';
			return tyr;
		case 20:
			salpha='L';
			return leu;
		case 21:
			salpha='F';
			return phe;
		case 22:
			salpha='L';
			return leu;
		case 23:
			salpha='F';
			return phe;
		case 24:
			salpha='X';
			return stop;
		case 25:
			salpha='C';
			return cys;
		case 26:
			salpha='W';
			return trp;
		case 27:
			salpha='C';
			return cys;
		case 28:
			salpha='S';
			return ser;
		case 29:
			salpha='S';
			return ser;
		case 30:
			salpha='S';
			return ser;
		case 31:
			salpha='S';
			return ser;
		case 32:
			salpha='E';
			return glu;
		case 33:
			salpha='D';
			return asp;
		case 34:
			salpha='E';
			return glu;
		case 35:
			salpha='D';
			return asp;
		case 36:
			salpha='V';
			return val;
		case 37:
			salpha='V';
			return val;
		case 38:
			salpha='V';
			return val;
		case 39:
			salpha='V';
			return val;
		case 40:
			salpha='G';
			return gly;
		case 41:
			salpha='G';
			return gly;
		case 42:
			salpha='G';
			return gly;
		case 43:
			salpha='G';
			return gly;
		case 44:
			salpha='A';
			return ala;
		case 45:
			salpha='A';
			return ala;
		case 46:
			salpha='A';
			return ala;
		case 47:
			salpha='A';
			return ala;
		case 48:
			salpha='Q';
			return gln;
		case 49:
			salpha='H';
			return his;
		case 50:
			salpha='Q';
			return gln;
		case 51:
			salpha='H';
			return his;
		case 52:
			salpha='L';
			return leu;
		case 53:
			salpha='L';
			return leu;
		case 54:
			salpha='L';
			return leu;
		case 55:
			salpha='L';
			return leu;
		case 56:
			salpha='R';
			return arg;
		case 57:
			salpha='R';
			return arg;
		case 58:
			salpha='R';
			return arg;
		case 59:
			salpha='R';
			return arg;
		case 60:
			salpha='P';
			return pro;
		case 61:
			salpha='P';
			return pro;
		case 62:
			salpha='P';
			return pro;
		case 63:
			salpha='P';
			return pro;
	}
}
// Main Program
main(int argc, char* argv[])
{
	if(argc!=3)
	{
		printf("\n usage <executable> <backtrans table file> <file containing amino acid sequence(s)>\n");
		exit(1);
	}
	else
	{
		char* btrans=argv[1];
		char* aaseq=argv[2];
		printf("btrans is %s aaseq is %s\n",btrans,aaseq);
		codonusage cu(aaseq);
		printf("after reading file\n");
		//reads the codonusage table for the organism 
		cu.readcodontable(btrans);
		printf("after reading backtransfile\n");
		cu.backtranslate();
		cu.printtofile("results.txt");
		printf("results are in results.txt\n");
	}
}
