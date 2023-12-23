#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define PWMLEN 50 // max length of motives
#define OLIGNUM 16

double pfm[PWMLEN][OLIGNUM];
double pwm[PWMLEN][OLIGNUM];

// delete symbol 'c' from input string
void DelChar(char *str,char c)
{
	int i, lens, size;

	size=0;
	lens=strlen(str);
	for(i=0;i<lens;i++)
	{
		if(str[i]!=c)str[size++]=str[i];
	}
	str[size]='\0';
}
int StrNStr(char *str,char c, int n)
{
	int i, len=strlen(str);
	int k=1;
	for(i=0;i<len;i++)
	{
		if(str[i]==c)
		{
			if(k==n)return i;
			k++;
		}
	}
	return -1;
}
int UnderStol(char* str, int nstol, char* ret, size_t size, char sep)
{
	memset(ret, 0, size);
	int p1, p2, len;
	if (nstol == 0)
	{
		p2 = StrNStr(str, sep, 1);
		if (p2 == -1)p2 = strlen(str);
		if (p2 == 0) return -1;
		strncpy(ret, str, p2);
		ret[p2] = '\0';
		return 1;
	}
	else
	{
		p1 = StrNStr(str, sep, nstol);
		p2 = StrNStr(str, sep, nstol + 1);
		if (p2 == -1)
		{
			p2 = strlen(str);
			if (p2 == 0) return -1;
		}
		if (p1 == -1 || p2 == -1) return -1;
		len = p2 - p1 - 1;
		strncpy(ret, &str[p1 + 1], len);
		ret[len] = '\0';
		return 1;
	}
}
int main(int argc, char *argv[])
{	
	char file_pfm[300], file_pwmt[300], file_log[]="pfm_to_pwm_mat.len";//file_pwmb[80],
	char d[500],head[500], s[500];	
    int i, j, len=0; 
   	if(argc!=3)
    {
        printf ("Syntax: 1file in pfm 2file out pwm txt");
        return -1;
    }
	double nseq;//=atoi(argv[3]);
	int shift_col;//=atoi(argv[4]);
	//3 int virtual_sample_size 4int shift_column(homer 0 cis-db 1) 5char Rr=row(Jaspar) Cc=column
	int mtype;
	//if(strchr("Rr",argv[5][0])!=NULL)mtype=1;//jaspar = 4 ryada L kolonok
	//else mtype=0;//Homer L ryadov  4 kolonki

	strcpy(file_pfm,argv[1]);
	strcpy(file_pwmt,argv[2]);	
	//strcpy(file_pwmb, argv[3]);
	//int nmat;
	//if (argc == 4)nmat = 0;
	//else nmat = atoi(argv[4]);
	FILE *in_pfm;
	if((in_pfm=fopen(file_pfm,"rt"))==NULL)
	{
		printf("Input file %s can't be opened!\n",file_pfm);
		return -1;
	}
	int alfabet=0;
	fgets(head,sizeof(head),in_pfm);			 
	fgets(head,sizeof(head),in_pfm);			 
	if(strchr("ATGC",head[0])!=NULL)
	{// jaspar
		nseq=1;
		mtype=1;
		shift_col=1;
		alfabet++;
		while(fgets(head,sizeof(head),in_pfm)!=NULL)
		{
			if(strchr("ATGC",head[0])!=NULL)alfabet++;
		}
		if(alfabet!=4 && alfabet!=16)
		{
			printf("Reading error in Jaspar matrix file %s !\n",file_pfm);
			return -1;
		}
	}
	else
	{//homer, cis-bp
		nseq=1000000000;
		mtype=0;
		//shift_col=0;
	}	
	if(mtype==0)
	{
		DelChar(head,'\n');
		DelChar(head,'\r');
		int headlen= strlen(head);
		if(head[headlen-1]=='\t')
		{
			head[headlen-1]='\0';
			headlen--;
		}
		int counttab=0;
		for(i=0;i<headlen;i++)
		{
			if(head[i]=='\t')counttab++;
		}
		if(counttab==4 || counttab==16)
		{
			shift_col=1;
			alfabet=counttab;
		}
		else 
		{
			if(counttab==3 || counttab==15)
			{
				shift_col=0;
				alfabet=counttab+1;
			}
			else
			{
				printf("Reading errorin Cisbp/Homer matrix file %s !\n",file_pfm);
				return -1;
			}
		}
	}
	rewind(in_pfm);
	int test;
	char sep = '\t';
	FILE *out_log;
	if((out_log=fopen(file_log,"at"))==NULL)
	{
		printf("Input file %s can't be opened!\n",file_log);
		return -1;
	}
	int olen=0;
	double olen_perf=0;
	fgets(head,sizeof(head),in_pfm);			 
	if(mtype==0)
	{		
		for(i=0;i<PWMLEN;i++)
		{
			if(fgets(d,sizeof(d),in_pfm)!=NULL)
			{
				DelChar(d,'\n');
				DelChar(d,'\r');
				char c=d[0];
				if(isdigit(c) || (strchr("-ATGC",c)!=0))
				{
					olen++;
					//	 printf("%s",d[i]);						
					for(j=0;j<alfabet;j++)
					{			
						/*if (NthColumn(d, s, j + shift_col) == -1)
						{
							break;			
						}*/	
						test = UnderStol(d, j + shift_col, s, sizeof(s), sep);
						if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
						double score=atof(s);
						pfm[i][j]=nseq*score;										
					}					
				}
				else break;
			}
			else break;
		}	
	}
	else
	{
		fgets(d,sizeof(d),in_pfm);
		int dlen=strlen(d);
		olen=0;
		for(i=1;i<dlen;i++)
		{
			if(d[i-1]=='\t' && isdigit(d[i]))olen++;
		}
		rewind(in_pfm);
		fgets(head,sizeof(head),in_pfm);
		for(i=0;i<alfabet;i++)
		{
			if(fgets(d,sizeof(d),in_pfm)!=NULL)
			{
				DelChar(d,'\n');
				char c=d[0];
				if(isdigit(c) || (strchr("-ATGC",c)!=0))
				{				
					//	 printf("%s",d[i]);	
					for(j=0;j<olen;j++)
					{			
						/*if (NthColumn(d, s, j + shift_col) == -1)
						{
							break;			
						}*/
						test = UnderStol(d, j + shift_col, s, sizeof(s), sep);
						if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
						double score=atof(s);
						pfm[j][i]=nseq*score;					
					}
				}
				else break;
			}
			else break;
		}
	}
	fclose(in_pfm);
	int alfabet_m1=alfabet-1;
	for(j=0;j<olen;j++)
	{		
		int zero_letter=0;
		for(i=0;i<alfabet;i++)
		{
			if(pfm[j][i]==0)zero_letter++;
		}
		if(zero_letter==alfabet_m1)olen_perf++;
	}
	fprintf(out_log,"%s\t%d\t%.f\n",argv[1],olen,olen_perf);
	fclose(out_log);
	//printf("matrix length %d\n",olen);
/*	for(i=0;i<4;i++)
	{
		for(j=0;j<olen;j++)
		{
			printf("%d\t",pfm[j][i]);
		}
		printf("\n");
	}*/	
/*	int nmat = 0;
	int mlen=strlen(argv[1]);
	for(i=0;i<mlen;i++)
	{
		if(isdigit(argv[1][i]))
		{
			nmat=atoi(&argv[1][i]);
			nmat--;
			break;
		}
	}*/
	FILE* out_pwmt;// , * out_pwmb;
	if((out_pwmt=fopen(file_pwmt,"wt"))==NULL)
	{
		printf("Ouput file %s can't be opened!\n",file_pwmt);
		return -1;
	}
/*	if ((out_pwmb = fopen(file_pwmb, "wb")) == NULL)
	{
		printf("Ouput file %s can't be opened!\n", file_pwmb);
		return -1;
	}*/
	/*char file_pwm2[] = "pwm.h";
	FILE *out_pwm2;
	if((out_pwm2=fopen(file_pwm2,"at"))==NULL)
	{
		printf("Ouput file %s can't be opened!\n",file_pwm2);
		return -1;
	}
	fprintf(out_pwm2,"//%s\n",file_pwmt);
	fprintf(out_pwm2,"case(%d):\t{\n",nmat+1);
	fprintf(out_pwm2,"\tmatrix->len=%d;\n",olen);
	fprintf(out_pwm2,"\tmatrix->mem_in(%d);\n",olen);

	fwrite(&olen, sizeof(int), 1, out_pwmb);
	if(*head!='>')fprintf(out_pwmt,"> ");*/
	fprintf(out_pwmt,"%s",head);
	{//logodds score
		double pse=0.25;
		double pse4=1;
		double vych=log10(pse);		
		int olen1=olen-1;
		int alfabet1=alfabet-1;
		for(i=0;i<olen;i++)
		{
			//fprintf(out_pwm2,"\tmatrix->init_wei(%d, ",i);			
			double sum=0;
			for(j=0;j<alfabet;j++)sum+=pfm[i][j];			
			for(j=0;j<alfabet;j++)
			{				
				double ves=(double)(pfm[i][j]+pse)/(sum+pse4);
				ves=log10(ves)-vych;
				pwm[i][j] = ves;
				fprintf(out_pwmt,"%.15f",ves);
				//fprintf(out_pwm2,"%.15f",ves);
				if(j==alfabet1)
				{
					fprintf(out_pwmt,"\n");					
				}
				else 
				{
					fprintf(out_pwmt,"\t");		
					//fprintf(out_pwm2,", ");										
				}
			}			
			//fprintf(out_pwm2,");\n");		
		}
		//fwrite(pwm, sizeof(double), 4*olen, out_pwmb);
		//fprintf(out_pwm2,"\tbreak;\t}\n");
		/*for (i = 0; i<olen; i++)
		{
			fprintf(out_pwm2,"\tmatrix->init_fre(%d, ",i);			
			int alfabet1=alfabet-1;
			for(j=0;j<alfabet;j++)
			{				
				fprintf(out_pwm2,"%.15f",(double)pfm[i][j]/nseq);
				if(j!=alfabet1)
				{
					fprintf(out_pwm2,", ");										
				}
			}			
			fprintf(out_pwm2,");\n");		
		}
		fprintf(out_pwm2,"\tbreak;\t}\n"); */
	}
	fclose(out_pwmt);
//	fclose(out_pwmb);
	//fclose(out_pwm2);
/*	FILE* out_len;
	char file_len[80];
	memset(file_len,'\0',sizeof(file_len));
	for(i=0;;i++)
	{
		if(file_pfm[i]=='\0')break;
		file_len[i]=file_pfm[i];
		if(file_pfm[i]=='.')break;
	}
	strcat(file_len,"len");
	if((out_len=fopen(file_len,"wt"))==NULL)
	{
		printf("Ouput file %s can't be opened!\n",file_len);
		exit(1);
	}
	fprintf(out_len,"%d",olen);
	fclose(out_len);
	printf("%d",olen);*/
	return 1;
}
