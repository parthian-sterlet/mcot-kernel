#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 12000 // max length of peak in input fasta
#define MATLEN 50 //max matrix length
#define SPACLEN 150 //max upper bound of spacer length
#define OLIGNUM 4// di 16 mono 4
#define NMAT_HS_CORE 403 // total count of matrices in hocomoco human core collection +1(anchor)
#define NMAT_MM_CORE 359 // total count of matrices in hocomoco mouse core collection +1(anchor)
#define NMAT_HS_FULL 772 // total count of matrices in hocomoco human full collection +1(anchor)
#define NMAT_MM_FULL 532 // total count of matrices in hocomoco mouse full collection +1(anchor)
#define NMAT_DAPSEQ 529 // total count of matrices in hocomoco mouse full collection +1(anchor)
#define NUM_THR 5 //4islo porogov
#define NUM_PVAL 30000 // max 4islo porogov v tablice Touzet
#define MOT_NAME_LEN 40 //max length of motif name

//return n-th occurrence of a certain symbol c in a string str
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
// menyaet registr stroki
char *TransStr(char *d)
{
   int i, c, lens;
   lens=strlen(d);
   for(i=0;i<lens;i++)
   {
	c=int(d[i]);
	if(c<97) d[i]=char(c+32);
	//else break;
   }
   return(d);
}
//udalenie simvola iz stroki
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
// ras4et 4asot oligonukleotidov po stroke (zdes' - nukleotidov)
void GetSostPro(char *d, int word, int *sost)
{
   int i, j, k, i_sost, let;   
   char letter[]="acgt";
   int ten[6]={1,4,16,64,256,1024};
   int lens=strlen(d);
   int size=1;
   for(k=0;k<word;k++)size*=4;
   for(i=0;i<size;i++)sost[i]=0;
   for(i=0;i<lens-word+1;i++)
   {
   	i_sost=0;
   	for(j=word-1;j>=0;j--)
       {
         for(k=0;k<4;k++)
         {
      		if(d[i+j]==letter[k]){let=k;break;}
         }
         i_sost+=ten[word-1-j]*let;
      }
      sost[i]=i_sost;
   }
}
//komplementaciya stroki
int ComplStr(char *d)
{
	char *d1;
	int i, len;
	len=strlen(d);
	d1=new char[len+1];
	if(d1==NULL) return 0;
	strcpy(d1,d);
//	memset(d,0,sizeof(d));
	for(i=0;i<len;i++)
	{
		switch(d1[len-i-1])
		{
			case 'a':{d[i]='t';break;}
			case 't':{d[i]='a';break;}
			case 'c':{d[i]='g';break;}
			case 'g':{d[i]='c';break;}
			case 'A':{d[i]='T';break;}
			case 'T':{d[i]='A';break;}
			case 'C':{d[i]='G';break;}
			case 'G':{d[i]='C';break;}
			case 'N':{d[i]='N';break;}
			case 'n':{d[i]='n';break;}
			default: d[i]='n';
		}
	}
	delete [] d1;
	return 1;
}
void Mix(int *a, int *b)
{
  int buf=*a;
  *a=*b;
  *b=buf;
}
void Mix(char *a, char *b)
{
  char buf=*a;
  *a=*b;
  *b=buf;
}
void BigMix1(char *d)//me6alka
{
    int r;
    int len=strlen(d);
 for(r=0;r<len-1;r++) Mix(&d[r], &d[1+r+(rand()%(len-1-r))]);
}
void BigMix1(int *d1, int len) // pereme6ivanie stroki
{
    int r, s;    
 for(r=0;r<len;r++)
 {
	 s = rand() % len;
	 if(s!=r)Mix(&d1[r], &d1[s]);	 
 }
}
#include "fasta_to_plain.h" //input = peaks (fasta), output = plain format OR lengths of peaks
#include "pwm_score_distr_granul.h" // touzet algorithm for PWM (PWM threshold <---> PWM p-value), PMID: 18072973 
#include "pwm_iz_pwm_thr_dist.h" // full list of touzet thresholds for PWM -> FP-rates po genomu
#include "select_thresholds_from_pvalues.h" //fp-rates po genomu -> vybor pyati porogov po pyati fixed FP rates

//thresholds
#include "hocomoco_thr_hs_core.h" //porogi po human hocomoco core collection of motifs = partnery
#include "hocomoco_thr_hs_full.h" //porogi po human hocomoco full collection of motifs = partnery
#include "hocomoco_thr_mm_core.h" //porogi po mouse hocomoco core collection of motifs = partnery
#include "hocomoco_thr_mm_full.h" //porogi po mouse hocomoco full collection of motifs = partnery
#include "dapseq_thr.h" //porogi po arabidopsis dapseq collection of motifs = partnery

// position frequency mattrix (PFM), position weight matrix (PWM)
struct matrices {
	int len;
	double min;
	double raz;
	double **wei;
	double **fre;
	void init_fre(int i, double a, double b, double c, double d);
	void init_wei(int i, double a, double b, double c, double d);
	void init_dapseq(int num);
	void init_hs_core(int num);
	void init_hs_full(int num);
	void init_mm_core(int num);
	void init_mm_full(int num);
	void get_copy(matrices *a);
	int mem_in(int len);
	void mem_out(int len);
	void norm(void);	
};

int matrices::mem_in(int length)
{
	int i;
	len=length;
	wei=new double * [len];
	if(wei==NULL) return -1;
	for(i=0;i<len;i++)
	{
		wei[i]=new double[OLIGNUM];
		if(wei[i]==NULL) return -1;	
	}
	fre=new double * [len];
	if(fre==NULL) return -1;
	for(i=0;i<len;i++)
	{
		fre[i]=new double[OLIGNUM];
		if(fre[i]==NULL) return -1;	
	}
	return 1;
}
void matrices::mem_out(int len)
{
	int i;
	for(i=0;i<len;i++) delete [] wei[i];
	delete [] wei;
	for(i=0;i<len;i++) delete [] fre[i];
	delete [] fre;
}
void matrices::norm(void)
{
	int i,j;
	min=raz=0;
	for(i=0;i<len;i++)
	{
		double pwmmin=100;
		double pwmmax=-100;
		for(j=0;j<OLIGNUM;j++)
		{							 
			double w=wei[i][j];
			if(w<pwmmin)pwmmin=w;
			if(w>pwmmax)pwmmax=w;
		}
		raz+=pwmmax;
		min+=pwmmin;
	}	
	raz-=min;
}

void matrices::init_wei(int i, double a, double b, double c, double d)
{
	wei[i][0]=a;
	wei[i][1]=b;
	wei[i][2]=c;
	wei[i][3]=d;
}
void matrices::init_fre(int i, double a, double b, double c, double d)
{
	fre[i][0]=a;
	fre[i][1]=b;
	fre[i][2]=c;
	fre[i][3]=d;
}
void matrices::get_copy(matrices *a)
{
	a->mem_in(len);
	a->len=len;
	a->min=min;
	a->raz=raz;
	int i,j;
	for(i=0;i<len;i++)
	{
		for(j=0;j<OLIGNUM;j++)
		{
			a->fre[i][j]=fre[i][j];
			a->wei[i][j]=wei[i][j];
		}
	}
}
#include "pfm_to_pwm.h" //conversion PFM -> PWM

#include "hocomoco_pwm_hs_core.h" //PWMs, human hocomoco core collection
#include "hocomoco_pwm_hs_full.h" //PWMs, human hocomoco full collection
#include "hocomoco_pwm_mm_core.h" //PWMs, mouse hocomoco core collection
#include "hocomoco_pwm_mm_full.h" //PWMs, mouse hocomoco full collection
#include "dapseq_pwm.h" //PWMs, arabidopsis dapseq collection

//profil' raspoznannuh saytov po pyati porogam i odin slitiy profil' (the most permissive threshold)
struct profile {
	int mot;//motif num
	int nam;//nomer poroga
	int nseq;
	int *nsit;//4islo saytov
	int nsit_all;
	int nseq_rec;
	int **sta;//nthr nsit
	char **cep;
	double **sco;	
	int mem_in_sta(void);
	int mem_in_cep(void);
	int mem_in_sco(void);
	int mem_in_nsit(void);
	void mem_out_sta(void);
	void mem_out_cep(void);
	void mem_out_sco(void);		
	void mem_out_nsit(void);	
	int get_copy_rand(profile *a, int height);
	int get_copy_real(profile *a);
	int fprintf_pro(char *mot_db, double thr, char *mode);//mot_db real/rand
	void count_sites(void);
	int test(void);
};
int profile::mem_in_sta(void)
{
	int i;
	sta=new int*[nseq];
	if(sta==NULL) return -1;
	for(i=0;i<nseq;i++)
	{
		if(nsit[i]>0)
		{
			sta[i]=new int[nsit[i]];
			if(sta[i]==NULL) return -1;
		}
	}
	return 1;
}
int profile::mem_in_cep(void)
{
	int i;
	cep=new char*[nseq];
	if(cep==NULL) return -1;
	for(i=0;i<nseq;i++)
	{
		if(nsit[i]>0)
		{
			cep[i]=new char[nsit[i]];	
			if(cep[i]==NULL) return -1;
		}
	}
	return 1;
}
int profile::mem_in_sco(void)
{
	int i;
	sco=new double*[nseq];
	if(sco==NULL) return -1;
	for(i=0;i<nseq;i++)
	{
		if(nsit[i]>0)
		{
			sco[i]=new double[nsit[i]];
			if(sco[i]==NULL) return -1;
		}
	}
	return 1;
}
void profile::mem_out_sta(void)
{
	int i;
	for(i=0;i<nseq;i++)if(nsit[i]>0)delete [] sta[i];		
	delete [] sta;
}
void profile::mem_out_cep(void)
{
	int i;
	for(i=0;i<nseq;i++)if(nsit[i]>0)delete [] cep[i];
	delete [] cep;
}
void profile::mem_out_sco(void)
{
	int i;
	for(i=0;i<nseq;i++)if(nsit[i]>0)delete [] sco[i];
	delete [] sco;
}
int profile::mem_in_nsit(void)
{	
	nsit=new int[nseq];
	if(nsit==NULL) return -1;
	int i;
	for(i=0;i<nseq;i++)nsit[i]=0;
	return 1;
}
void profile::mem_out_nsit(void)
{	
	delete [] nsit;
}
int profile::get_copy_rand(profile *a, int height)
{
	int i,j,h,z, ini;
	a->mot=mot;
	a->nam=nam;
	a->nseq=height*nseq;
	a->nsit_all=height*nsit_all;
	a->nseq_rec=height*nseq_rec;
	z=0;
	for(i=0;i<nseq;i++)
	{
		int nsi=nsit[i];
		for(h=0;h<height;h++)
		{
			a->nsit[z++]=nsi;
		}
	}
	ini=a->mem_in_sta();
	if(ini==-1){puts("Not enough memory...\n");return -1;}
	ini=a->mem_in_cep();
	if(ini==-1){puts("Not enough memory...\n");return -1;}
	z=0;
	for(i=0;i<nseq;i++)
	{		
		for(h=0;h<height;h++)
		{		
			for(j=0;j<a->nsit[z];j++)
			{
				a->sta[z][j]=sta[i][j];			
				a->cep[z][j]=cep[i][j];			
			}
			z++;
		}		
	}
	return 1;
}
int profile::get_copy_real(profile *a)
{
	int i,j,z=0, ini;
	a->mot=mot;
	a->nam=nam;
	a->nseq=nseq;
	a->nsit_all=nsit_all;
	a->nseq_rec=nseq_rec;
	for(i=0;i<nseq;i++)
	{
		a->nsit[i]=nsit[i];
	}
	ini=a->mem_in_sta();
	if(ini==-1){puts("Not enough memory...\n");return -1;}
	ini=a->mem_in_cep();
	if(ini==-1){puts("Not enough memory...\n");return -1;}
	ini=a->mem_in_sco();
	if(ini==-1){puts("Not enough memory...\n");return -1;}
	for(i=0;i<nseq;i++)
	{		
		for(j=0;j<a->nsit[i];j++)
		{
			a->sta[i][j]=sta[i][j];			
			a->cep[i][j]=cep[i][j];			
			a->sco[i][j]=sco[i][j];			
		}	
	}
	return 1;
}
int profile::fprintf_pro(char *mot_db, double thr,char *mode)
{
	int print_sco;
	if(strncmp(mode,"real",4)==0)print_sco=1;
	else print_sco=0;
	int i,j;
	char fileo[80];
	FILE *out;
	memset(fileo,'\0',sizeof(fileo));
	strcpy(fileo,mode);
	strcat(fileo,"_");
	strcat(fileo,mot_db);//real or random	
	char buf[10];
	sprintf(buf,"%d",mot);
	strcat(fileo,buf);//hocomoco or dapseq		
	strcat(fileo,"_thr");
	memset(buf,'\0',sizeof(buf));
	sprintf(buf,"%d",nam);
	strcat(fileo,buf);//nomer poroga	
	if((out=fopen(fileo,"wt"))==NULL)
	{
 		printf("Input file %s can't be opened!\n",fileo);
		return -1;
	}	
	for(i=0;i<nseq;i++)
	{		
		fprintf(out,">Seq %d\tThr %f\tNsites %d\n",i+1,thr,nsit[i]);
		for(j=0;j<nsit[i];j++)
		{
			fprintf(out,"%d\t",sta[i][j]);			
			if(print_sco==1)fprintf(out,"%f",sco[i][j]);			
			fprintf(out,"\t");
			fprintf(out,"%c\n",cep[i][j]);						;			
		}	
	}
	fclose(out);
	return 1;
}
void profile::count_sites(void)
{
	nsit_all=nseq_rec=0;
	int i;
	for(i=0;i<nseq;i++)
	{
		int ns=nsit[i];
		if(ns>0)
		{
			nsit_all+=ns;
			nseq_rec++;
		}
	}
}
int profile::test(void)
{
	int i, j;
	for(i=0;i<nseq;i++)
	{
		for(j=0;j<nsit[i];j++)
		{
			if(sta[i][j]>5*SEQLEN || sta[i][j]<0)return -1;
		}
	}
	return 1;
}
#include "pwm_rec.h" //raspoznavanie matricey

//dlya pods4eta zna4imostey CE, full_overlap, partial_overlap, spacer (overlap = partial_overlap OR full_overlap)
struct count {
	int two_sites;
	int any;
	int partial;
	int full;
	int over;
	int spacer;
	void ini(void);
};
void count::ini(void)
{
	two_sites=any=partial=full=over=spacer=0;
}

//dlya vyvoda rezultata po zna4imostyam CE b asymmetrii CE
struct result {
	struct count observed;
	struct count expected;
	double fold_a;//any
	double pvalue_a;
	double fold_p;//partial
	double pvalue_p;
	double fold_f;//full
	double pvalue_f;
	double fold_s;//spacer
	double pvalue_s;
	double fold_o;//overlap
	double pvalue_o;
	void ini(void);
} ;
void result::ini(void)
{
	fold_a=pvalue_a=1;
	fold_f=pvalue_f=1;
	fold_p=pvalue_p=1;
	fold_o=pvalue_o=1;
	fold_s=pvalue_s=1;	
	observed.ini();
	expected.ini();
}

// dlya vyvoda histogram of CE distribution as fanction of mutual orientation and location of anchor/partner motifs
struct combi {	
//	int n_partial;
//	int n_full;
	//int n_tot;
	double freq[4][MATLEN+SPACLEN];	
//	void ini(int len_a, int len_p, int len_sp);
//	void mem_out(void);
	int fprintf_all(char *file, int mot, char *motif, int len_a, int len_p, int len_sp);
};
int combi::fprintf_all(char *file, int mot, char *motif, int len_a, int len_p, int len_sp)
{
	char head[4][10];
	strcpy(head[0],"Direct AP");
	strcpy(head[1],"Direct PA");
	strcpy(head[2],"Inverted");
	strcpy(head[3],"Everted");	
	int n_partial;
	n_partial=Min(len_a,len_p);//partial overlap
	n_partial--;
	int n_full=1+abs(len_a-len_p)/2;// full overlap	
	int n_tot=n_full+n_partial+len_sp+1;
	char sp='S';//spacer
	char bo='P';//partial
	char in='F';//full
	FILE *out;
	if((out=fopen(file,"at"))==NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		return -1;
	}
	fprintf(out,"%d\t%s\t",mot,motif);
	int i,j;	
	for(j=n_full;j>=1;j--)fprintf(out,"\t%d%c",j,in);		
	for(j=n_partial;j>=1;j--)fprintf(out,"\t%d%c",j,bo);		
	for(j=0;j<=len_sp;j++)fprintf(out,"\t%d%c",j,sp);
	fprintf(out,"\n");
	for(i=3;i>=0;i--)
	{
		fprintf(out,"\t\t%s",head[i]);
		for(j=0;j<n_tot;j++)fprintf(out,"\t%f",100*freq[i][j]);
		fprintf(out,"\n");
	}	
	fclose(out);
	return 1;
}

#include "projoin.h" // soedinenie dvuh profiley saytov dlya odnogo pika
#include "throw_predictions.h" //permutation of sites in a peak
#include "fisher_exact_test.h" // exact fisher test

//spisok imen motivov-partnerov
struct pfm_list {
int num;
char nam[MOT_NAME_LEN];
};
#include "pfm_list.h" //spiski imen motivov-partnerov
#include "pfm_similarity.h" //permutation test for anchor/partner motif comparison (separate task for all algorithm)

int main(int argc, char *argv[])
{
	int i, j, k, m, n, n_motifs, mot, *bad_matrix; 
	char file_fasta[80], mot_db[30], mypath_data[200], prom[200], partner_db[30], file_pfm_anchor[50];	
	char ***seq;// peaks
	
	char file_hist[80], file_pval[5][80], file_pval_table[80];
	char name_anchor[MOT_NAME_LEN], name_partner[MOT_NAME_LEN];
	char xreal[]="real", xrand[]="rand", xreal_one[]="real_one";

	if(argc!=7)
	{
		printf ("%s 1file_fasta",argv[0]);//1int thresh_num_min 2int thresh_num_max
	//	printf ("4int height_permut 5int size_min_permut 6int size_max_permut 7double pvalue 8double pvalue_mult");
		printf(" 2char anchor_motif 3char partner_db 4int spacer_min 5int spacer_max 6char path_genome\n");//9char mot_anchor 
        return -1;
	}
	int thresh_num_min = 1, thresh_num_max  = 5;	// 1 5    or 5 5
	strcpy(file_fasta,argv[1]);
	int height_permut = 100, size_min_permut= 200000, size_max_permut =300000; //50000 150000 25  parametry permutacii
	double pvalue = 0.0005, pvalue_mult = 1.5; // 0.0005 1.5 parametry dlya porogov matric
	int mot_anchor = 0;// 0 = pwm from file >0 pwm from pre-computed database	
	int s_overlap_min=6, s_ncycle_small=2000, s_ncycle_large=20000;//for permutation(motif_comparison) min_length_of_alignment, no. of permutation (test & detailed)
	double s_granul=0.001;//for permutation(motif_comparison) okruglenie 4astotnyh matric	
	strcpy(file_pfm_anchor,argv[2]);		
	strcpy(partner_db,argv[3]); //hs_core, hs_full, mm_core, mm_full, dapseq
	int shift_min = atoi(argv[4]); // minimal spacer length
	int shift_max = atoi(argv[5]); // maximal spacer length
	strcpy(mypath_data,argv[6]); // ./.../hs, mm, at

	strcpy(prom,mypath_data);
	int nseq_genome, len_genome;
	if(strstr(partner_db,"core")!=NULL || strstr(partner_db,"full")!=NULL)
	{
		strcpy(mot_db,"hocomoco");
		strcat(prom,"ups2kb.plain");
		len_genome=2000;
	}
	else
	{		
		if(strstr(partner_db,"dapseq")!=NULL)
		{
			len_genome=1500;
			strcat(prom,"ups1500.plain");
			strcpy(mot_db,"dapseq");
		}
		else
		{
			printf("Error partner database %s\n", partner_db);
			return -1;
		}
	}
	if(strstr(partner_db,"hs")!=NULL)
	{		
		nseq_genome=19795;
	}
	else 
	{
		if(strstr(partner_db,"mm")!=NULL)
		{			
			nseq_genome=19991;
		}
		else 
		{
			if(strstr(partner_db,"dapseq")!=NULL)
			{				
				nseq_genome=27202;
			}
			else
			{
				printf("Partner database %s is wrong\n",partner_db);
			}
		}
	}
	strcpy(file_hist,"out_hist");
	strcpy(file_pval[0],"fisher_any_mot");
	strcpy(file_pval[1],"fisher_full_mot");
	strcpy(file_pval[2],"fisher_part_mot");
	strcpy(file_pval[3],"fisher_over_mot");
	strcpy(file_pval[4],"fisher_spac_mot");
	strcpy(file_pval_table,"out_pval");

	// nomera plohih matric, ih v ignor
    int bad_matrix_mm358[] = {172, 186, 192, 261, 287, -1};
    int bad_matrix_hs402[] = {170, 183, 190, 253, 280, 324, -1};
    int bad_matrix_at528[] = {  2, 102, 183, 184, 186, 207, 212, 213, 217, 286, 299, 300, 303, 439, -1};
    int bad_matrix_mm531[] = { 40,  61,  64,  99, 107, 108, 158, 169, 243, 253, 264, 272, 277, 283, 330, 380, 381, 424, 428, 429, 499, 529, -1};
    int bad_matrix_hs771[] = { 80,  84, 129, 139, 141, 142, 212, 230, 341, 360, 369, 379, 389, 394, 408, 443, 509, 510, 558, 564, 565, 634, 647, 719, -1};
    // Vanya vybrosil human full evx1 141, evx2 142, nfkb1 394  mouse full cdx4 40 evx1 107,evx2 108
	
	memset(name_anchor,'\0',sizeof(name_anchor));	
	if(strstr(partner_db,"hs_core") !=NULL) {
		n_motifs=NMAT_HS_CORE; bad_matrix = bad_matrix_hs402; 		
		strcpy(name_anchor,name_hs_core[0].nam);
		int nlen=strlen(name_hs_core[0].nam);name_anchor[nlen]='\0';
	}
	else
	{
		if(strstr(partner_db,"mm_core") !=NULL) {
			n_motifs=NMAT_MM_CORE;bad_matrix = bad_matrix_mm358;
			strcpy(name_anchor,name_mm_core[0].nam);
			int nlen=strlen(name_mm_core[0].nam);name_anchor[nlen]='\0';
		}
		else
		{
			if(strstr(partner_db,"mm_full") !=NULL) {
				n_motifs=NMAT_MM_FULL; bad_matrix = bad_matrix_mm531; 				
				strcpy(name_anchor,name_mm_full[0].nam);
				int nlen=strlen(name_mm_full[0].nam);name_anchor[nlen]='\0';
			}
			else
			{
				if(strstr(partner_db,"hs_full") !=NULL) {
					n_motifs=NMAT_HS_FULL; bad_matrix = bad_matrix_hs771; 					
					strcpy(name_anchor,name_hs_full[0].nam);
					int nlen=strlen(name_hs_full[0].nam);name_anchor[nlen]='\0';
				}
				else
				{
					if(strstr(partner_db,"dapseq") !=NULL) {
						n_motifs=NMAT_DAPSEQ; bad_matrix = bad_matrix_at528;						
						strcpy(name_anchor,name_dapseq[0].nam);
						int nlen=strlen(name_dapseq[0].nam);name_anchor[nlen]='\0';
					}
					else
					{
						printf("Error partner database %s\n", partner_db);
						return -1;
					}
				}
			}
		}
	}	
	double pvalue_equal = 0.01;	
	double pvalue_similarity_tot;

	int length_fasta_max=0, nseq_real=0;
	seq=NULL;
	int ftp=fasta_to_plain0(file_fasta, length_fasta_max, nseq_real);
	if(ftp==-1)
	{
		printf("File %s error 1st stage\n",file_fasta);
		return -1;
	}
	int *peak_len_real;
	peak_len_real=new int[nseq_real];
	if(peak_len_real==NULL){puts("Out of memory...");return -1;}

	seq = new char**[2];
	if(seq==NULL){puts("Out of memory...");return -1;}
	for(k=0;k<2;k++)
	{
		seq[k] = new char*[nseq_real];				
		if(seq[k]==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<nseq_real;i++)
		{			
			seq[k][i] = new char[length_fasta_max+1];
			if(seq[k][i]==NULL){puts("Out of memory...");return -1;}
		}	
	}	
	ftp=fasta_to_plain1(file_fasta, length_fasta_max, nseq_real,seq,peak_len_real);
	if(ftp==-1)
	{
		printf("File %s error 2nd stage\n",file_fasta);
		return -1;
	}		
		
	double thr[NUM_THR], thr_anchor[NUM_THR];	

	profile real[2][NUM_THR], rand[2][NUM_THR], rand_hom[NUM_THR];
	profile real_one[2], rand_one[2], rand_hom_one;
	matrices matrix, matrix0;	
	result table[NUM_THR][NUM_THR], table_one;
	combi hist_obs[NUM_THR][NUM_THR], hist_obs_one;		
	combi hist_exp[NUM_THR][NUM_THR];		
	
	//for real
	for(j=0;j<2;j++)
	{
		for(k=0;k<NUM_THR;k++)
		{
			real[j][k].nseq=nseq_real;
			real[j][k].nam=k+1;
			int ini=real[j][k].mem_in_nsit();
			if(ini==-1){puts("Not enough memory...\n");return -1;}
		}
		real_one[j].nseq=nseq_real;
		real_one[j].nam=1;
		int ini=real_one[j].mem_in_nsit();
		if(ini==-1){puts("Not enough memory...\n");return -1;}		
	}
	
	//for permutation
	srand( (unsigned)time( NULL ) );
	int nseq_rand=nseq_real*height_permut;
	if(nseq_rand<size_min_permut)height_permut=size_min_permut/nseq_real;
	if(nseq_rand>size_max_permut)height_permut=size_max_permut/nseq_real;
	nseq_rand=nseq_real*height_permut;	

	for(j=0;j<2;j++)
	{
		for(k=0;k<NUM_THR;k++)
		{
			rand[j][k].nseq=nseq_rand;
			rand[j][k].nam=k+1;
			int ini=rand[j][k].mem_in_nsit();
			if(ini==-1){puts("Not enough memory...\n");return -1;}
		}		
		rand_one[j].nseq=nseq_rand;
		rand_one[j].nam=1;
		int ini=rand_one[j].mem_in_nsit();
		if(ini==-1){puts("Not enough memory...\n");return -1;}		
	}
	for(j=0;j<NUM_THR;j++)
	{
		rand_hom[j].nseq=nseq_rand;
		rand_hom[j].nam=j+1;
		rand_hom[j].mot=0;
		int ini=rand_hom[j].mem_in_nsit();
		if(ini==-1){puts("Not enough memory...\n");return -1;}
	}		
	rand_hom_one.nseq=nseq_rand;
	rand_hom_one.nam=1;
	rand_hom_one.mot=0;
	int ini=rand_hom_one.mem_in_nsit();
	if(ini==-1){puts("Not enough memory...\n");return -1;}	

	int *peak_len_rand;
	peak_len_rand=new int[nseq_rand];
	if(peak_len_rand==NULL){puts("Out of memory...");return -1;}
	k=0;
	for(i=0;i<nseq_real;i++)
	{
		for(m=0;m<height_permut;m++)peak_len_rand[k++]=peak_len_real[i];
	}
	int *thr_err_real, *thr_err_rand, *thr_err_real_one;
	thr_err_real=new int[nseq_real];
	if(thr_err_real==NULL){puts("Out of memory...");return -1;}
	thr_err_rand=new int[nseq_rand];
	if(thr_err_rand==NULL){puts("Out of memory...");return -1;}			
	thr_err_real_one=new int[nseq_real];
	if(thr_err_real_one==NULL){puts("Out of memory...");return -1;}
	for(i=0;i<nseq_real;i++)thr_err_real_one[i]=0;

	FILE *out_hist;
	if((out_hist=fopen(file_hist,"wt"))==NULL)
	{
		printf("Input file %s can't be opened!\n", file_hist);
		return -1;
	}
	fclose(out_hist);
	FILE *out_pval_table;
	if((out_pval_table=fopen(file_pval_table,"wt"))==NULL)
	{
		printf("Input file %s can't be opened!\n", file_pval_table);
		return -1;
	}
//	fprintf(out_pval_table,"%s\n",file_fasta);
	fprintf(out_pval_table,"Motif Num");
	fprintf(out_pval_table,"\tMotif Name");
	fprintf(out_pval_table,"\tAny Pv");
	fprintf(out_pval_table,"\tFull Pv");
	fprintf(out_pval_table,"\tPart Pv");
	fprintf(out_pval_table,"\tOver Pv");
	fprintf(out_pval_table,"\tSpac Pv");
	fprintf(out_pval_table,"\t");
	fprintf(out_pval_table,"\tAny Asy");
	fprintf(out_pval_table,"\tFull Asy");
	fprintf(out_pval_table,"\tPart Asy");
	fprintf(out_pval_table,"\tOver Asy");
	fprintf(out_pval_table,"\tSpac Asy");	
	fprintf(out_pval_table,"\t");
	fprintf(out_pval_table,"\tSim Pv");	
	/*fprintf(out_pval_table,"\tSimEDpre");	
	fprintf(out_pval_table,"\tSimPCpre");	
	fprintf(out_pval_table,"\tSimED");	
	fprintf(out_pval_table,"\tSimPC");	*/
	fprintf(out_pval_table,"\n");
	fclose(out_pval_table);

	FILE *out_pval[5];
	double pval_sim[4];

	FILE *out_stat;
	if((out_stat=fopen("rec_pos.txt","wt"))==NULL)    
	{
		printf("Input file can't be opened!\n");
		return -1;
	}	
	fprintf(out_stat,"Motif Num\tMotif Name\t# Threshold\tThreshold\t%% of peaks\tRec. peaks\tTotal peaks\tRate of hits\tRec. hits\tTotal positions\n");
	for(mot=0;mot<n_motifs;mot++)
	//for(mot=0;mot<n_motifs;mot+=68)
	{	
		{
			int bad=0;
			for(i=0;bad_matrix[i]!=-1;i++)
			{
				if(mot==bad_matrix[i])
				{
					bad=1;
					break;
				}
			}
			if(bad==1)continue;
		}		
		printf("Mot %d\n",mot); 
		int len_anchor, len_partner;
		if(mot==0)
		{
			double thr_touzet[NUM_PVAL];
			double pwm_anchor[MATLEN][OLIGNUM];
			int length=pfm_to_pwm(file_pfm_anchor,&matrix);						
			if (length<=0 || length>=MATLEN)
			{
				printf("PFM to PWM conversion error, file %s\n",file_pfm_anchor);
				return -1;
			}
			for(i=0;i<length;i++)for(j=0;j<OLIGNUM;j++)pwm_anchor[i][j]=matrix.wei[i][j];			
			matrix.get_copy(&matrix0);
			int n_thr_touzet=0;
			int touzet = pwm_score_distr_granul(pvalue_equal, length, pwm_anchor, n_thr_touzet, NUM_PVAL, thr_touzet);
			if (touzet==-1)
			{
				printf("Threshold Pvalue table Touzet error\n");
				return -1;
			}
			double *fp_rate;
			fp_rate=new double [n_thr_touzet];
			if(fp_rate==NULL){puts("Out of memory...");return -1;}
			int piptd=pwm_iz_pwm_thr_dist(pwm_anchor,length,prom,n_thr_touzet,thr_touzet,fp_rate,mypath_data,nseq_genome,len_genome);
			if(piptd==-1)
			{
				printf("FP rate table error\n");
				return -1;
			}
			double fpr_select[NUM_THR];
			int stfp=select_thresholds_from_pvalues(n_thr_touzet,thr_touzet,fp_rate,pvalue,pvalue_mult,fpr_select,thr);
			delete [] fp_rate;							
			memset(name_partner,'\0',sizeof(name_partner));
			strcpy(name_partner,name_anchor);
			int nlen=strlen(name_anchor);
			name_partner[nlen]='\0';
			pvalue_similarity_tot=1E-300;
			for(i=0;i<4;i++)pval_sim[i]=pvalue_similarity_tot;
		}				
		else
		{					
			int max;
			int check=0;				
			memset(name_partner,'\0',sizeof(name_partner));
			if(strcmp(partner_db,"hs_core")==0)
			{
				max=NMAT_HS_CORE;
				if(mot>max)
				{
					printf("Max number %d of partner matrix is %d\n", mot,max);
					return -1;
				}				
				matrix.init_hs_core(mot);			
				for(j=0;j<NUM_THR;j++)thr[j]=hocomoco_thr_hs_core[mot-1][j];
				strcpy(name_partner,name_hs_core[mot].nam);
				int nlen=strlen(name_hs_core[mot].nam);
				name_partner[nlen]='\0';
				check=1;
			}
			else
			{
				if(strcmp(partner_db,"mm_core")==0)
				{
					max=NMAT_MM_CORE;
					if(mot>max)
					{
						printf("Max number %d of partner matrix is %d\n",mot,max);
						return -1;
					}
					matrix.init_mm_core(mot);			
					for(j=0;j<NUM_THR;j++)thr[j]=hocomoco_thr_mm_core[mot-1][j];
					strcpy(name_partner,name_mm_core[mot].nam);
					int nlen=strlen(name_mm_core[mot].nam);
					name_partner[nlen]='\0';
					check=1;					
				}
				else
				{
					if(strcmp(partner_db,"hs_full")==0)
					{
						max=NMAT_HS_FULL;
						if(mot>max)
						{
							printf("Max number %d of partner matrix is %d\n",mot,max);
							return -1;
						}						
						for(j=0;j<NUM_THR;j++)thr[j]=hocomoco_thr_hs_full[mot-1][j];					
						matrix.init_hs_full(mot);
						strcpy(name_partner,name_hs_full[mot].nam);
						int nlen=strlen(name_hs_full[mot].nam);
						name_partner[nlen]='\0';
						check=1;
					}
					else
					{
						if(strcmp(partner_db,"mm_full")==0)
						{
							max=NMAT_MM_FULL;
							if(mot>max)
							{
								printf("Max number %d of partner matrix is %d\n",mot,max);
								return -1;
							}							
							for(j=0;j<NUM_THR;j++)thr[j]=hocomoco_thr_mm_full[mot-1][j];						
							matrix.init_mm_full(mot);	
							strcpy(name_partner,name_mm_full[mot].nam);
							int nlen=strlen(name_mm_full[mot].nam);
							name_partner[nlen]='\0';
							check=1;
						}
						else
						{
							if(strcmp(partner_db,"dapseq")==0)
							{
								max=NMAT_DAPSEQ;
								if(mot>max)
								{
									printf("Max number %d of partner matrix is %d\n",mot,max);
									return -1;
								}								
								for(j=0;j<NUM_THR;j++)thr[j]=dapseq_thr[mot-1][j];							
								matrix.init_dapseq(mot);						
								strcpy(name_partner,name_dapseq[mot].nam);
								int nlen=strlen(name_dapseq[mot].nam);
								name_partner[nlen]='\0';
								check=1;
							}
						}
					}
				}
			}		
			if(check!=1)			
			{
				printf("Error partner database %s\n", partner_db);
				return -1;
			}			
			for(i=0;i<4;i++)pval_sim[i]=1;
			pvalue_similarity_tot=pfm_similarity(&matrix,&matrix0,s_granul,s_overlap_min,s_ncycle_small,s_ncycle_large,pval_sim);		
		}
		matrix.norm();				
		int ap;
		if(mot==0)
		{
			ap=0;
			for(j=0;j<NUM_THR;j++)thr_anchor[j]=thr[j];
		}
		else ap=1;
//int pwm_rec0(matrices *mat, double thr, int len_pro, int nseq_pro, char ***seq, profile *real)  count all sites
		////recognition 1st		
		int all_pos=0;//total number of available positions
		int wm_rec=pwm_rec0(&matrix,thr[NUM_THR-1],length_fasta_max,nseq_real,seq,&real_one[ap],all_pos);
		if(wm_rec==-1)
		{
			printf("Motif %d recognition 1st stage error\n", mot);
			return -1;
		}	
		//memory allocation for all sites
		int ini=real_one[ap].mem_in_sta();
		if(ini==-1){puts("Not enough memory...\n");return -1;}
		ini=real_one[ap].mem_in_cep();
		if(ini==-1){puts("Not enough memory...\n");return -1;}
		ini=real_one[ap].mem_in_sco();	
		if(ini==-1){puts("Not enough memory...\n");return -1;}
		real_one[ap].mot=mot;
		real_one[ap].nam=NUM_THR;
		real_one[ap].count_sites();
		//recognition 2nd
		wm_rec=pwm_rec1(&matrix,thr[NUM_THR-1],length_fasta_max,nseq_real,seq,&real_one[ap]);	
		if(wm_rec==-1)
		{
			printf("Motif %d recognition 2nd stage error\n", mot);
			return -1;
		}
		int fprint_pro=real_one[ap].fprintf_pro(mot_db,thr[NUM_THR-1],xreal);
		if(fprint_pro==-1)
		{
			printf("Real print profile error, motif %d\n",mot);
			return -1;
		}
		//count nsites for various thresholds
		for(j=0;j<NUM_THR;j++)real[ap][j].nsit[i]=0;
		for(i=0;i<nseq_real;i++)
		{
			for(k=0;k<real_one[ap].nsit[i];k++)
			{
				double sco=real_one[ap].sco[i][k];
				for(j=0;j<NUM_THR;j++)
				{
					if(sco>=thr[j])
					{
						real[ap][j].nsit[i]++;
						break;
					}
				}
			}
		}
		//memory allcation for various thresholds
		for(j=0;j<NUM_THR;j++)
		{				
			int ini=real[ap][j].mem_in_sta();
			if(ini==-1){puts("Not enough memory...\n");return -1;}
			ini=real[ap][j].mem_in_cep();
			if(ini==-1){puts("Not enough memory...\n");return -1;}
			ini=real[ap][j].mem_in_sco();	
			if(ini==-1){puts("Not enough memory...\n");return -1;}
			real[ap][j].mot=mot;
			real[ap][j].nam=j+1;
			real[ap][j].count_sites();
		}
		//initiation of profiles for various thresholds
		int rec_seq[NUM_THR], rec_pos[NUM_THR];
		for(i=0;i<NUM_THR;i++)rec_seq[i]=rec_pos[i]=0;
		for(i=0;i<nseq_real;i++)
		{
			int inx[NUM_THR];
			for(j=0;j<NUM_THR;j++)inx[j]=0;
			for(k=0;k<real_one[ap].nsit[i];k++)
			{
				double sco=real_one[ap].sco[i][k];
				for(j=0;j<NUM_THR;j++)
				{
					if(sco>=thr[j])
					{
						rec_pos[j]++;						
						real[ap][j].sta[i][inx[j]]=real_one[ap].sta[i][k];
						real[ap][j].cep[i][inx[j]]=real_one[ap].cep[i][k];
						real[ap][j].sco[i][inx[j]]=sco;
						inx[j]++;
						if(inx[j]==1)rec_seq[j]++;
						break;
					}
				}
			}
		}
		for(j=0;j<NUM_THR;j++)
		{
			fprintf(out_stat,"%d\t%s\t%d\t%f\t",mot,name_partner,j+1,thr[j]);
			fprintf(out_stat,"%f\t%d\t%d\t",100*(double)rec_seq[j]/nseq_real,rec_seq[j],nseq_real);
			fprintf(out_stat,"%g\t%d\t%d\n",(double)rec_pos[j]/all_pos,rec_pos[j],all_pos);
		}
		for(j=0;j<NUM_THR;j++)
		{
			real[ap][j].count_sites();
			int fprint_pro=real[ap][j].fprintf_pro(mot_db,thr[j],xreal);
			if(fprint_pro==-1)
			{
				printf("Real print profile error, motif %d\n",mot);
				return -1;
			}
		}
		if(ap==0){len_anchor=len_partner=matrix.len;}
		else len_partner=matrix.len;
		matrix.mem_out(matrix.len);
		//int len_ap=matrix.len;					
		/*int test;
		test = real_one[ap].test();
		if(test==-1)
		{
			printf("Real One %d error Mot %d Thr One\n",ap, mot);
			return -1;
		}*/
	
//one thresh  rand - &rand_hom_one,&rand_one[ap]   real - real_one[0],real_one[ap]
//many thresh rand - &rand_hom[j],&rand[ap][k]     real - real[0][j],real[ap][k]
		//anchor		
		{
			int cop=real_one[0].get_copy_rand(&rand_hom_one,height_permut);
			if(cop==-1)
			{
				printf("Rand Copy error Mot %d Thr One\n",mot);
				return -1;
			}
			for(j=0;j<NUM_THR;j++)
			{				
				cop=real[0][j].get_copy_rand(&rand_hom[j],height_permut);
				if(cop==-1)
				{
					printf("Rand Copy error Mot %d Thr%d\n",mot,j);
					return -1;
				}
			}
			/*
			test = rand_hom_one.test();
			if(test==-1)
			{
				printf("Rand Hom One error Mot %d Thr One\n", mot);
				return -1;
			}
			for(j=0;j<NUM_THR;j++)
			{
				test = rand_hom[j].test();
				if(test==-1)
				{
					printf("Rand Hom error Mot %d Thr %d\n",mot,j+1);
					return -1;
				}
			}*/

		}				
		//partner
		{
			int cop=real_one[ap].get_copy_rand(&rand_one[ap],height_permut);
			if(cop==-1)
			{
				printf("Rand Copy error Mot %d Thr One\n",mot);
				return -1;
			}
			for(j=0;j<NUM_THR;j++)
			{				
				cop=real[ap][j].get_copy_rand(&rand[ap][j],height_permut);
				if(cop==-1)
				{
					printf("Rand Copy error Mot %d Thr%d\n",mot,j);
					return -1;
				}
			}
			/*
			test = rand_one[ap].test();
			if(test==-1)
			{
				printf("Rand One %d error Mot %d Thr One\n",ap, mot);
				return -1;
			}
			for(j=0;j<NUM_THR;j++)
			{
				test = rand[ap][j].test();
				if(test==-1)
				{
					printf("Rand %d error Mot %d Thr %d\n",ap, mot,j+1);
					return -1;
				}
			}*/
		}
		/*for(j=0;j<NUM_THR;j++)
		{
			int fprint_pro=rand[ap][j].fprintf_pro(mot_db,thr[j],"rand_do");
			if(fprint_pro==-1)
			{
				printf("Real print profile error, motif %d\n",mot);
				return -1;
			}
		}*/		
//		char file_throw_err[80], file_throw_err0[80];
	//	FILE *out_nsit_throw;
		{// one threshold
			for(m=0;m<nseq_real;m++)thr_err_real[m]=0;			
			rand_one[ap].mot=mot;			
			int throwp=throw_predictions(peak_len_rand,&rand_hom_one,&rand_one[ap],len_anchor,len_partner,0,thr_err_real,nseq_real,nseq_rand,seq[0],height_permut);			
			if(throwp==-1)
			{
				printf("Throw Prediction error One - Anc 0 Par %d\n", mot);
				return -1;
			}
		/*	int fprint_pro=rand_hom_one.fprintf_pro(mot_db,thr_anchor[NUM_THR-1],"rand_po");
			if(fprint_pro==-1)
			{
				printf("Real print profile error, motif %d\n",mot);
				return -1;
			}
			fprint_pro=rand_one[ap].fprintf_pro(mot_db,thr[NUM_THR-1],"rand_po");
			if(fprint_pro==-1)
			{
				printf("Real print profile error, motif %d\n",mot);
				return -1;
			}*/
		/*	memset(file_throw_err,'\0',sizeof(file_throw_err));
			memset(file_throw_err0,'\0',sizeof(file_throw_err0));
			strcpy(file_throw_err,"throw_nsit_mot0");
			char buf[10];
			sprintf(buf,"%d",mot);
			strcat(file_throw_err,buf);
			strcpy(file_throw_err0,file_throw_err);
			strcat(file_throw_err,"_one.txt");
			if((out_nsit_throw=fopen(file_throw_err,"wt"))==NULL)
			{
				printf("Input file %s can't be opened!\n", file_throw_err);
				return -1;
			}
			for(m=0;m<nseq_real;m++)fprintf(out_nsit_throw,"%d\n",thr_err_real[m]);
			fclose(out_nsit_throw);*/
			table_one.ini();
			real_one[ap].mot=mot;
			k=0;
			for(i=0;i<nseq_real;i++)
			{
				for(j=0;j<height_permut;j++)
				{
					thr_err_rand[k++]=thr_err_real[i];
				}		
			}					
			int proj=projoin(xreal_one,mot_db,real_one[0],real_one[ap],shift_min,shift_max,len_anchor,len_partner,thr_err_real,nseq_real,seq,&table_one.observed,&hist_obs_one,peak_len_real);		
			if(proj==-1)
			{
				printf("Projoin Real error Anc 0 Par %d\n", mot);
				return -1;
			}
			hist_obs_one.fprintf_all(file_hist,mot,name_partner,len_anchor,len_partner,shift_max);
		}// one threshold
		//many thresholds
		//strcat(file_throw_err0,"_thr");
		for(j=0;j<NUM_THR;j++)
		{
			for(k=0;k<NUM_THR;k++)
			{
				//printf("Anchor %d thresh, Partner %d thresh\n",j+1,k+1);
				for(m=0;m<nseq_real;m++)thr_err_real[m]=0;				
				int throwp=throw_predictions(peak_len_rand,&rand_hom[j],&rand[ap][k],len_anchor,len_partner,0,thr_err_real,nseq_real,nseq_rand,seq[0],height_permut);
				if(throwp==-1)
				{
					printf("Throw Prediction error Anc 0(%d) Par %d(%d)\n", j, mot, k);
					return -1;
				}
				n=0;
				for(i=0;i<nseq_real;i++)
				{
					for(m=0;m<height_permut;m++)
					{
						thr_err_rand[n++]=thr_err_real[i];
					}		
				}	
			/*	strcpy(file_throw_err,file_throw_err0);
				char buf[10];
				sprintf(buf,"%d",j+1);
				strcat(file_throw_err,buf);
				sprintf(buf,"%d",k+1);
				strcat(file_throw_err,buf);
				strcat(file_throw_err,".txt");	
				if((out_nsit_throw=fopen(file_throw_err,"wt"))==NULL)
				{
					printf("Input file %s can't be opened!\n", file_throw_err);
					return -1;
				}
				for(m=0;m<nseq_real;m++)fprintf(out_nsit_throw,"%d\n",thr_err_real[m]);
				fclose(out_nsit_throw);*/
				/*	int fprint_pro=rand_hom[j].fprintf_pro(mot_db,thr_anchor[j],"rand_po");
				if(fprint_pro==-1)
				{
					printf("Real print profile error, motif %d\n",mot);
					return -1;
				}
				fprint_pro=rand[ap][k].fprintf_pro(mot_db,thr[k],"rand_po");
				if(fprint_pro==-1)
				{
					printf("Real print profile error, motif %d\n",mot);
					return -1;
				}*/
				table[j][k].ini();
//				hist_exp[j][k].ini(len_anchor,len_partner,shift);
				int proj=projoin(xrand,mot_db,rand_hom[j],rand[ap][k],shift_min,shift_max,len_anchor,len_partner,thr_err_rand,nseq_rand,seq,&table[j][k].expected,&hist_exp[j][k],peak_len_rand);		
				if(proj==-1)
				{
					printf("Projoin Rand error Anc 0 Par %d\n", mot);
					return -1;
				}
	//			hist_obs[j][k].ini(len_anchor,len_partner,shift);
				proj=projoin(xreal,mot_db,real[0][j],real[ap][k],shift_min,shift_max,len_anchor,len_partner,thr_err_real,nseq_real,seq,&table[j][k].observed,&hist_obs[j][k],peak_len_real);		
				if(proj==-1)
				{
					printf("Projoin Real error Anc 0 Par %d\n", mot);
					return -1;
				}
			}
		}
		for(j=0;j<NUM_THR;j++)
		{
			for(k=0;k<NUM_THR;k++)
			{			
				//printf("Mot %d J %d K %d\n",mot,j+1,k+1);
				int fisher;
				fisher=fisher_exact_test(table[j][k].observed.any, table[j][k].observed.two_sites,table[j][k].expected.any,table[j][k].expected.two_sites,table[j][k].fold_a,table[j][k].pvalue_a);					
				if(fisher==-1)
				{
					printf("Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher=fisher_exact_test(table[j][k].observed.full, table[j][k].observed.two_sites,table[j][k].expected.full,table[j][k].expected.two_sites,table[j][k].fold_f,table[j][k].pvalue_f);					
				if(fisher==-1)
				{
					printf("Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher=fisher_exact_test(table[j][k].observed.partial, table[j][k].observed.two_sites,table[j][k].expected.partial,table[j][k].expected.two_sites,table[j][k].fold_p,table[j][k].pvalue_p);					
				if(fisher==-1)
				{
					printf("Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher=fisher_exact_test(table[j][k].observed.over, table[j][k].observed.two_sites,table[j][k].expected.over,table[j][k].expected.two_sites,table[j][k].fold_o,table[j][k].pvalue_o);					
				if(fisher==-1)
				{
					printf("Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher=fisher_exact_test(table[j][k].observed.spacer, table[j][k].observed.two_sites,table[j][k].expected.spacer,table[j][k].expected.two_sites,table[j][k].fold_s,table[j][k].pvalue_s);					
				if(fisher==-1)
				{
					printf("Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
			}
		}
		for(i=0;i<5;i++)
		{
			char file_pval0[80];
			strcpy(file_pval0,file_pval[i]);
			char buf[10];
			sprintf(buf,"%d",mot);					
			strcat(file_pval0,buf);
			if((out_pval[i]=fopen(file_pval0,"wt"))==NULL)
			{
				printf("Input file %s can't be opened!\n", file_pval0);
				return -1;
			}
			fprintf(out_pval[i],"Anchor Thr\tPartner Thr\t\tReal CE+\tReal Total\tRand CE+\tRand Total\tFold\tP-value\n");
		}
		for(j=0;j<NUM_THR;j++)
		{
			for(k=0;k<NUM_THR;k++)
			{		
				for(i=0;i<5;i++)fprintf(out_pval[i],"A %d\tP %d\t\t",j+1,k+1);												
				fprintf(out_pval[0],"%d\t%d\t%d\t%d\t%f\t%g\n",table[j][k].observed.any, table[j][k].observed.two_sites,table[j][k].expected.any,table[j][k].expected.two_sites,table[j][k].fold_a,table[j][k].pvalue_a);
				fprintf(out_pval[1],"%d\t%d\t%d\t%d\t%f\t%g\n",table[j][k].observed.full, table[j][k].observed.two_sites,table[j][k].expected.full,table[j][k].expected.two_sites,table[j][k].fold_f,table[j][k].pvalue_f);
				fprintf(out_pval[2],"%d\t%d\t%d\t%d\t%f\t%g\n",table[j][k].observed.partial, table[j][k].observed.two_sites,table[j][k].expected.partial,table[j][k].expected.two_sites,table[j][k].fold_p,table[j][k].pvalue_p);
				fprintf(out_pval[3],"%d\t%d\t%d\t%d\t%f\t%g\n",table[j][k].observed.over, table[j][k].observed.two_sites,table[j][k].expected.over,table[j][k].expected.two_sites,table[j][k].fold_o,table[j][k].pvalue_o);
				fprintf(out_pval[4],"%d\t%d\t%d\t%d\t%f\t%g\n",table[j][k].observed.spacer, table[j][k].observed.two_sites,table[j][k].expected.spacer,table[j][k].expected.two_sites,table[j][k].fold_s,table[j][k].pvalue_s);
			}
		}
		for(i=0;i<5;i++)fclose(out_pval[i]);									
		double pval_tot_min[5]={0,0,0,0,0};
		double limit=300;
		double pv_limit=1E-300;
		double lgpv[5][NUM_THR][NUM_THR];
		for(j=0;j<NUM_THR;j++)
		{
			for(k=0;k<NUM_THR;k++)
			{
				for(i=0;i<5;i++)lgpv[i][k][j]=0;
				if(table[k][j].fold_a>1)
				{
					if(table[k][j].pvalue_a<=pv_limit)lgpv[0][k][j]=limit;
					else lgpv[0][k][j]=-log10(table[k][j].pvalue_a);
				}
				if(table[k][j].fold_f>1)
				{
					if(table[k][j].pvalue_f<=pv_limit)lgpv[1][k][j]=limit;
					else lgpv[1][k][j]=-log10(table[k][j].pvalue_f);
				}
				if(table[k][j].fold_p>1)
				{
					if(table[k][j].pvalue_p<=pv_limit)lgpv[2][k][j]=limit;
					else lgpv[2][k][j]=-log10(table[k][j].pvalue_p);
				}
				if(table[k][j].fold_o>1)
				{
					if(table[k][j].pvalue_o<=pv_limit)lgpv[3][k][j]=limit;
					else lgpv[3][k][j]=-log10(table[k][j].pvalue_o);
				}
				if(table[k][j].fold_s>1)
				{
					if(table[k][j].pvalue_s<=pv_limit)lgpv[4][k][j]=limit;
					else lgpv[4][k][j]=-log10(table[k][j].pvalue_s);
				}
			}
		}
		for(i=0;i<5;i++)
		{
			int limit_found=0;
			int mnoj=NUM_THR*NUM_THR;
			for(j=0;j<NUM_THR;j++)
			{
				for(k=0;k<NUM_THR;k++)
				{				
					double val=lgpv[i][k][j];
					if(val==limit)
					{
						limit_found=1;
						break;
					}
				}
				if(limit_found==1)break;
			}
			if(limit_found==1)
			{
				pval_tot_min[i]=pv_limit;
			}
			else
			{
				for(j=0;j<NUM_THR;j++)
				{
					for(k=0;k<NUM_THR;k++)
					{				
						double val=lgpv[i][k][j];							
						if(val>pval_tot_min[i])pval_tot_min[i]=val;				
					}
				}
				pval_tot_min[i]=pow(10,-pval_tot_min[i]);
			}
		}
		//double asy[5]={0,0,0,0,0};
/*		for(j=0;j<NUM_THR;j++)
		{
			for(k=j+1;k<NUM_THR;k++)
			{
				for(i=0;i<5;i++)asy[i]+=(lgpv[i][j][k]-lgpv[i][k][j]);
			}
		}
		{
			int mnoj=NUM_THR*(NUM_THR-1)/2;
			for(i=0;i<5;i++)
			{
				asy[i]/=mnoj;					
			}
		}*/
		if((out_pval_table=fopen(file_pval_table,"at"))==NULL)
		{
			printf("Input file %s can't be opened!\n", file_pval_table);
			return -1;
		}
		fprintf(out_pval_table,"0,%d",mot);
		fprintf(out_pval_table,"\t%s",name_partner);
		for(i=0;i<5;i++)fprintf(out_pval_table,"\t%g",pval_tot_min[i]);
		fprintf(out_pval_table,"\t");
		//for(i=0;i<5;i++)fprintf(out_pval_table,"\t%f",asy[i]);
		//fprintf(out_pval_table,"\t\t");
		if(mot!=0)fprintf(out_pval_table,"%g",pvalue_similarity_tot);
		//for(i=0;i<4;i++)fprintf(out_pval_table,"\t%g",pval_sim[i]);
		fprintf(out_pval_table,"\n");
		fclose(out_pval_table);
		if(ap==1)
		{
			for(j=0;j<NUM_THR;j++)
			{
				real[ap][j].mem_out_sta();
				real[ap][j].mem_out_cep();
				real[ap][j].mem_out_sco();
				for(i=0;i<nseq_real;i++)real[ap][j].nsit[i]=0;
				rand[ap][j].mem_out_sta();
				rand[ap][j].mem_out_cep();									
				for(i=0;i<nseq_rand;i++)rand[ap][j].nsit[i]=0;
			}
			real_one[ap].mem_out_sta();
			real_one[ap].mem_out_cep();
			real_one[ap].mem_out_sco();
			for(i=0;i<nseq_real;i++)real_one[ap].nsit[i]=0;
			rand_one[ap].mem_out_sta();
			rand_one[ap].mem_out_cep();			
			for(i=0;i<nseq_rand;i++)rand_one[ap].nsit[i]=0;
		}
	}		
	delete [] thr_err_real;
	delete [] thr_err_rand;
	delete [] peak_len_real;
	delete [] peak_len_rand;
	for(k=0;k<2;k++)
	{
		for(i=0;i<nseq_real;i++)
		{			
			delete [] seq[k][i];			
		}	
		delete [] seq[k];
	}
	delete [] seq;
	real_one[0].mem_out_sta();
	real_one[0].mem_out_cep();
	real_one[0].mem_out_sco();
	rand_one[0].mem_out_sta();
	rand_one[0].mem_out_cep();
	for(j=0;j<NUM_THR;j++)
	{
		real[0][j].mem_out_sta();
		real[0][j].mem_out_cep();
		real[0][j].mem_out_sco();
	}
	for(j=0;j<NUM_THR;j++)
	{
		rand[0][j].mem_out_sta();
		rand[0][j].mem_out_cep();	
	}
	for(j=0;j<NUM_THR;j++)
	{
		rand_hom[j].mem_out_sta();
		rand_hom[j].mem_out_cep();
	}
	rand_hom_one.mem_out_sta();
	rand_hom_one.mem_out_cep();
	for(j=0;j<2;j++)
	{
		real_one[j].mem_out_nsit();
		for(k=0;k<NUM_THR;k++)
		{
			real[j][k].mem_out_nsit();
		}
	}
	for(j=0;j<2;j++)
	{
		rand_one[j].mem_out_nsit();
		for(k=0;k<NUM_THR;k++)
		{
			rand[j][k].mem_out_nsit();
		}
	}
	matrix0.mem_out(matrix0.len);
	return 1;
}