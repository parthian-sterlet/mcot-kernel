#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 12000
#define MATLEN 50 //max matrix length
#define SPACLEN 100 //max spacer length
#define OLIGNUM 4// di 16 mono 4
#define NUM_THR 5 //4islo porogov

int StrNStr(char *str, char c, int n)
{
	int i, len = strlen(str);
	int k = 1;
	for (i = 0; i < len; i++)
	{
		if (str[i] == c)
		{
			if (k == n)return i;
			k++;
		}
	}
	return -1;
}
char *TransStr(char *d)
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i < lens; i++)
	{
		c = int(d[i]);
		if (c < 97) d[i] = char(c + 32);
		//else break;
	}
	return(d);
}
char *TransStrBack(char *d)
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i < lens; i++)
	{
		c = int(d[i]);
		if (c >= 97) d[i] = char(c - 32);
		//else break;
	}
	return(d);
}
void DelChar(char *str, char c)
{
	int i, lens, size;

	size = 0;
	lens = strlen(str);
	for (i = 0; i < lens; i++)
	{
		if (str[i] != c)str[size++] = str[i];
	}
	str[size] = '\0';
}
int ComplStr(char *d)
{
	char *d1;
	int i, len;
	len = strlen(d);
	d1 = new char[len + 1];
	if (d1 == NULL) return 0;
	strcpy(d1, d);
	//	memset(d,0,sizeof(d));
	for (i = 0; i < len; i++)
	{
		switch (d1[len - i - 1])
		{
		case 'a': {d[i] = 't'; break; }
		case 't': {d[i] = 'a'; break; }
		case 'c': {d[i] = 'g'; break; }
		case 'g': {d[i] = 'c'; break; }
		case 'A': {d[i] = 'T'; break; }
		case 'T': {d[i] = 'A'; break; }
		case 'C': {d[i] = 'G'; break; }
		case 'G': {d[i] = 'C'; break; }
		case 'N': {d[i] = 'N'; break; }
		case 'n': {d[i] = 'n'; break; }
		default: d[i] = 'n';
		}
	}
	delete[] d1;
	return 1;
}
void Mix(int *a, int *b)
{
	int buf = *a;
	*a = *b;
	*b = buf;
}
void Mix(char *a, char *b)
{
	char buf = *a;
	*a = *b;
	*b = buf;
}
void BigMix1(char *d)
{
	int r;
	int len = strlen(d);
	for (r = 0; r < len - 1; r++) Mix(&d[r], &d[1 + r + (rand() % (len - 1 - r))]);
}
void BigMix1(int *d1, int len) // pereme6ivanie stroki
{
	int r, s;
	for (r = 0; r < len; r++)
	{
		s = rand() % len;
		if (s != r)Mix(&d1[r], &d1[s]);
	}
}
#include "fasta_to_plain.h"
#include "select_thresholds_from_pvalues.h"

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
	int **cel;// popadanie v interval porogov
	double **sco;
	double **pv;
	int mem_in_sta(void);
	int mem_in_cep(void);
	int mem_in_cel(void);
	int mem_in_sco(void);
	int mem_in_pv(void);
	int mem_in_nsit(void);
	void mem_out_sta(void);
	void mem_out_cep(void);
	void mem_out_cel(void);
	void mem_out_sco(void);
	void mem_out_pv(void);
	void mem_out_nsit(void);
	int get_copy_rand(profile *a, int height);
	int clear_real(void);
	int fprintf_pro(char *mot_db, double thr, char *mode);//mot_db real/rand
	void count_sites(void);
	int test(void);
};
int profile::mem_in_sta(void)
{
	int i;
	sta = new int*[nseq];
	if (sta == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		if (nsit[i] > 0)
		{
			sta[i] = new int[nsit[i]];
			if (sta[i] == NULL) return -1;
		}
	}
	return 1;
}
int profile::mem_in_cep(void)
{
	int i;
	cep = new char*[nseq];
	if (cep == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		if (nsit[i] > 0)
		{
			cep[i] = new char[nsit[i]];
			if (cep[i] == NULL) return -1;
		}
	}
	return 1;
}
int profile::mem_in_cel(void)
{
	int i;
	cel = new int*[nseq];
	if (cep == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		if (nsit[i] > 0)
		{
			cel[i] = new int[nsit[i]];
			if (cel[i] == NULL) return -1;
		}
	}
	return 1;
}
int profile::mem_in_sco(void)
{
	int i;
	sco = new double*[nseq];
	if (sco == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		if (nsit[i] > 0)
		{
			sco[i] = new double[nsit[i]];
			if (sco[i] == NULL) return -1;
		}
	}
	return 1;
}
int profile::mem_in_pv(void)
{
	int i;
	pv = new double*[nseq];
	if (pv == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		if (nsit[i] > 0)
		{
			pv[i] = new double[nsit[i]];
			if (pv[i] == NULL) return -1;
		}
	}
	return 1;
}
void profile::mem_out_sta(void)
{
	int i;
	for (i = 0; i < nseq; i++)
	{
		if (nsit[i] > 0)delete[] sta[i];
	}
	delete[] sta;
	sta = NULL;
}
void profile::mem_out_cep(void)
{
	int i;
	for (i = 0; i < nseq; i++)if (nsit[i] > 0)delete[] cep[i];
	delete[] cep;
	cep = NULL;
}
void profile::mem_out_cel(void)
{
	int i;
	for (i = 0; i < nseq; i++)if (nsit[i] > 0)delete[] cel[i];
	delete[] cel;
	cel = NULL;
}
void profile::mem_out_sco(void)
{
	int i;
	for (i = 0; i < nseq; i++)if (nsit[i] > 0)delete[] sco[i];
	delete[] sco;
	sco = NULL;
}
void profile::mem_out_pv(void)
{
	int i;
	for (i = 0; i < nseq; i++)if (nsit[i] > 0)delete[] pv[i];
	delete[] pv;
	pv = NULL;
}
int profile::mem_in_nsit(void)
{
	nsit = new int[nseq];
	if (nsit == NULL) return -1;
	int i;
	for (i = 0; i < nseq; i++)nsit[i] = 0;
	return 1;
}
void profile::mem_out_nsit(void)
{
	delete[] nsit;
}
int profile::get_copy_rand(profile *a, int height)
{
	int i, j, h, z, ini;
	a->mot = mot;
	a->nam = nam;
	a->nseq = height * nseq;
	a->nsit_all = height * nsit_all;
	a->nseq_rec = height * nseq_rec;
	z = 0;
	for (i = 0; i < nseq; i++)
	{
		int nsi = nsit[i];
		for (h = 0; h < height; h++)
		{
			a->nsit[z++] = nsi;
		}
	}
	ini = a->mem_in_sta();
	if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	ini = a->mem_in_cep();
	if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	ini = a->mem_in_cel();
	if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	ini = a->mem_in_pv();
	if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	z = 0;
	for (i = 0; i < nseq; i++)
	{
		for (h = 0; h < height; h++)
		{
			for (j = 0; j < a->nsit[z]; j++)
			{
				a->sta[z][j] = sta[i][j];
				a->cep[z][j] = cep[i][j];
				a->cel[z][j] = cel[i][j];
				a->pv[z][j] = pv[i][j];
			}
			z++;
		}
	}
	return 1;
}
int profile::clear_real(void)
{
	/*	if (mem_clear!=0)
	{
	mem_out_sta();
	mem_out_cep();
	mem_out_cel();
	}*/
	int ini = mem_in_sta();
	if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	ini = mem_in_cep();
	if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	ini = mem_in_cel();
	if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	ini = mem_in_sco();
	if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	ini = mem_in_pv();
	if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	return 1;
}
int profile::fprintf_pro(char *mot_db, double thr, char *mode)
{
	//int print_sco=1;
	//if(strncmp(mode,"real",4)==0)print_sco=1;
	//	else print_sco=0;
	int i, j;
	char fileo[80];
	FILE *out;
	memset(fileo, '\0', sizeof(fileo));
	strcpy(fileo, mode);
	strcat(fileo, "_");
	strcat(fileo, mot_db);//real or random	
	char buf[10];
	sprintf(buf, "%d", mot);
	strcat(fileo, buf);//hocomoco or dapseq		
	strcat(fileo, "_thr");
	memset(buf, '\0', sizeof(buf));
	sprintf(buf, "%d", nam);
	strcat(fileo, buf);//nomer poroga	
	if ((out = fopen(fileo, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", fileo);
		return -1;
	}
	for (i = 0; i < nseq; i++)
	{
		fprintf(out, ">Seq %d\tThr %f\tNsites %d\n", i + 1, thr, nsit[i]);
		for (j = 0; j < nsit[i]; j++)
		{
			fprintf(out, "%d\t", sta[i][j]);
			fprintf(out, "%f", pv[i][j]);
			fprintf(out, "\t");
			fprintf(out, "%c\n", cep[i][j]);;
		}
	}
	fclose(out);
	return 1;
}
void profile::count_sites(void)
{
	nsit_all = nseq_rec = 0;
	int i;
	for (i = 0; i < nseq; i++)
	{
		int ns = nsit[i];
		if (ns > 0)
		{
			nsit_all += ns;
			nseq_rec++;
		}
	}
}
int profile::test(void)
{
	int i, j;
	for (i = 0; i < nseq; i++)
	{
		for (j = 0; j < nsit[i]; j++)
		{
			if (sta[i][j] > 5 * SEQLEN || sta[i][j] < 0)return -1;
		}
	}
	return 1;
}
int UnderStol(char *str, int nstol, char *ret, char sep)
{
	memset(ret, 0, sizeof(ret));
	int p1, p2, len;
	if (nstol == 0)
	{
		p2 = StrNStr(str, sep, 1);
		if (p2 == -1)p2 = strlen(str);
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
		}
		if (p1 == -1 || p2 == -1) return 0;
		len = p2 - p1 - 1;
		strncpy(ret, &str[p1 + 1], len);
		ret[len] = '\0';
		return 1;
	}
}
#include "rec.h"

//dlya pods4eta zna4imostey CE, full_overlap, partial_overlap, spacer (overlap = partial_overlap OR full_overlap)
struct count {
	int two_sites;
	int any;
	int partial;
	int full;
	int overlap;
	int spacer;
	void ini(void);
};
void count::ini(void)
{
	two_sites = any = partial = full = overlap = spacer = 0;
}
struct result {
	count cell[NUM_THR][NUM_THR];
	count anc;
	count par;
	count sit;
	count anc_sit;
	count par_sit;
	void ini(void);
} observed, expected;
void result::ini(void)
{
	int j, k;
	for (j = 0; j < NUM_THR; j++)for (k = 0; k < NUM_THR; k++)cell[j][k].ini();
	anc.ini();
	par.ini();
	sit.ini();
	anc_sit.ini();
	par_sit.ini();
}

/*
int i,j;
*/
double pvalue_a[NUM_THR][NUM_THR];
double pvalue_p[NUM_THR][NUM_THR];
double pvalue_f[NUM_THR][NUM_THR];
double pvalue_s[NUM_THR][NUM_THR];
double pvalue_o[NUM_THR][NUM_THR];

struct pval {
	double anchor;
	double partner;
	double anc_par;
	double fold_anc_par;
	void ini(void);
} pv_any, pv_full, pv_partial, pv_overlap, pv_spacer;
void pval::ini(void)
{
	anchor = partner = anc_par = 1;
	fold_anc_par = 1;
}
// dlya vyvoda histogram of CE distribution as fanction of mutual orientation and location of anchor/partner motifs
struct combi {
	//	int n_partial;
	//	int n_full;
	//int n_tot;
	double freq[4][MATLEN + SPACLEN];
	double freqa[MATLEN + SPACLEN];// sum of all orientation
	double freqc[MATLEN + SPACLEN];// sum of all orientation
	//	void ini(int len_a, int len_p, int len_sp);
	//	void mem_out(void);
	int fprintf_all(char *file, int mot, char *motif, int len_a, int len_p, int len_sp, char *mode);
};
int combi::fprintf_all(char *file, int mot, char *motif, int len_a, int len_p, int len_sp, char *mode)
{
	char head[6][12];
	strcpy(head[5], "Cumulative");
	strcpy(head[4], "Any");
	strcpy(head[0], "DirectAP");
	strcpy(head[1], "DirectPA");
	strcpy(head[2], "Inverted");
	strcpy(head[3], "Everted");	
	char head0[10];
	strcpy(head0, "DirectAA");
	int n_partial;
	n_partial = Min(len_a, len_p);//partial overlap
	n_partial--;
	int n_full = 1 + abs(len_a - len_p) / 2;// full overlap	
	int n_tot = n_full + n_partial + len_sp + 1;
	char sp = 'S';//spacer
	char bo = 'P';//partial
	char in = 'F';//full
	FILE *out;
	if ((out = fopen(file, mode)) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		return -1;
	}
	fprintf(out, "%d\t%s\t", mot, motif);
	int i, j;
	for (j = n_full; j >= 1; j--)fprintf(out, "\t%d%c", j - 1, in);
	for (j = n_partial; j >= 1; j--)fprintf(out, "\t%d%c", j, bo);
	for (j = 0; j <= len_sp; j++)fprintf(out, "\t%d%c", j, sp);
	fprintf(out, "\n");
	if (mot > 0)
	{
		for (i = 3; i >= 0; i--)
		{
			fprintf(out, "\t\t%s", head[i]);
			for (j = 0; j < n_tot; j++)fprintf(out, "\t%f", 100 * freq[i][j]);
			fprintf(out, "\n");
		}
	}
	else
	{
		for (i = 3; i >= 1; i--)
		{
			if (i != 1)fprintf(out, "\t\t%s", head[i]);
			else fprintf(out, "\t\t%s", head0);
			for (j = 0; j < n_tot; j++)fprintf(out, "\t%f", 100 * freq[i][j]);
			fprintf(out, "\n");
		}
	}
	fprintf(out, "\t\t%s", head[4]);
	for (j = 0; j < n_tot; j++)fprintf(out, "\t%f", 100 * freqa[j]);
	fprintf(out, "\n");
	fprintf(out, "\t\t%s", head[5]);
	for (j = 0; j < n_tot; j++)fprintf(out, "\t%f", 100 * freqc[j]);
	fprintf(out, "\n");
	fclose(out);
	return 1;
}
struct asy_plot {
	int **full;
	int **partial;
	int **overlap;
	int **spacer;
	int **any;
	int total_full;
	int total_partial;
	int total_overlap;
	int total_spacer;
	int total_any;
	double min;
	double max;
	double inter;
	int n_karman;
	int mem_in(int karman, double xmin, double xmax, double xinter);
	void mem_out(void);
	void zero(void);
	void sum(void);
};
int asy_plot::mem_in(int karman, double xmin, double xmax, double xinter)
{
	n_karman = karman - 1;//max karman
	min = xmin, max = xmax, inter = xinter;
	int i;
	full = new int*[karman];
	if (full == NULL) return -1;
	for (i = 0; i < karman; i++)
	{
		full[i] = new int[karman];
		if (full[i] == NULL) return -1;
	}
	partial = new int*[karman];
	if (partial == NULL) return -1;
	for (i = 0; i < karman; i++)
	{
		partial[i] = new int[karman];
		if (partial[i] == NULL) return -1;
	}
	overlap = new int*[karman];
	if (overlap == NULL) return -1;
	for (i = 0; i < karman; i++)
	{
		overlap[i] = new int[karman];
		if (overlap[i] == NULL) return -1;
	}
	spacer = new int*[karman];
	if (spacer == NULL) return -1;
	for (i = 0; i < karman; i++)
	{
		spacer[i] = new int[karman];
		if (spacer[i] == NULL) return -1;
	}
	any = new int*[karman];
	if (any == NULL) return -1;
	for (i = 0; i < karman; i++)
	{
		any[i] = new int[karman];
		if (any[i] == NULL) return -1;
	}
	return 1;
}
void asy_plot::zero(void)
{
	int i, j;
	for (i = 0; i <= n_karman; i++)for (j = 0; j <= n_karman; j++)full[i][j] = partial[i][j] = overlap[i][j] = spacer[i][j] = any[i][j] = 0;
}
void asy_plot::sum(void)
{
	int i, j;
	total_full = total_partial = total_overlap = total_spacer = total_any = 0;
	for (i = 0; i <= n_karman; i++)
	{
		for (j = 0; j <= n_karman; j++)
		{
			total_full += full[i][j];
			total_partial += partial[i][j];
			total_overlap += overlap[i][j];
			total_spacer += spacer[i][j];
			total_any += any[i][j];
		}
	}
}
void asy_plot::mem_out(void)
{
	int i;
	for (i = 0; i <= n_karman; i++) delete[] full[i];
	delete[] full;
	for (i = 0; i <= n_karman; i++) delete[] partial[i];
	delete[] partial;
	for (i = 0; i <= n_karman; i++) delete[] overlap[i];
	delete[] overlap;
	for (i = 0; i <= n_karman; i++) delete[] spacer[i];
	delete[] spacer;
	for (i = 0; i <= n_karman; i++) delete[] any[i];
	delete[] any;
}
#include "projoin.h"
#include "throw_predictions.h"
#include "fisher_exact_test.h"

int main(int argc, char *argv[])
{
	int i, j, k, m, mot;
	char file_fasta[120], file_profile[2][120], file_table[2][120];
	char ***seq;// peaks
	char file_hist[120], file_pval[5][120], file_pval_table[120], file_hist_rand[120];
	char name[2][120], name_fasta[120];
	char xreal[] = "real", xrand[] = "rand", xreal_one[] = "real_one";
	char file_fpr[2][80];
	strcpy(file_fpr[0], "fpr_anchor.txt");
	strcpy(file_fpr[1], "fpr_partner.txt");

	if (argc != 10)
	{
		printf("%s 1file_fasta", argv[0]);//1int thresh_num_min 2int thresh_num_max
		printf("2 motif1.profile 3motif2.profile 4int motif1.length 5int motif1.length 6int motif1.table_thr_fpr 7int motif1.table_thr_fpr");
		printf(" 8int spacer_min 9int spacer_max\n");//9char mot_anchor 
		return -1;
	}
	int thresh_num_min = 1, thresh_num_max = 5;	// 1 5    or 5 5
	strcpy(file_fasta, argv[1]);
	strcpy(file_profile[0], argv[2]);	
	strcpy(file_profile[1], argv[3]);
	int len_motif[2];
	len_motif[0] = atoi(argv[4]);
	len_motif[1] = atoi(argv[5]);
	strcpy(file_table[0], argv[6]);
	strcpy(file_table[1], argv[7]);
	int shift_min = atoi(argv[8]); // minimal spacer length
	int shift_max = atoi(argv[9]); // upper bound of spacer length

	int height_permut = 100, size_min_permut = 200000, size_max_permut = 300000; //50000 150000 25
	//	int height_permut = 10, size_min_permut = 200, size_max_permut = 300; //50000 150000 25
	double pvalue = 0.0005, pvalue_mult = 1.5;// , dpvalue = 0.0000005; // 0.0005 1.5
	int mot_anchor = 0;// 0 = pwm from file >0 pwm from pre-computed database
	
	memset(name_fasta, '\0', sizeof(name_fasta));
	{
		k = 0;
		int len = strlen(file_fasta);
		for (j = 0; j < len; j++)
		{
			char cc = file_fasta[j];
			if (cc == '.')
			{
				name_fasta[k] = '\0';
				break;
			}
			name_fasta[k++] = cc;
		}
	}
	for (i = 0; i < 2; i++)
	{
		memset(name[i], '\0', sizeof(name[i]));
		int len = strlen(file_profile[i]);
		k = 0;
		char raz = '_';
		int count = 0, kmin=3, kmax = 4;
		for (j = 0; j < len; j++)
		{
			char cc = file_profile[i][j];
			if (cc == raz)
			{
				count++;
				if (count == kmin)continue;
				if (count == kmax)break;
			}	
			if (count < kmin)continue;
			if (cc == '.' || cc == '\0')
			{
				name[i][k] = '\0';
				break;
			}
			name[i][k++] = cc;
		}
		TransStrBack(name[i]);
	}
	strcpy(file_hist, name_fasta);
	strcat(file_hist, "_");
	strcat(file_hist, name[0]);
	strcat(file_hist, "_");
	strcat(file_hist, name[1]);
	strcat(file_hist, "_");
	strcat(file_hist, "out_hist");
	for (i = 0; i < 5; i++)
	{
		strcpy(file_pval[i], name_fasta);
		strcat(file_pval[i], "_");
		strcat(file_pval[i], name[0]);
		strcat(file_pval[i], "_");
		strcat(file_pval[i], name[1]);
		strcat(file_pval[i], "_");		
	}
	strcat(file_pval[0], "fisher_any_mot0");
	strcat(file_pval[1], "fisher_full_mot0");
	strcat(file_pval[2], "fisher_part_mot0");
	strcat(file_pval[3], "fisher_over_mot0");
	strcat(file_pval[4], "fisher_spac_mot0");
	
	strcpy(file_pval_table, name_fasta);
	strcat(file_pval_table, "_");
	strcat(file_pval_table, name[0]);
	strcat(file_pval_table, "_");
	strcat(file_pval_table, name[1]);
	strcat(file_pval_table, "_");
	strcat(file_pval_table, "out_pval");

	strcpy(file_hist_rand, name_fasta);
	strcat(file_hist_rand, "_");
	strcat(file_hist_rand, name[0]);
	strcat(file_hist_rand, "_");
	strcat(file_hist_rand, name[1]);
	strcat(file_hist_rand, "_");
	strcat(file_hist_rand, "out_hist_rand");
	
	int length_fasta_max = 0, nseq_real = 0;
	seq = NULL;
	int ftp = fasta_to_plain0(file_fasta, length_fasta_max, nseq_real);
	if (ftp == -1)
	{
		printf("File %s error 1st stage\n", file_fasta);
		return -1;
	}
	int *peak_len_real;
	peak_len_real = new int[nseq_real];
	if (peak_len_real == NULL) { puts("Out of memory..."); return -1; }

	seq = new char**[2];
	if (seq == NULL) { puts("Out of memory..."); return -1; }
	for (k = 0; k < 2; k++)
	{
		seq[k] = new char*[nseq_real];
		if (seq[k] == NULL) { puts("Out of memory..."); return -1; }
		for (i = 0; i < nseq_real; i++)
		{
			seq[k][i] = new char[length_fasta_max + 1];
			if (seq[k][i] == NULL) { puts("Out of memory..."); return -1; }
		}
	}
	ftp = fasta_to_plain1(file_fasta, length_fasta_max, nseq_real, seq, peak_len_real);
	if (ftp == -1)
	{
		printf("File %s error 2nd stage\n", file_fasta);
		return -1;
	}
	double thr[2][NUM_THR];

	profile real_one[2], rand_one[2], rand_hom_one;
	combi hist_obs_one, hist_exp_one;	

	//for real
	for (j = 0; j < 2; j++)
	{
		real_one[j].nseq = nseq_real;
		real_one[j].nam = 1;
		int ini = real_one[j].mem_in_nsit();
		if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	}

	//for permutation
	srand((unsigned)time(NULL));
	int nseq_rand = nseq_real * height_permut;
	if (nseq_rand < size_min_permut)height_permut = size_min_permut / nseq_real;
	if (nseq_rand > size_max_permut)height_permut = size_max_permut / nseq_real;
	nseq_rand = nseq_real * height_permut;
	double bonferroni_corr, bonferroni_corr_ap, bonferroni_corr_asy;
	for (j = 0; j < 2; j++)
	{
		rand_one[j].nseq = nseq_rand;
		rand_one[j].nam = 1;
		int ini = rand_one[j].mem_in_nsit();
		if (ini == -1) { puts("Not enough memory...\n"); return -1; }
	}
	rand_hom_one.nseq = nseq_rand;
	rand_hom_one.nam = 1;
	rand_hom_one.mot = 0;
	int ini = rand_hom_one.mem_in_nsit();
	if (ini == -1) { puts("Not enough memory...\n"); return -1; }

	int *peak_len_rand;
	peak_len_rand = new int[nseq_rand];
	if (peak_len_rand == NULL) { puts("Out of memory..."); return -1; }
	k = 0;
	for (i = 0; i < nseq_real; i++)
	{
		for (m = 0; m < height_permut; m++)peak_len_rand[k++] = peak_len_real[i];
	}
	int *thr_err_real, *thr_err_rand;
	thr_err_real = new int[nseq_real];
	if (thr_err_real == NULL) { puts("Out of memory..."); return -1; }
	thr_err_rand = new int[nseq_rand];
	if (thr_err_rand == NULL) { puts("Out of memory..."); return -1; }

	FILE *out_hist, *out_hist_rand;
	if ((out_hist = fopen(file_hist, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file_hist);
		return -1;
	}
	fclose(out_hist);
	if ((out_hist_rand = fopen(file_hist_rand, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file_hist_rand);
		return -1;
	}
	fclose(out_hist_rand);
	FILE *out_pval_table;
	if ((out_pval_table = fopen(file_pval_table, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file_pval_table);
		return -1;
	}
	{
		fprintf(out_pval_table, "# Motif");
		fprintf(out_pval_table, "\tMotif Name");
		fprintf(out_pval_table, "\tFull overlap, -Log10[P-value]");
		fprintf(out_pval_table, "\tPartial overlap,-Log10[P-value]");
		fprintf(out_pval_table, "\tOverlap, -Log10[P-value]");
		fprintf(out_pval_table, "\tSpacer, -Log10[P-value]");
		fprintf(out_pval_table, "\tAny, -Log10[P-value]");
		fprintf(out_pval_table, "\tSimilarity to Anchor, -Log10[P-value]");
		fprintf(out_pval_table, "\tSimilarity to Anchor, SSD");
		fprintf(out_pval_table, "\tSimilarity to Anchor, PCC\t");
		fprintf(out_pval_table, "Full overlap, Conservative Anchor, -Log10[P-value]\t");
		fprintf(out_pval_table, "Partial overlap, Conservative Anchor, -Log10[P-value]\t");
		fprintf(out_pval_table, "Overlap, Conservative Anchor, -Log10[P-value]\t");
		fprintf(out_pval_table, "Spacer, Conservative Anchor, -Log10[P-value]\t");
		fprintf(out_pval_table, "Any, Conservative Anchor, -Log10[P-value]\t");
		fprintf(out_pval_table, "Full overlap, Conservative Partner, -Log10[P-value]\t");
		fprintf(out_pval_table, "Partial overlap, Conservative Partner, -Log10[P-value]\t");
		fprintf(out_pval_table, "Overlap, Conservative Partner, -Log10[P-value]\t");
		fprintf(out_pval_table, "Spacer, Conservative Partner, -Log10[P-value]\t");
		fprintf(out_pval_table, "Any, Conservative Partner, -Log10[P-value]\t");
		fprintf(out_pval_table, "Full overlap, Asymmetry to Anchor+/Partner-, -Log10[P-value]\t");
		fprintf(out_pval_table, "Partial overlap, Asymmetry to Anchor+/Partner-, -Log10[P-value]\t");
		fprintf(out_pval_table, "Overlap, Asymmetry to Anchor+/Partner-, -Log10[P-value]\t");
		fprintf(out_pval_table, "Spacer, Asymmetry to Anchor+/Partner-, -Log10[P-value]\t");
		fprintf(out_pval_table, "Any, Asymmetry to Anchor+/Partner-, -Log10[P-value]\t");
		fprintf(out_pval_table, "Bonferroni_CE\tBonferroni_CE(AncPar)\tBonferroni_Asym\n");
		fclose(out_pval_table);
	}
	FILE *out_pval[5];
	double *thr_all, *fp_rate;	
	
	double thr_asy_min = 3.5, thr_asy_max = 5.5, dthr_asy = 0.2;
	int nthr_asy = (int)((thr_asy_max - thr_asy_min) / dthr_asy) + 2;
	asy_plot real_plot, rand_plot;
	real_plot.mem_in(nthr_asy, thr_asy_min, thr_asy_max, dthr_asy);
	rand_plot.mem_in(nthr_asy, thr_asy_min, thr_asy_max, dthr_asy);

	/*FILE *out_stat;
	if ((out_stat = fopen("rec_pos.txt", "wt")) == NULL)
	{
		printf("Input file can't be opened!\n");
		return -1;
	}
	fprintf(out_stat, "# Motif\tMotif Name\t# Threshold\tThreshold\t%% of peaks\tRec. peaks\tTotal peaks\tRate of hits\tRec. hits\tTotal positions\n");
	*/
	FILE *out_stat;
	char file_stat[120];
	memset(file_stat, '\0', sizeof(file_stat));
	strcpy(file_stat, name_fasta);
	strcat(file_stat, "_");
	strcat(file_stat, name[0]);
	strcat(file_stat, "_");
	strcat(file_stat, name[1]);	
	strcat(file_stat, "_stat.txt");
	if ((out_stat = fopen(file_stat, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n",file_stat);
		return -1;
	}
	fprintf(out_stat, "Motif Num\tMotif Name\t# Threshold\tThreshold\tFPR\t%% of peaks\tRec. peaks\tTotal peaks\tRate of hits\tRec. hits\tTotal positions with N\n");
	for (mot = 0; mot < 2; mot++)
	{
		//printf("Mot %d\n", mot);
		int n_thresh = 0; 
		FILE *in_tab;
		if ((in_tab = fopen(file_table[mot], "rt")) == NULL)
		{
			printf("Input file %s can't be opened!", file_table[mot]);
			return -1;
		}
		char d[100];		
		fgets(d, sizeof(d), in_tab);//header
		while (fgets(d, sizeof(d), in_tab) != NULL)
		{
			char c = d[0];
			char sep = '\t';
			if (c == '-' || isdigit(c))
			{
				char s[30];
				int test = UnderStol(d, 1, s, sep);
				if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
				n_thresh++;
				double fprx = atof(s);
				if (fprx > pvalue)break;				
			}
		}
		rewind(in_tab);
		fgets(d, sizeof(d), in_tab);//header		
		thr_all = new double[n_thresh];
		if (thr_all == NULL) { puts("Out of memory..."); return -1; }
		fp_rate = new double[n_thresh];
		if (fp_rate == NULL) { puts("Out of memory..."); return -1; }
		k = 0;
		while (fgets(d, sizeof(d), in_tab) != NULL)
		{
			char c = d[0];
			if (c == '-' || isdigit(c))
			{
				char s[30];
				int test = UnderStol(d, 1, s, '\t');
				if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
				thr_all[k] = atof(d);
				fp_rate[k] = atof(s);
				if (fp_rate[k] > pvalue)break;
				k++;
			}
		}
		fclose(in_tab);
		double fpr_select[NUM_THR];
		int index[NUM_THR];
		{
			int stfp = select_thresholds_from_pvalues(n_thresh, thr_all, fp_rate, pvalue, pvalue_mult, fpr_select, thr[mot], index);
			if (stfp == -1)
			{
				printf("Threshold selection of %d motif is wrong\n", mot);
				return -1;
			}
		}		
		////count sites per peaks				
		int rec = rec0(file_profile[mot], nseq_real,&real_one[mot]);
		if (rec == -1)
		{
			printf("Motif %d recognition 1st stage error\n", mot);
			return -1;
		}
		//memory allocation for all sites
		real_one[mot].clear_real();
		real_one[mot].mot = mot;
		real_one[mot].nam = NUM_THR;
		real_one[mot].count_sites();		
		//takes profiles
		rec = rec1(file_profile[mot], nseq_real, &real_one[mot]);
		if (rec == -1)
		{
			printf("Motif %d recognition 2nd stage error\n", mot);
			return -1;
		}
		int nhits[NUM_THR], npeaks[NUM_THR];
		for (j = 0; j < NUM_THR; j++)nhits[j] = npeaks[j] = 0;
		//count nsites for various thresholds
		for (i = 0; i < nseq_real; i++)
		{
			int jbest=10;
			for (k = 0; k < real_one[mot].nsit[i]; k++)
			{
				double sco = real_one[mot].sco[i][k];
				for (j = 0; j < NUM_THR; j++)
				{
					if (sco >= thr[mot][j])
					{
						if(j<jbest)jbest=j;
						real_one[mot].cel[i][k] = j;
						nhits[j]++;
						break;
					}
				}
			}
			if (jbest != 10)
			{
				for (j = NUM_THR - 1; j >= jbest; j--)npeaks[j]++;
			}
		}
		{
			int all_pos = 0;
			for (i = 0; i < nseq_real; i++)all_pos += (peak_len_real[i]- len_motif[mot]+1);
			for (i = 0; i < NUM_THR; i++)
			{
				fprintf(out_stat, "%d\t%s\t%d\t", mot, file_profile[mot],i+1);
				fprintf(out_stat, "%.12f\t%.13g\t",thr[mot][i], fpr_select[i]);
				fprintf(out_stat, "%.6f\t%d\t%d\t%d\t%g\t%d\n",100*(double)npeaks[i] / nseq_real, npeaks[i], nseq_real, nhits[i], (double)nhits[i] / all_pos, all_pos);
			}
		}
		//transform score to fprate
		for (i = 0; i < nseq_real; i++)
		{			
			for (k = 0; k < real_one[mot].nsit[i]; k++)
			{
				double sco = real_one[mot].sco[i][k];
				int inter = real_one[mot].cel[i][k];
				int inx2, inx1 = index[inter];
				if (inter == 0)inx2 = 0;
				else inx2 = index[inter - 1];
				double pv_sc = 0;
				for (j = inx1; j >= inx2; j--)
				{
					if (thr_all[j + 1] <= sco && thr_all[j] >= sco)
					{
						pv_sc = -log10(fp_rate[j]);
						break;
					}
				}
				real_one[mot].pv[i][k] = pv_sc;
			}
		}
		{
			char name_print[120];
			memset(name_print, '\0', sizeof(name_print));
			strcpy(name_print, name_fasta);
			strcat(name_print, "_");
			strcat(name_print, name[mot]);
			int fprint_pro = real_one[mot].fprintf_pro(name_print, thr[mot][NUM_THR - 1], xreal);
			if (fprint_pro == -1)
			{
				printf("Real print profile error, motif %d\n", mot);
				return -1;
			}
		}
		int rec_seq[NUM_THR], rec_pos[NUM_THR];
		for (i = 0; i < NUM_THR; i++)rec_seq[i] = rec_pos[i] = 0;
		for (i = 0; i < nseq_real; i++)
		{
			int inx[NUM_THR];
			for (j = 0; j < NUM_THR; j++)inx[j] = 0;
			for (k = 0; k < real_one[mot].nsit[i]; k++)
			{
				double sco = real_one[mot].sco[i][k];
				for (j = 0; j < NUM_THR; j++)
				{
					if (sco >= thr[mot][j])
					{
						rec_pos[j]++;
						inx[j]++;
						if (inx[j] == 1)rec_seq[j]++;
						break;
					}
				}
			}
		}
		delete[] fp_rate;
		delete[] thr_all;
		/*for (j = 0; j < NUM_THR; j++)
		{
			if (mot == 0)fprintf(out_stat, "Anchor");
			else fprintf(out_stat, "Partner %d", mot);
			fprintf(out_stat, "\t%s\t%d\t%f\t", name[mot], j + 1, thr[mot][j]);
			fprintf(out_stat, "%f\t%d\t%d\t", 100 * (double)rec_seq[j] / nseq_real, rec_seq[j], nseq_real);
			fprintf(out_stat, "%g\t%d\t%d\n", (double)rec_pos[j] / all_pos, rec_pos[j], all_pos);
		}*/
	}
	fclose(out_stat);		
	int len_anchor, len_partner;
	int mot_a = 0, mot_p;
	len_anchor = len_motif[mot_a];
	//anchor		
	int cop = real_one[mot_a].get_copy_rand(&rand_hom_one, height_permut);
	if (cop == -1)
	{
		printf("Rand Copy error Mot %d Thr One\n", mot);
		return -1;
	}

	char file_err[120];// = "throw_prediction.txt";	
	for (mot_p = 0; mot_p < 2; mot_p++)
	{
		len_partner = len_motif[mot_p];
		//one thresh  rand - &rand_hom_one,&rand_one[ap]   real - real_one[0],real_one[ap]
		//partner
		{
			int cop = real_one[mot_p].get_copy_rand(&rand_one[mot_p], height_permut);
			if (cop == -1)
			{
				printf("Rand Copy error Mot %d Thr One\n", mot);
				return -1;
			}
		}
		//		char file_throw_err[80], file_throw_err0[80];
		//	FILE *out_nsit_throw;
		{// one threshold
			//				int fprint_pro;
			/*int fprint_pro=rand_hom_one.fprintf_pro(name[mot_a],thr[mot_a][NUM_THR-1],"rand_do");
			if(fprint_pro==-1)
			{
			printf("Real print profile error, motif %d\n",mot);
			return -1;
			}
			fprint_pro=rand_one[mot_p].fprintf_pro(name[mot_p],thr[mot_p][NUM_THR-1],"rand_do");
			if(fprint_pro==-1)
			{
			printf("Real print profile error, motif %d\n",mot);
			return -1;
			}*/
			memset(file_err, '\0', sizeof(file_err));
			strcpy(file_err, name_fasta);
			strcat(file_err, "_");
			strcat(file_err, name[1-mot_p]);
			strcat(file_err, "_");
			strcat(file_err, name[mot_p]);
			strcat(file_err, "_throw_prediction.txt");
			FILE *out_err;
			if ((out_err = fopen(file_err, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!\n", file_err);
				return -1;
			}
			fclose(out_err);
			for (m = 0; m < nseq_real; m++)thr_err_real[m] = 0;
			rand_one[mot_p].mot = mot_p;
			int throwp = throw_predictions(peak_len_rand, &rand_hom_one, &rand_one[mot_p], len_anchor, len_partner, 0, thr_err_real, nseq_real, nseq_rand, seq[0], height_permut, file_err);
			if (throwp == -1)
			{
				printf("Throw Prediction error One - Anc 0 Par %d\n", mot);
				return -1;
			}
			/*	if(mot_p==1)
			{
			int fprint_pro=rand_hom_one.fprintf_pro(name[mot_a],thr[mot_a][NUM_THR-1],"rand_po");
			if(fprint_pro==-1)
			{
			printf("Real print profile error, motif %d\n",mot);
			return -1;
			}
			fprint_pro=rand_one[mot_p].fprintf_pro(name[mot_p],thr[mot_p][NUM_THR-1],"rand_po");
			if(fprint_pro==-1)
			{
			printf("Real print profile error, motif %d\n",mot);
			return -1;
			}
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
			real_plot.zero();
			rand_plot.zero();
			observed.ini();
			expected.ini();
			real_one[mot_p].mot = mot_p;
			k = 0;
			for (i = 0; i < nseq_real; i++)
			{
				for (j = 0; j < height_permut; j++)
				{
					thr_err_rand[k++] = thr_err_real[i];
				}
			}
			int nseq_two_sites_real = 0, nseq_two_sites_rand = 0;
			int proj = projoin(xrand, name[mot_p], rand_hom_one, rand_one[mot_p], shift_min, shift_max, len_anchor, len_partner, thr_err_rand, nseq_rand, seq, &expected, &hist_exp_one, peak_len_rand, &rand_plot, nseq_two_sites_rand);
			if (proj == -1)
			{
				printf("Projoin Rand error Anc 0 Par %d\n", mot);
				return -1;
			}
			char name_proj[120];
			memset(name_proj, '\0', sizeof(name_proj));
			strcpy(name_proj, name_fasta);
			strcat(name_proj, "_");
			strcat(name_proj, name[mot_a]);
			strcat(name_proj, "_");
			strcat(name_proj, name[mot_p]);
			proj = projoin(xreal, name_proj, real_one[0], real_one[mot_p], shift_min, shift_max, len_anchor, len_partner, thr_err_real, nseq_real, seq, &observed, &hist_obs_one, peak_len_real, &real_plot, nseq_two_sites_real);
			if (proj == -1)
			{
				printf("Projoin Real error Anc 0 Par %d\n", mot);
				return -1;
			}
			{
				double pv_standard = -log10(0.05);
				bonferroni_corr = (double)nseq_two_sites_real*nseq_two_sites_rand;
				bonferroni_corr *= 5;//potoki		
				bonferroni_corr_ap = bonferroni_corr;
				bonferroni_corr_asy = bonferroni_corr;
				bonferroni_corr *= NUM_THR;
				bonferroni_corr *= NUM_THR;
				bonferroni_corr_ap *= 2;
				bonferroni_corr = pv_standard + log10(bonferroni_corr);
				bonferroni_corr_ap = pv_standard + log10(bonferroni_corr_ap);
				bonferroni_corr_asy = pv_standard + log10(bonferroni_corr_asy);
			}
			char modew[] = "wt", modea[] = "at";
			char file_hist_one[80];
			strcpy(file_hist_one, file_hist);
			char buf[4];
			memset(buf, '\0', sizeof(buf));
			sprintf(buf, "%d", mot_p);
			strcat(file_hist_one, buf);
			hist_obs_one.fprintf_all(file_hist, mot_p, name[mot_p], len_anchor, len_partner, shift_max, modea);
			hist_obs_one.fprintf_all(file_hist_one, mot_p, name[mot_p], len_anchor, len_partner, shift_max, modew);
			hist_exp_one.fprintf_all(file_hist_rand, mot_p, name[mot_p], len_anchor, len_partner, shift_max, modea);
			real_plot.sum();
			rand_plot.sum();
			if (mot_p != 0)
			{
				char flow[5][8] = { "Full", "Partial", "Overlap", "Spacer", "Any" };
				FILE *out_plot[5];
				char file_plot[5][80];
				for (i = 0; i < 5; i++)
				{
					strcpy(file_plot[i], name_fasta);
					strcat(file_plot[i], "_");
					strcat(file_plot[i], name[0]);
					strcat(file_plot[i], "_");
					strcat(file_plot[i], name[1]);
					strcat(file_plot[i], "_");
					strcat(file_plot[i], "plot_");
					strcat(file_plot[i], flow[i]);
					strcat(file_plot[i], "_Anchor");
					strcat(file_plot[i], "_Partner");
					strcat(file_plot[i], buf);
					if ((out_plot[i] = fopen(file_plot[i], "wt")) == NULL)
					{
						printf("Input file %s can't be opened!\n", file_plot[i]);
						return -1;
					}
				}
				char rzd = ',';
				double val = real_plot.min;
				fprintf(out_plot[0], "%c<%.1f", rzd, val);
				for (i = 1; i < real_plot.n_karman; i++)
				{
					fprintf(out_plot[0], "%c%.1f", rzd, val);
					val += real_plot.inter;
					fprintf(out_plot[0], "..%.1f", val);
				}
				fprintf(out_plot[0], "%c>%.1f", rzd, val);
				fprintf(out_plot[0], "\n");
				val = real_plot.min;
				for (i = 0; i <= real_plot.n_karman; i++)
				{
					if (i == 0)fprintf(out_plot[0], "<%.1f", val);
					else
					{
						if (i == real_plot.n_karman)fprintf(out_plot[0], ">%.1f", val);
						else
						{
							fprintf(out_plot[0], "%.1f", val);
							val += real_plot.inter;
							fprintf(out_plot[0], "..%.1f", val);
						}
					}
					for (j = 0; j <= real_plot.n_karman; j++)
					{
						fprintf(out_plot[0], "%c", rzd);
						if (real_plot.total_full != 0 && rand_plot.total_full != 0)
						{
							double mill = 1000 * ((double)real_plot.full[i][j] / real_plot.total_full - (double)rand_plot.full[i][j] / rand_plot.total_full);
							if (fabs(mill) > 0.5)fprintf(out_plot[0], "%.f", mill);
						}
					}
					fprintf(out_plot[0], "\n");
				}
				fclose(out_plot[0]);
				val = real_plot.min;
				fprintf(out_plot[1], "%c<%.1f", rzd, val);
				for (i = 1; i < real_plot.n_karman; i++)
				{
					fprintf(out_plot[1], "%c%.1f", rzd, val);
					val += real_plot.inter;
					fprintf(out_plot[1], "..%.1f", val);
				}
				fprintf(out_plot[1], "%c>%.1f", rzd, val);
				fprintf(out_plot[1], "\n");
				val = real_plot.min;
				for (i = 0; i <= real_plot.n_karman; i++)
				{
					if (i == 0)fprintf(out_plot[1], "<%.1f", val);
					else
					{
						if (i == real_plot.n_karman)fprintf(out_plot[1], ">%.1f", val);
						else
						{
							fprintf(out_plot[1], "%.1f", val);
							val += real_plot.inter;
							fprintf(out_plot[1], "..%.1f", val);
						}
					}
					for (j = 0; j <= real_plot.n_karman; j++)
					{
						fprintf(out_plot[1], "%c", rzd);
						if (real_plot.total_partial != 0 && rand_plot.total_partial != 0)
						{
							double mill = 1000 * ((double)real_plot.partial[i][j] / real_plot.total_partial - (double)rand_plot.partial[i][j] / rand_plot.total_partial);
							if (fabs(mill) > 0.5)fprintf(out_plot[1], "%.f", mill);
						}
					}
					fprintf(out_plot[1], "\n");
				}
				fclose(out_plot[1]);
				val = real_plot.min;
				fprintf(out_plot[2], "%c<%.1f", rzd, val);
				for (i = 1; i < real_plot.n_karman; i++)
				{
					fprintf(out_plot[2], "%c%.1f", rzd, val);
					val += real_plot.inter;
					fprintf(out_plot[2], "..%.1f", val);
				}
				fprintf(out_plot[2], "%c>%.1f", rzd, val);
				fprintf(out_plot[2], "\n");
				val = real_plot.min;
				for (i = 0; i <= real_plot.n_karman; i++)
				{
					if (i == 0)fprintf(out_plot[2], "<%.1f", val);
					else
					{
						if (i == real_plot.n_karman)fprintf(out_plot[2], ">%.1f", val);
						else
						{
							fprintf(out_plot[2], "%.1f", val);
							val += real_plot.inter;
							fprintf(out_plot[2], "..%.1f", val);
						}
					}
					for (j = 0; j <= real_plot.n_karman; j++)
					{
						fprintf(out_plot[2], "%c", rzd);
						if (real_plot.total_overlap != 0 && rand_plot.total_overlap != 0)
						{
							double mill = 1000 * ((double)real_plot.overlap[i][j] / real_plot.total_overlap - (double)rand_plot.overlap[i][j] / rand_plot.total_overlap);
							if (fabs(mill) > 0.5)fprintf(out_plot[2], "%.f", mill);
						}
					}
					fprintf(out_plot[2], "\n");
				}
				fclose(out_plot[2]);
				val = real_plot.min;
				fprintf(out_plot[3], "%c<%.1f", rzd, val);
				for (i = 1; i < real_plot.n_karman; i++)
				{
					fprintf(out_plot[3], "%c%.1f", rzd, val);
					val += real_plot.inter;
					fprintf(out_plot[3], "..%.1f", val);
				}
				fprintf(out_plot[3], "%c>%.1f", rzd, val);
				fprintf(out_plot[3], "\n");
				val = real_plot.min;
				for (i = 0; i <= real_plot.n_karman; i++)
				{
					if (i == 0)fprintf(out_plot[3], "<%.1f", val);
					else
					{
						if (i == real_plot.n_karman)fprintf(out_plot[3], ">%.1f", val);
						else
						{
							fprintf(out_plot[3], "%.1f", val);
							val += real_plot.inter;
							fprintf(out_plot[3], "..%.1f", val);
						}
					}
					for (j = 0; j <= real_plot.n_karman; j++)
					{
						fprintf(out_plot[3], "%c", rzd);
						if (real_plot.total_spacer != 0 && rand_plot.total_spacer != 0)
						{
							double mill = 1000 * ((double)real_plot.spacer[i][j] / real_plot.total_spacer - (double)rand_plot.spacer[i][j] / rand_plot.total_spacer);
							if (fabs(mill) > 0.5)fprintf(out_plot[3], "%.f", mill);
						}
					}
					fprintf(out_plot[3], "\n");
				}
				fclose(out_plot[3]);
				val = real_plot.min;
				fprintf(out_plot[4], "%c<%.1f", rzd, val);
				for (i = 1; i < real_plot.n_karman; i++)
				{
					fprintf(out_plot[4], "%c%.1f", rzd, val);
					val += real_plot.inter;
					fprintf(out_plot[4], "..%.1f", val);
				}
				fprintf(out_plot[4], "%c>%.1f", rzd, val);
				fprintf(out_plot[4], "\n");
				val = real_plot.min;
				for (i = 0; i <= real_plot.n_karman; i++)
				{
					if (i == 0)fprintf(out_plot[4], "<%.1f", val);
					else
					{
						if (i == real_plot.n_karman)fprintf(out_plot[4], ">%.1f", val);
						else
						{
							fprintf(out_plot[4], "%.1f", val);
							val += real_plot.inter;
							fprintf(out_plot[4], "..%.1f", val);
						}
					}
					for (j = 0; j <= real_plot.n_karman; j++)
					{
						fprintf(out_plot[4], "%c", rzd);
						if (real_plot.total_any != 0 && rand_plot.total_any != 0)
						{
							double mill = 1000 * ((double)real_plot.any[i][j] / real_plot.total_any - (double)rand_plot.any[i][j] / rand_plot.total_any);
							if (fabs(mill) > 0.5)fprintf(out_plot[4], "%.f", mill);
						}
					}
					fprintf(out_plot[4], "\n");
				}
				fclose(out_plot[4]);
			}
		}// one threshold
		//many thresholds
		for (i = 0; i < NUM_THR; i++)for (j = 0; j < NUM_THR; j++)pvalue_a[i][j] = pvalue_f[i][j] = pvalue_p[i][j] = pvalue_o[i][j] = pvalue_s[i][j] = 1;
		for (j = 0; j < NUM_THR; j++)
		{
			for (k = 0; k < NUM_THR; k++)
			{
				//printf("Mot %d J %d K %d\n",mot,j+1,k+1);								
				int fisher = fisher_exact_test(observed.cell[j][k].any, observed.cell[j][k].two_sites, expected.cell[j][k].any, expected.cell[j][k].two_sites, pvalue_a[j][k], 0);
				if (fisher == -1)
				{
					printf("Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher = fisher_exact_test(observed.cell[j][k].full, observed.cell[j][k].two_sites, expected.cell[j][k].full, expected.cell[j][k].two_sites, pvalue_f[j][k], 0);
				if (fisher == -1)
				{
					printf("Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher = fisher_exact_test(observed.cell[j][k].partial, observed.cell[j][k].two_sites, expected.cell[j][k].partial, expected.cell[j][k].two_sites, pvalue_p[j][k], 0);
				if (fisher == -1)
				{
					printf("Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher = fisher_exact_test(observed.cell[j][k].overlap, observed.cell[j][k].two_sites, expected.cell[j][k].overlap, expected.cell[j][k].two_sites, pvalue_o[j][k], 0);
				if (fisher == -1)
				{
					printf("Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher = fisher_exact_test(observed.cell[j][k].spacer, observed.cell[j][k].two_sites, expected.cell[j][k].spacer, expected.cell[j][k].two_sites, pvalue_s[j][k], 0);
				if (fisher == -1)
				{
					printf("Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
			}
		}
		{
			pv_any.ini();
			pv_full.ini();
			pv_partial.ini();
			pv_overlap.ini();
			pv_spacer.ini();
			//any
			int fisher;
			fisher = fisher_exact_test(observed.anc.any, observed.anc.two_sites, expected.anc.any, expected.anc.two_sites, pv_any.anchor, 0);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.par.any, observed.par.two_sites, expected.par.any, expected.par.two_sites, pv_any.partner, 0);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.anc_sit.any, observed.sit.any, expected.anc_sit.any, expected.sit.any, pv_any.anc_par, 1);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			//full
			fisher = fisher_exact_test(observed.anc.full, observed.anc.two_sites, expected.anc.full, expected.anc.two_sites, pv_full.anchor, 0);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.par.full, observed.par.two_sites, expected.par.full, expected.par.two_sites, pv_full.partner, 0);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.anc_sit.full, observed.sit.full, expected.anc_sit.full, expected.sit.full, pv_full.anc_par, 1);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			//partial
			fisher = fisher_exact_test(observed.anc.partial, observed.anc.two_sites, expected.anc.partial, expected.anc.two_sites, pv_partial.anchor, 0);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.par.partial, observed.par.two_sites, expected.par.partial, expected.par.two_sites, pv_partial.partner, 0);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.anc_sit.partial, observed.sit.partial, expected.anc_sit.partial, expected.sit.partial, pv_partial.anc_par, 1);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			//overlap
			fisher = fisher_exact_test(observed.anc.overlap, observed.anc.two_sites, expected.anc.overlap, expected.anc.two_sites, pv_overlap.anchor, 0);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.par.overlap, observed.par.two_sites, expected.par.overlap, expected.par.two_sites, pv_overlap.partner, 0);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.anc_sit.overlap, observed.sit.overlap, expected.anc_sit.overlap, expected.sit.overlap, pv_overlap.anc_par, 1);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			//spacer
			fisher = fisher_exact_test(observed.anc.spacer, observed.anc.two_sites, expected.anc.spacer, expected.anc.two_sites, pv_spacer.anchor, 0);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.par.spacer, observed.par.two_sites, expected.par.spacer, expected.par.two_sites, pv_spacer.partner, 0);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.anc_sit.spacer, observed.sit.spacer, expected.anc_sit.spacer, expected.sit.spacer, pv_spacer.anc_par, 1);
			if (fisher == -1)
			{
				printf("Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
		}
		for (i = 0; i < 5; i++)
		{
			char file_pval0[80];
			strcpy(file_pval0, file_pval[i]);
			char buf[10];
			sprintf(buf, "%d", mot_p);
			strcat(file_pval0, buf);
			if ((out_pval[i] = fopen(file_pval0, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!\n", file_pval0);
				return -1;
			}
			fprintf(out_pval[i], "Anchor Thr\tPartner Thr\t\tReal CE+\tReal Total\tRand CE+\tRand Total\tFold\tP-value\n");
		}
		for (j = 0; j < NUM_THR; j++)
		{
			for (k = 0; k < NUM_THR; k++)
			{
				for (i = 0; i < 5; i++)fprintf(out_pval[i], "A %d\tP %d\t\t", j + 1, k + 1);
				fprintf(out_pval[0], "%d\t%d\t%d\t%d\t\t%g\n", observed.cell[j][k].any, observed.cell[j][k].two_sites, expected.cell[j][k].any, expected.cell[j][k].two_sites, pvalue_a[j][k]);
				fprintf(out_pval[1], "%d\t%d\t%d\t%d\t\t%g\n", observed.cell[j][k].full, observed.cell[j][k].two_sites, expected.cell[j][k].full, expected.cell[j][k].two_sites, pvalue_f[j][k]);
				fprintf(out_pval[2], "%d\t%d\t%d\t%d\t\t%g\n", observed.cell[j][k].partial, observed.cell[j][k].two_sites, expected.cell[j][k].partial, expected.cell[j][k].two_sites, pvalue_p[j][k]);
				fprintf(out_pval[3], "%d\t%d\t%d\t%d\t\t%g\n", observed.cell[j][k].overlap, observed.cell[j][k].two_sites, expected.cell[j][k].overlap, expected.cell[j][k].two_sites, pvalue_o[j][k]);
				fprintf(out_pval[4], "%d\t%d\t%d\t%d\t\t%g\n", observed.cell[j][k].spacer, observed.cell[j][k].two_sites, expected.cell[j][k].spacer, expected.cell[j][k].two_sites, pvalue_s[j][k]);
			}
		}
		double rat_a, rat_p;
		for (i = 0; i < 5; i++)fprintf(out_pval[i], "\n");
		//any vs rand
		fprintf(out_pval[0], "Anchor\t\t\t%d\t%d\t%d\t%d\t\t%g\n", observed.anc.any, observed.anc.two_sites, expected.anc.any, expected.anc.two_sites, pv_any.anchor);
		fprintf(out_pval[0], "Partner\t\t\t%d\t%d\t%d\t%d\t\t%g\n", observed.par.any, observed.par.two_sites, expected.par.any, expected.par.two_sites, pv_any.partner);
		//any vs real
		rat_a = (double)observed.anc_sit.any / observed.sit.any;
		rat_p = (double)expected.anc_sit.any / expected.sit.any;
		if (rat_p == 0)pv_any.fold_anc_par = 1000;
		else pv_any.fold_anc_par = rat_a / rat_p;
		fprintf(out_pval[0], "Anchor_Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc_sit.any, observed.sit.any, expected.anc_sit.any, expected.sit.any, pv_any.fold_anc_par, pv_any.anc_par);
		//full vs rand		
		fprintf(out_pval[1], "Anchor\t\t\t%d\t%d\t%d\t%d\t\t%g\n", observed.anc.full, observed.anc.two_sites, expected.anc.full, expected.anc.two_sites, pv_full.anchor);
		fprintf(out_pval[1], "Partner\t\t\t%d\t%d\t%d\t%d\t\t%g\n", observed.par.full, observed.par.two_sites, expected.par.full, expected.par.two_sites, pv_full.partner);
		//full vs real
		rat_a = (double)observed.anc_sit.full / (observed.anc_sit.full + observed.par_sit.full);
		rat_p = (double)expected.anc_sit.full / (expected.anc_sit.full + expected.par_sit.full);
		if (rat_p == 0)pv_full.fold_anc_par = 1000;
		else pv_full.fold_anc_par = rat_a / rat_p;
		fprintf(out_pval[1], "Anchor_Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc_sit.full, observed.sit.full, expected.anc_sit.full, expected.sit.full, pv_full.fold_anc_par, pv_full.anc_par);
		//partial vs rand
		fprintf(out_pval[2], "Anchor\t\t\t%d\t%d\t%d\t%d\t\t%g\n", observed.anc.partial, observed.anc.two_sites, expected.anc.partial, expected.anc.two_sites, pv_partial.anchor);
		fprintf(out_pval[2], "Partner\t\t\t%d\t%d\t%d\t%d\t\t%g\n", observed.par.partial, observed.par.two_sites, expected.par.partial, expected.par.two_sites, pv_partial.partner);
		//partial vs real
		rat_a = (double)observed.anc_sit.partial / observed.sit.partial;
		rat_p = (double)expected.anc_sit.partial / expected.sit.partial;
		if (rat_p == 0)pv_partial.fold_anc_par = 1000;
		else pv_partial.fold_anc_par = rat_a / rat_p;
		fprintf(out_pval[2], "Anchor_Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc_sit.partial, observed.sit.partial, expected.anc_sit.partial, expected.sit.partial, pv_partial.fold_anc_par, pv_partial.anc_par);
		//overlap vs rand
		fprintf(out_pval[3], "Anchor\t\t\t%d\t%d\t%d\t%d\t\t%g\n", observed.anc.overlap, observed.anc.two_sites, expected.anc.overlap, expected.anc.two_sites, pv_overlap.anchor);
		fprintf(out_pval[3], "Partner\t\t\t%d\t%d\t%d\t%d\t\t%g\n", observed.par.overlap, observed.par.two_sites, expected.par.overlap, expected.par.two_sites, pv_overlap.partner);
		//overlap vs real
		rat_a = (double)observed.anc_sit.overlap / observed.sit.overlap;
		rat_p = (double)expected.anc_sit.overlap / expected.sit.overlap;
		if (rat_p == 0)pv_overlap.fold_anc_par = 1000;
		else pv_overlap.fold_anc_par = rat_a / rat_p;
		fprintf(out_pval[3], "Anchor_Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc_sit.overlap, observed.sit.overlap, expected.anc_sit.overlap, expected.sit.overlap, pv_overlap.fold_anc_par, pv_overlap.anc_par);
		//spacer vs rand
		fprintf(out_pval[4], "Anchor\t\t\t%d\t%d\t%d\t%d\t\t%g\n", observed.anc.spacer, observed.anc.two_sites, expected.anc.spacer, expected.anc.two_sites, pv_spacer.anchor);
		fprintf(out_pval[4], "Partner\t\t\t%d\t%d\t%d\t%d\t\t%g\n", observed.par.spacer, observed.par.two_sites, expected.par.spacer, expected.par.two_sites, pv_spacer.partner);
		//spacer vs real
		rat_a = (double)observed.anc_sit.spacer / observed.sit.spacer;
		rat_p = (double)expected.anc_sit.spacer / expected.sit.spacer;
		if (rat_p == 0)pv_spacer.fold_anc_par = 1000;
		else pv_spacer.fold_anc_par = rat_a / rat_p;
		fprintf(out_pval[4], "Anchor_Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc_sit.spacer, observed.sit.spacer, expected.anc_sit.spacer, expected.sit.spacer, pv_spacer.fold_anc_par, pv_spacer.anc_par);
		for (i = 0; i < 5; i++)fclose(out_pval[i]);

		double pval_tot_min[5] = { 0, 0, 0, 0, 0 };
		double limit = 300;
		double pv_limit = 1E-300;
		double lgpv[5][NUM_THR][NUM_THR];
		for (j = 0; j < NUM_THR; j++)
		{
			for (k = 0; k < NUM_THR; k++)
			{
				for (i = 0; i < 5; i++)lgpv[i][k][j] = 0;
				if (pvalue_a[k][j] <= pv_limit)lgpv[0][k][j] = limit;
				else lgpv[0][k][j] = -log10(pvalue_a[k][j]);
				if (pvalue_f[k][j] <= pv_limit)lgpv[1][k][j] = limit;
				else lgpv[1][k][j] = -log10(pvalue_f[k][j]);
				if (pvalue_p[k][j] <= pv_limit)lgpv[2][k][j] = limit;
				else lgpv[2][k][j] = -log10(pvalue_p[k][j]);
				if (pvalue_o[k][j] <= pv_limit)lgpv[3][k][j] = limit;
				else lgpv[3][k][j] = -log10(pvalue_o[k][j]);
				if (pvalue_s[k][j] <= pv_limit)lgpv[4][k][j] = limit;
				else lgpv[4][k][j] = -log10(pvalue_s[k][j]);
			}
		}
		for (i = 0; i < 5; i++)
		{
			int limit_found = 0;
			int mnoj = NUM_THR * NUM_THR;
			for (j = 0; j < NUM_THR; j++)
			{
				for (k = 0; k < NUM_THR; k++)
				{
					double val = lgpv[i][k][j];
					if (val == limit)
					{
						limit_found = 1;
						break;
					}
				}
				if (limit_found == 1)break;
			}
			if (limit_found == 1)
			{
				pval_tot_min[i] = limit;
			}
			else
			{
				for (j = 0; j < NUM_THR; j++)
				{
					for (k = 0; k < NUM_THR; k++)
					{
						double val = lgpv[i][k][j];
						if (val > pval_tot_min[i])pval_tot_min[i] = val;
					}
				}
				//pval_tot_min[i]=pow(10,-pval_tot_min[i]);
			}
		}
		{
			pv_any.anc_par = -log10(pv_any.anc_par);
			pv_full.anc_par = -log10(pv_full.anc_par);
			pv_partial.anc_par = -log10(pv_partial.anc_par);
			pv_overlap.anc_par = -log10(pv_overlap.anc_par);
			pv_spacer.anc_par = -log10(pv_spacer.anc_par);
			if (pv_any.anc_par > bonferroni_corr_asy)
			{
				if (pv_any.fold_anc_par < 1)pv_any.anc_par *= -1;
			}
			else pv_any.anc_par = 0;
			if (pv_full.anc_par > bonferroni_corr_asy)
			{
				if (pv_full.fold_anc_par < 1)pv_full.anc_par *= -1;
			}
			else pv_full.anc_par = 0;
			if (pv_partial.anc_par > bonferroni_corr_asy)
			{
				if (pv_partial.fold_anc_par < 1)pv_partial.anc_par *= -1;
			}
			else pv_partial.anc_par = 0;
			if (pv_overlap.anc_par > bonferroni_corr_asy)
			{
				if (pv_overlap.fold_anc_par < 1)pv_overlap.anc_par *= -1;
			}
			else pv_overlap.anc_par = 0;
			if (pv_spacer.anc_par > bonferroni_corr_asy)
			{
				if (pv_spacer.fold_anc_par < 1)pv_spacer.anc_par *= -1;
			}
			else pv_spacer.anc_par = 0;
		}
		if ((out_pval_table = fopen(file_pval_table, "at")) == NULL)
		{

			printf("Input file %s can't be opened!\n", file_pval_table);
			return -1;
		}
		if (mot_p == 0)fprintf(out_pval_table, "Anchor");
		else fprintf(out_pval_table, "Partner %d", mot_p);
		fprintf(out_pval_table, "\t%s", name[mot_p]);
		for (i = 1; i < 5; i++)
		{
			if (pval_tot_min[i] > bonferroni_corr)fprintf(out_pval_table, "\t%.2f", pval_tot_min[i]);
			else fprintf(out_pval_table, "\t0");
		}
		if (pval_tot_min[0] > bonferroni_corr)fprintf(out_pval_table, "\t%.2f", pval_tot_min[0]);
		else fprintf(out_pval_table, "\t0");				
		if (mot_a != mot_p)
		{
			pv_full.anchor = -log10(pv_full.anchor);
			pv_full.partner = -log10(pv_full.partner);
			pv_partial.anchor = -log10(pv_partial.anchor);
			pv_partial.partner = -log10(pv_partial.partner);
			pv_overlap.anchor = -log10(pv_overlap.anchor);
			pv_overlap.partner = -log10(pv_overlap.partner);
			pv_spacer.anchor = -log10(pv_spacer.anchor);
			pv_spacer.partner = -log10(pv_spacer.partner);
			pv_any.anchor = -log10(pv_any.anchor);
			pv_any.partner = -log10(pv_any.partner);
			//pvalue_similarity_tot = pfm_similarity(&matrix[mot_a], &matrix[mot_p], s_granul, s_overlap_min, s_ncycle_small, s_ncycle_large, pval_sim);
			//fprintf(out_pval_table,"\t%+.2f\t%+.2f\t%+.2f",pv_full.anc_par,pv_overlap.anc_par,pv_any.anc_par);
			fprintf(out_pval_table, "\t\t\t");
			if (pv_full.anchor > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_full.anchor);
			else fprintf(out_pval_table, "\t0");
			if (pv_partial.anchor > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_partial.anchor);
			else fprintf(out_pval_table, "\t0");
			if (pv_overlap.anchor > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_overlap.anchor);
			else fprintf(out_pval_table, "\t0");
			if (pv_spacer.anchor > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_spacer.anchor);
			else fprintf(out_pval_table, "\t0");
			if (pv_any.anchor > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_any.anchor);
			else fprintf(out_pval_table, "\t0");
			if (pv_full.partner > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_full.partner);
			else fprintf(out_pval_table, "\t0");
			if (pv_partial.partner > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_partial.partner);
			else fprintf(out_pval_table, "\t0");
			if (pv_overlap.partner > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_overlap.partner);
			else fprintf(out_pval_table, "\t0");
			if (pv_spacer.partner > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_spacer.partner);
			else fprintf(out_pval_table, "\t0");
			if (pv_any.partner > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_any.partner);
			else fprintf(out_pval_table, "\t0");
			if (pv_full.anc_par != 0)fprintf(out_pval_table, "\t%+.2f", pv_full.anc_par);
			else fprintf(out_pval_table, "\t0");
			if (pv_partial.anc_par != 0)fprintf(out_pval_table, "\t%+.2f", pv_partial.anc_par);
			else fprintf(out_pval_table, "\t0");
			if (pv_overlap.anc_par != 0)fprintf(out_pval_table, "\t%+.2f", pv_overlap.anc_par);
			else fprintf(out_pval_table, "\t0");
			if (pv_spacer.anc_par != 0)fprintf(out_pval_table, "\t%+.2f", pv_spacer.anc_par);
			else fprintf(out_pval_table, "\t0");
			if (pv_any.anc_par != 0)fprintf(out_pval_table, "\t%+.2f", pv_any.anc_par);
			else fprintf(out_pval_table, "\t0");
			fprintf(out_pval_table, "\t%.2f", bonferroni_corr);
			fprintf(out_pval_table, "\t%.2f", bonferroni_corr_ap);
			fprintf(out_pval_table, "\t%.2f", bonferroni_corr_asy);
		}
		else
		{
			char slash = '/';
			for (i = 0; i < 18; i++)fprintf(out_pval_table, "\tn%ca", slash);
			fprintf(out_pval_table, "\t%.2f\t\t", bonferroni_corr);
		}
		fprintf(out_pval_table, "\n");
		fclose(out_pval_table);
		rand_one[mot_p].mem_out_sta();
		rand_one[mot_p].mem_out_cep();
		rand_one[mot_p].mem_out_cel();
		rand_one[mot_p].mem_out_pv();
		for (i = 0; i < nseq_rand; i++)rand_one[mot_p].nsit[i] = 0;
	}
	real_plot.mem_out();
	rand_plot.mem_out();
	rand_hom_one.mem_out_sta();
	rand_hom_one.mem_out_cep();
	rand_hom_one.mem_out_cel();
	rand_hom_one.mem_out_pv();
	rand_hom_one.mem_out_nsit();
	for (i = 0; i < 2; ++i)
	{
		real_one[i].mem_out_sta();
		real_one[i].mem_out_cep();
		real_one[i].mem_out_cel();
		real_one[i].mem_out_sco();
		real_one[i].mem_out_pv();
		real_one[i].mem_out_nsit();
		rand_one[i].mem_out_nsit();
	}
	delete[] thr_err_real;
	delete[] thr_err_rand;
	delete[] peak_len_real;
	delete[] peak_len_rand;
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < nseq_real; i++)
		{
			delete[] seq[k][i];
		}
		delete[] seq[k];
	}
	delete[] seq;
	return 0;
}
