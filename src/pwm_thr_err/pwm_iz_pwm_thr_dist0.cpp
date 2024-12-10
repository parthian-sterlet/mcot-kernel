#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 5050
#define PWMLEN 50 // max length of motifs
#define OLIGNUM 4// di 16 mono 4
double pfm[SEQLEN][OLIGNUM];
double pwm[SEQLEN][OLIGNUM];

char *TransStr(char *d)
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i<lens; i++)
	{
		c = int(d[i]);
		if (c<97) d[i] = char(c + 32);
		//else break;
	}
	return(d);
}
char *TransStrBack(char *d)//a->A
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i<lens; i++)
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
	for (i = 0; i<lens; i++)
	{
		if (str[i] != c)str[size++] = str[i];
	}
	str[size] = '\0';
}
int IdeLet(char c)
{
	int ret;
	switch (c){
	case 'n': return -1000;
	case 'a': ret = 0; break;
	case 't': ret = 1; break;
	case 'g': ret = 2; break;
	case 'c': ret = 3; break;
	default: ret = -1;
	}
	return(ret);
}
/*
void GetSost(char *d, int word, int *sost, char *letter)
{
int i, j, k, i_sost, let;
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
sost[i_sost]++;
}
}
*/
int GetSostPro(char *d, int word, int *sost, char *letter)
{
	int i, j, k, i_sost, let;
	int ten[6] = { 1, 4, 16, 64, 256, 1024 };
	int lens = strlen(d);
	int size = 1;
	for (k = 0; k<word; k++)size *= 4;
	for (i = 0; i<size; i++)sost[i] = 0;
	for (i = 0; i<lens - word + 1; i++)
	{
		i_sost = 0;
		let = -1;
		for (j = word - 1; j >= 0; j--)
		{
			for (k = 0; k<4; k++)
			{
				if (d[i + j] == letter[k]){ let = k; break; }
			}
			if (let == -1)return -1;
			i_sost += ten[word - 1 - j] * let;
		}
		sost[i] = i_sost;
	}
	return 0;
}
void PWMScore(double &min, double &raz, int len1, int olig)
{
	int i, j;
	for (i = 0; i<len1; i++)
	{
		double pwmmin = 100;
		double pwmmax = -100;
		for (j = 0; j<olig; j++)
		{
			if (pwm[i][j]<pwmmin)pwmmin = pwm[i][j];
			if (pwm[i][j]>pwmmax)pwmmax = pwm[i][j];
		}
		raz += pwmmax;
		min += pwmmin;
	}
	raz -= min;
}
int ComplStr(char *d)
{
	char d1[SEQLEN];
	int i, len;
	len = strlen(d);
	strcpy(d1, d);
	//	memset(d,0,sizeof(d));
	for (i = 0; i<len; i++)
	{
		switch (d1[len - i - 1])
		{
		case 'a':{d[i] = 't'; break; }
		case 't':{d[i] = 'a'; break; }
		case 'c':{d[i] = 'g'; break; }
		case 'g':{d[i] = 'c'; break; }
		case 'A':{d[i] = 'T'; break; }
		case 'T':{d[i] = 'A'; break; }
		case 'C':{d[i] = 'G'; break; }
		case 'G':{d[i] = 'C'; break; }
		case 'N':{d[i] = 'N'; break; }
		case 'n':{d[i] = 'n'; break; }
		default: d[i] = 'n';
		}
	}	
	return 1;
}
int CheckStr(char *d)
{
	int i, len, ret, size;
	ret = size = 0;
	len = strlen(d);
	for (i = 0; i<len; i++)
	{
		if (strchr("atgcn", (int)d[i]) != NULL){ continue; }
		else { ret++; }
	}
	return(ret);
}
int ConvertSym(int c)
{
	char four[] = "atgc";
	char di[6][3] = { "ag", "tc", "at", "ac", "gt", "gc" };
	char tri[4][4] = { "agt", "agc", "gtc", "tca" };

	if ((c>64) && (c<91)) c = char(c + 32);

	switch (c)
	{
	case 'r':{
		c = (int)di[1][rand() % 2];
		break; }
	case 'y':{
		c = (int)di[1][rand() % 2];
		break; }
	case 'w':{
		c = (int)di[2][rand() % 2];
		break; }
	case 'm':{
		c = (int)di[3][rand() % 2];
		break; }
	case 'k':{
		c = (int)di[4][rand() % 2];
		break; }
	case 's':{
		c = (int)di[5][rand() % 2];
		break; }
	case 'd':{
		c = (int)tri[0][rand() % 3];
		break; }
	case 'v':{
		c = (int)tri[1][rand() % 3];
		break; }
	case 'b':{
		c = (int)tri[2][rand() % 3];
		break; }
	case 'h':{
		c = (int)tri[3][rand() % 3];
		break; }
	case 'n':{
		c = (int)four[rand() % 4];
		break; }
	default:{
		c = (int)'n';
		return -1; }
	}
	return 1;
}

int ReadSeq(char *file, int &n, int &len1, int &all_pos)
{
	char head[1000];
	int fl = 0, len;
	char symbol;
	int c;
	//	char cyfr[]="0123456789";
	FILE  *in;
	len1 = len = n = 0;
	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		return -1;
	}
	symbol = fgetc(in);
	rewind(in);
	while ((c = fgetc(in)) != -1)
	{
		if ((char)c == symbol)
		{
			if (fgets(head, sizeof(head), in) == NULL)return -1;
			if (len>len1)len1 = len;
			len = 0;
			n++;
			continue;
		}
		if (strchr("\t\n ", c) != NULL)continue;
		if (strchr("ATGCNatgc", c) != NULL)all_pos++;
		if (strchr("ATGCNatgcn", c) != NULL)
		{
			len++;						
			continue;
		}
		else
		{
			printf("Unusual base: sequence N %d, letter position %d\n symbol %c\n%s", n, len, c, head);
			exit(1);
		}
	}
	if (len>len1)len1 = len;
	fclose(in);
	return 1;
}
void ReplaceChar(char *str, char c1, char c2)
{
	int i, len = strlen(str);
	for (i = 0; i<len; i++)
	{
		if (str[i] == c1) str[i] = c2;
	}
}
int StrNStr(char *str, char c, int n)
{
	int i, len = strlen(str);
	int k = 1;
	for (i = 0; i<len; i++)
	{
		if (str[i] == c)
		{
			if (k == n)return i;
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
void Mix(double *a, double *b)
{
	double buf = *a;
	*a = *b;
	*b = buf;
}
int main(int argc, char *argv[])
{
	int i, j, k, n;
	char head[1000], d1[SEQLEN], file_out_distt[300], file_out_distb[300], letter[5], binary_ext[] = ".binary", binary_mode[3];
	//char file_out_cpp_arr[80], file_out_cpp_struct_many[80], name[20];
	char file_pfm[300], file_pwm[300], file_seq[300], file_sta[300];
	FILE *in,* in_pfm, *in_pwm, * out_distt, * out_distb;// , * out_cpp_arr;
	
	if (argc != 12)
	{
		printf("%s 1 pfm 2pwm 3file_profile_fasta 4file out_dist_text 5file out_dist_binary 6double pvalue_large 7double bin 8file_sta 9double thr_best_pval 10int 1/0 check/don't check Bad&Good PWM 11char wb OR ab mode for output binary", argv[0]);
		return -1;
	}	
	strcpy(file_pfm, argv[1]);
	strcpy(file_pwm, argv[2]);
	strcpy(file_seq, argv[3]);
	strcpy(letter, "acgt");//atgc
	letter[4] = '\0';
	strcpy(file_out_distt, argv[4]);
	strcpy(file_out_distb, argv[5]);
	double pvalue_large = atof(argv[6]);
	double log_pvalue_large = -log10(pvalue_large);
	double bin =atof(argv[7]);	
	strcpy(file_sta, argv[8]);
	double thr_best_pval = atof(argv[9]);//2E-5=0.00002
	double log_thr_best_pval = -log10(thr_best_pval);	
	int check_bad = atoi(argv[10]); // 0 all results are written, 1 results passing p-value threshold are written only
	strcpy(binary_mode, argv[11]);// wb / ab create new binary file for output data or add output data to this file
	int nseq = 0;
	int len1 = 0;
	int word;
	int olig;
	if ((in_pwm = fopen(file_pwm, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file_pwm);
		return -1;
	}
	int test;
	char val[50];
	char sep = '\t';
	int lenp = 0;// dlina matricy
	fgets(d1, sizeof(d1), in_pwm);
	if (fgets(d1, sizeof(d1), in_pwm) != NULL)
	{
		olig = 0;
		for (i = 0;; i++)
		{
			if (d1[i] == '.')olig++;
			if (d1[i] == '\n')
			{
				break;
			}
		}
		DelChar(d1, '\n');
		char razdelitel;
		if (strchr(d1, sep) != NULL)
		{
			razdelitel = sep;
		}
		else
		{
			if (strchr(d1, ' ') != NULL)razdelitel = ' ';
			else
			{
				printf("Unknown razdelitel v stroke\n%s", d1);
				return -1;
			}
		}
		if (olig != 4)
		{
			printf("Wrong no. of columns in the input matrix!\n");
			return -1;
		}
		if (olig == 4)word = 1;
		if (olig == 16)word = 2;
		j = 0;//pozicii
		int fl = 0;
		do
		{
			for (i = 0; i<olig; i++)
			{
				test = UnderStol(d1, i, val, sizeof(val), razdelitel);
				if (test == -1) { printf("Wrong format %s\n", d1); exit(1); }
				pwm[j][i] = atof(val);
			}
			j++;
			if (fgets(d1, sizeof(d1), in_pwm) == NULL)break;
			if (*d1 == '\n')fl = 1;
		} 
		while (fl == 0);
		lenp = j;
	}
	fclose(in_pwm);
	int shift_col, mtype;
	if ((in_pfm = fopen(file_pfm, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file_pfm);
		return -1;
	}
	int alfabet = 0;
	fgets(head, sizeof(head), in_pfm);
	fgets(head, sizeof(head), in_pfm);
	if (strchr("ATGC", head[0]) != NULL)
	{// jaspar
		nseq = 1;
		mtype = 1;
		shift_col = 1;
		alfabet++;
		while (fgets(head, sizeof(head), in_pfm) != NULL)
		{
			if (strchr("ATGC", head[0]) != NULL)alfabet++;
		}
		if (alfabet != 4 && alfabet != 16)
		{
			printf("Reading error in Jaspar matrix file %s !\n", file_pfm);
			return -1;
		}
	}
	else
	{//homer, cis-bp
		nseq = 1000000000;
		mtype = 0;
		//shift_col=0;
	}
	if (mtype == 0)
	{
		DelChar(head, '\n');
		DelChar(head, '\r');
		int headlen = strlen(head);
		if (head[headlen - 1] == sep)
		{
			head[headlen - 1] = '\0';
			headlen--;
		}
		int counttab = 0;
		for (i = 0; i < headlen; i++)
		{
			if (head[i] == '\t')counttab++;
		}
		if (counttab == 4)// || counttab == 16
		{
			shift_col = 1;
			alfabet = counttab;
		}
		else
		{
			if (counttab == 3)// || counttab == 15
			{
				shift_col = 0;
				alfabet = counttab + 1;
			}
			else
			{
				printf("Reading errorin Cisbp/Homer matrix file %s !\n", file_pfm);
				return -1;
			}
		}
	}
	rewind(in_pfm);
	int olen = 0;
	double olen_perf = 0;
	fgets(head, sizeof(head), in_pfm);
	if (mtype == 0)
	{
		for (i = 0; i < PWMLEN; i++)
		{
			if (fgets(d1, sizeof(d1), in_pfm) != NULL)
			{
				DelChar(d1, '\n');
				DelChar(d1, '\r');
				char c = d1[0];
				if (isdigit(c) || (strchr("-ATGC", c) != 0))
				{
					olen++;
					//	 printf("%s",d[i]);						
					for (j = 0; j < alfabet; j++)
					{
						/*if (NthColumn(d, s, j + shift_col) == -1)
						{
							break;
						}*/
						test = UnderStol(d1, j + shift_col, val, sizeof(val), sep);
						if (test == -1) { printf("Wrong format %s\n", d1); exit(1); }
						double score = atof(val);
						pfm[i][j] = nseq * score;
					}
				}
				else break;
			}
			else break;
		}
	}
	else
	{
		fgets(d1, sizeof(d1), in_pfm);
		int dlen = strlen(d1);
		olen = 0;
		for (i = 1; i < dlen; i++)
		{
			if (d1[i - 1] == sep && isdigit(d1[i]))olen++;
		}
		rewind(in_pfm);
		fgets(head, sizeof(head), in_pfm);
		for (i = 0; i < alfabet; i++)
		{
			if (fgets(d1, sizeof(d1), in_pfm) != NULL)
			{
				DelChar(d1, '\n');
				char c = d1[0];
				if (isdigit(c) || (strchr("-ATGC", c) != 0))
				{
					//	 printf("%s",d[i]);	
					for (j = 0; j < olen; j++)
					{
						/*if (NthColumn(d, s, j + shift_col) == -1)
						{
							break;
						}*/
						test = UnderStol(d1, j + shift_col, val, sizeof(val), sep);
						if (test == -1) { printf("Wrong format %s\n", d1); exit(1); }
						double score = atof(val);
						pfm[j][i] = nseq * score;
					}
				}
				else break;
			}
			else break;
		}
	}
	fclose(in_pfm);
	len1 = lenp + word - 1;//dlina vyborki obu4eniya
	int ten[6] = { 1, 4, 16, 64, 256, 1024 };
	double min = 0, raz = 0;
	PWMScore(min, raz, lenp, olig);
	double min0 = min, raz0 = raz;	
	int nseq_pro = 0, len_pro = 0;
	int all_pos = 0;
	ReadSeq(file_seq, nseq_pro, len_pro,all_pos);
	double pvalue_large1 = pvalue_large*1.001;
	int nthr = (int)(pvalue_large1*all_pos);
	double *thr;
	thr = new double[nthr];
	if (thr == NULL) { puts("Out of memory..."); return -1; }
	for (i = 0; i < nthr; i++)thr[i] = 0;
	int nthr_max = nthr - 1;
	char **dp;
	dp = new char*[2];
	if (dp == NULL) { puts("Out of memory..."); return -1; }
	for (n = 0; n < 2; n++)
	{
		dp[n] = new char[len_pro + 1];
		if (dp[n] == NULL) { puts("Out of memory..."); return -1; }
	}
	if ((in = fopen(file_seq, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", file_seq);
		return -1;
	}	
	double all_pos_rec = 0;
	double score_min = 0;	
	int count_val= 0;
	double thr_bot = 0;
	for (i = 0; i < nthr; i++)thr[i] = thr_bot;
	for (n = 0; n<nseq_pro; n++)
	{		
		fgets(head, sizeof(head), in);
		memset(dp[0], 0, len_pro + 1);
		fgets(dp[0], len_pro + 1, in);
		DelChar(dp[0], '\n');
		TransStr(dp[0]);
		int len_pro1 = strlen(dp[0]);
		strcpy(dp[1], dp[0]);
		ComplStr(dp[1]);
		int len21 = len_pro1 - len1;
		if (n % 100 == 0)
		{
			int di = nthr_max / 10;
			printf("%d\t", n + 1);
			for (i = 0; i < nthr_max; i += di)printf("%d %f ", i + 1, thr[i]);
			printf("\n");
			//printf("%5d %f\n", n, thr[nthr_max]);
		}
		for (i = 0; i <= len21; i++)
		{						
			char d2[PWMLEN];
			int cod[PWMLEN];
			int gom = 0;
			double sco2 = -1000;
			for (k = 0; k < 2; k++)
			{
				int kpos;
				if (k == 0)kpos = i;
				else kpos = len21 - i;
				strncpy(d2, &dp[k][i], len1);
				d2[len1] = '\0';
				if (strstr(d2, "n") != NULL) { gom = -1; break; }
				gom = GetSostPro(d2, word, cod, letter);
				if (gom == -1)break;
				double score = 0;
				for (j = 0; j < lenp; j++)
				{
					score += pwm[j][cod[j]];
				}
				score -= min0;
				score /= raz0;
				if (score > sco2)sco2 = score;
			}
			if (gom == 0)
			{
				all_pos_rec++;
				double thr_check = Max(thr_bot, thr[nthr_max]);
				if (sco2 >= thr_check)
				{
					int gomt = 0;
					for (j = 0; j < nthr; j++)
					{
						if (sco2 >= thr[j])
						{
							int ksta = Min(nthr_max, count_val);
							for (k = ksta; k > j; k--)
							{
								//									if (thr[k] == 0)continue;
								Mix(&thr[k - 1], &thr[k]);
							}
							thr[j] = sco2;
							gomt = 1;
							break;
						}
						if (gomt == 1)break;
					}
					count_val++;
				}
				all_pos_rec++;
			}
		}
		//fprintf(out_distt, "%d\t%.1f\n", n+1,all_pos_rec);
	}
	fclose(in);
	if ((out_distt = fopen(file_out_distt, "wt")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_out_distt);
		return -1;
	}
	int nthr_dist = 0;
	int nthr_final = nthr-1;
	double fpr_pred = (double)1 / all_pos_rec;
	double fpr_pred_step = fpr_pred;
	//double thr_pred = thr[0];
	for (j = 1; j < nthr; j++)
	{
		double fpr = (double)(j+1) / all_pos_rec;
		double thr_pred = thr[j - 1];
		if ((thr[j] != thr_pred && fpr - fpr_pred_step > bin) || j == nthr_final)
		{
			nthr_dist++;
		//	fprintf(out_distt, "%.18f\t%.18g\n", thr_pred, fpr_pred);
			if (fpr_pred >= pvalue_large)
			{
				break;
			}
			fpr_pred_step = fpr;
		}
		fpr_pred = fpr;
	}
	double *thr_dist, *fpr_dist;
	thr_dist = new double[nthr_dist];
	if (thr_dist == NULL) { puts("Out of memory..."); return -1; }
	fpr_dist = new double[nthr_dist];
	if (fpr_dist == NULL) { puts("Out of memory..."); return -1; }
	int count = 0;
	fpr_pred = (double)1 / all_pos_rec;
	fpr_pred_step = fpr_pred;
	for (j = 1; j < nthr; j++)
	{
		double fpr = (double)(j + 1) / all_pos_rec;
		double thr_pred = thr[j - 1];
		if ((thr[j] != thr_pred && fpr - fpr_pred_step > bin) || j == nthr_final)
		{						
			if (j != nthr_final)
			{
				thr_dist[count] = thr_pred;
				fpr_dist[count] = -log10(fpr_pred);
			}
			else
			{
				thr_dist[count] = thr[j];
				fpr_dist[count] = -log10(fpr);					
			}
			fprintf(out_distt, "%.18f\t%.18g\n", thr_dist[count], fpr_dist[count]);			
			if (fpr_dist[count] <= log_pvalue_large)
			{
				count++;
				break;
			}				
			count++;
			fpr_pred_step = fpr;				
		}		
		fpr_pred = fpr;
	}
	fclose(out_distt);
	for (i = 0; i < 4; i++)for (j = 0; j < lenp; j++)pfm[j][i] /= nseq;
	fpr_pred = (double)1 / all_pos_rec;
	double fpr1st = 0;
	for (j = 1; j < nthr; j++)
	{
		if ((thr[j] != thr[j-1]) || j == nthr_final)
		{
			fpr1st = -log10((double)j / all_pos_rec);
			break;
		}
	}
	if ((check_bad == 1 && fpr1st > log_thr_best_pval) || check_bad == 0)
	{
		int lenp4 = 4 * lenp;
		/*if ((out_distb = fopen(file_out_distb, "ab")) == NULL)
		{
			printf("Out file %s can't be opened!\n", file_out_distb);
			return -1;
		}
		fwrite(&lenp, sizeof(int), 1, out_distb);		
		fwrite(pfm, sizeof(double), lenp4, out_distb);
		fwrite(pwm, sizeof(double), lenp4, out_distb);
		fwrite(&count, sizeof(int), 1, out_distb);
		fwrite(thr_dist, sizeof(double), count, out_distb);
		fwrite(fpr_dist, sizeof(double), count, out_distb);
		fclose(out_distb);*/
		//char  file_out_distb_one[300];
	//	strcpy(file_out_distb_one, file_pfm);
	//	strcat(file_out_distb_one, binary_ext);
		if ((out_distb = fopen(file_out_distb, binary_mode)) == NULL)
		{
			printf("Out file %s can't be opened!\n", file_out_distb);
			return -1;
		}
		fwrite(&lenp, sizeof(int), 1, out_distb);		
		fwrite(pfm, sizeof(double), lenp4, out_distb);
		fwrite(pwm, sizeof(double), lenp4, out_distb);
		fwrite(&count, sizeof(int), 1, out_distb);
		fwrite(thr_dist, sizeof(double), count, out_distb);
		fwrite(fpr_dist, sizeof(double), count, out_distb);
		fclose(out_distb);
	}
	FILE *out_sta;
	if ((out_sta = fopen(file_sta, "at")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_sta);
		return -1;
	}
	fprintf(out_sta, "%s\t%d\t%f\t%f\t", file_out_distt, nthr_dist,thr[0], fpr1st);
	if (fpr1st > log_thr_best_pval)fprintf(out_sta, "Good");
	else fprintf(out_sta, "Bad");
	fprintf(out_sta, "\n");
	fclose(out_sta);
	delete[] thr_dist;
	delete[] fpr_dist;
	delete[] thr;	
	for(n=0;n<2;n++)delete[] dp[n];
	delete[] dp;
	return 1;
}
