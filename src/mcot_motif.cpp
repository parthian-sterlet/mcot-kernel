#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define PWMLEN 100 // 2 * (max length of motives)

int ma_ca[PWMLEN][4];
int ma_cp[PWMLEN][4];
int mp_ca[PWMLEN][4];
int mp_cp[PWMLEN][4];
int m_ca[PWMLEN][4];//both
int m_cp[PWMLEN][4];//both

// razbor of asymmetric CEs with overlap
int StrNStr(char *str, char c, int n)
{
	if (n == 0)return -1;
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
void DelHole(char *str)
{
	char *hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
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
// delete symbol 'c' from input string
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
int IdeLet(char c)
{
	int ret;
	switch (c) {
	case 'a': ret = 0; break;
	case 'c': ret = 1; break;
	case 'g': ret = 2; break;
	case 't': ret = 3; break;
	default: ret = -1;
	}
	return(ret);
}
int IdeLetCom(char c)
{
	int ret;
	switch (c) {
	case 'a': ret = 3; break;
	case 'c': ret = 2; break;
	case 'g': ret = 1; break;
	case 't': ret = 0; break;
	default: ret = -1;
	}
	return(ret);
}
int main(int argc, char *argv[])
{
	int i, j, k;
	char file_hist[120] = { 0 }, file_ce[120] = { 0 }, fileo[120] = { 0 }, fileo_mata[120] = { 0 }, fileo_matp[120] = { 0 };
	char name[10] = { 0 }, namea[10] = { 0 }, namep[10] = { 0 }, nameap[10] = { 0 }, partner_num[5] = { 0 }, partner_name[50] = { 0 };
	char d[4][5000] = { 0 }, s[500] = { 0 }, ori[3][15] = { 0 }, loc[5] = { 0 }, val[20] = { 0 }, loc_best[5] = { 0 }, ori_best[15] = { 0 };
	//	double *val;
	FILE *in, *out, *outa;// , *outp;

	if (argc != 5)
	{
		printf("%s 1file_input_hist 2file_input_CEs 3file out_sta 4file out_mat_anc", argv[0]);// f f 9 2 10 5 5
		return -1;
	}
	strcpy(file_hist, argv[1]);//in_file
	strcpy(file_ce, argv[2]);//in_file
	strcpy(fileo, argv[3]);//out_file_stat
	strcpy(fileo_mata, argv[4]);//out_file_matrix anchor
	//strcpy(fileo_matp, argv[5]);//out_file_matrix partner
	double fr, fr_max = -1;
	double thr = 0.9;
	if ((in = fopen(file_hist, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", file_hist);
		return -1;
	}
	for (i = 0; i < 4; i++)
	{
		fgets(d[i], sizeof(d[i]), in);
		DelChar(d[i], '\n');
	}
	int test;
	test = UnderStol(d[0], 0, partner_num, sizeof(partner_num), '\t');
	if (test == -1) { printf("Wrong format %s\n", d[0]); exit(1); }
	test = UnderStol(d[0], 1, partner_name, sizeof(partner_name), '\t');
	if (test == -1) { printf("Wrong format %s\n", d[0]); exit(1); }
	for (i = 1; i < 4; i++)
	{
		int i1 = i - 1;
		memset(ori[i1], '\0', sizeof(ori[i1]));
		test = UnderStol(d[i], 2, ori[i1], sizeof(ori[i1]), '\t');
		if (test == -1) { printf("Wrong format %s\n", d[i]); exit(1); }
	}
	for (j = 3;; j++)
	{
		memset(loc, '\0', sizeof(loc));
		test = UnderStol(d[0], j, loc, sizeof(loc), '\t');
		if (test == -1) { printf("Wrong format %s\n", d[0]); exit(1); }
		if (strstr(loc, "S") != NULL)break;
		for (i = 1; i < 4; i++)
		{
			test = UnderStol(d[i], j, val, sizeof(val), '\t');
			if (test == -1) { printf("Wrong format %s\n", d[i]); exit(1); }
			fr = atof(val);
			if (fr > fr_max)
			{
				fr_max = fr;
				strcpy(ori_best, ori[i - 1]);
				strcpy(loc_best, loc);
			}
		}
	}
	fclose(in);
	strcpy(namea, "no");
	strcpy(namep, "no");
	strcpy(nameap, "no");
	int npeak_p = 0, npeak_a = 0, nce_a = 0, nce_p = 0, npeak = 0;
	if ((in = fopen(file_ce, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", file_ce);
		return -1;
	}
	int lena = 0, lenp = 0, asta = 0, aend = 0, psta = 0, pend = 0;
	double ca, cp;
	char sa[PWMLEN], sp[PWMLEN];
	int overlap, len_tot;
	fgets(s, sizeof(s), in);
	while (fgets(s, sizeof(s), in) != NULL)
	{
		test = UnderStol(s, 5, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		if (strstr(val, loc_best) == NULL)continue;
		test = UnderStol(s, 8, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		if (strncmp(val, ori_best,6) != 0)continue;
		test = UnderStol(s, 11, sa, sizeof(sa), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		lena = strlen(sa);
		test = UnderStol(s, 12, sp, sizeof(sp), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		DelHole(sp);
		lenp = strlen(sp);
		test = UnderStol(s, 1, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		asta = atoi(val);
		test = UnderStol(s, 2, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		aend = atoi(val);
		test = UnderStol(s, 3, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		psta = atoi(val);
		test = UnderStol(s, 4, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		pend = atoi(val);
		int min = Min(asta, psta);
		min--;
		asta -= min;
		aend -= min;
		psta -= min;
		pend -= min;
		int minsta = Min(asta, psta);
		int minend = Min(aend, pend);
		int maxsta = Max(asta, psta);
		int maxend = Max(aend, pend);
		len_tot = maxend - minsta + 1;
		overlap = minend - maxsta + 1;
		test = UnderStol(s, 7, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		if (val[0] == '-')
		{
			int asta1 = len_tot - aend + 1;
			int aend1 = len_tot - asta + 1;
			int psta1 = len_tot - pend + 1;
			int pend1 = len_tot - psta + 1;
			asta = asta1;
			aend = aend1;
			psta = psta1;
			pend = pend1;
		}
		break;
	}
	int rev;
	if (strcmp(ori_best, "Everted") == 0)rev = 1;
	else rev = 0;
	{
	}
	for (i = 0; i < lena; i++)for (j = 0; j < 4; j++)ma_ca[i][j] = ma_cp[i][j] = 0;
	for (i = 0; i < lenp; i++)for (j = 0; j < 4; j++)mp_ca[i][j] = mp_cp[i][j] = 0;
	for (i = 0; i < len_tot; i++)for (j = 0; j < 4; j++)m_ca[i][j] = m_ca[i][j] = 0;
	rewind(in);
	fgets(s, sizeof(s), in);
	while (fgets(s, sizeof(s), in) != NULL)
	{
		test = UnderStol(s, 0, name, sizeof(name), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		test = UnderStol(s, 5, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		if (strstr(val, loc_best) == NULL)continue;
		test = UnderStol(s, 8, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		if (strncmp(val, ori_best, 6) != 0)continue;
		test = UnderStol(s, 9, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		ca = atof(val);
		test = UnderStol(s, 10, val, sizeof(val), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		cp = atof(val);
		test = UnderStol(s, 11, sa, sizeof(sa), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		if (strstr(sa, "n") != NULL)continue;
		test = UnderStol(s, 12, sp, sizeof(sp), '\t');
		if (test == -1) { printf("Wrong format %s\n", s); exit(1); }
		if (strstr(sp, "n") != NULL)continue;
		DelHole(sp);
		int r;
		if (strcmp(name, nameap) != 0)
		{
			npeak++;
			strcpy(nameap, name);
		}
		//if (ca > cp)
		{
			nce_a++;
			if (strcmp(name, namea) != 0)
			{
				npeak_a++;
				strcpy(namea, name);
			}
			for (j = 0; j < lena; j++)
			{
				r = IdeLet(sa[j]);
				ma_ca[j][r]++;
			}
			if (ori_best[0] == 'D')//DirectPA, DirectAP
			{
				for (j = 0; j < lenp; j++)
				{
					r = IdeLet(sp[j]);
					mp_ca[j][r]++;
				}
			}
			else
			{
				k = lenp - 1;
				for (j = 0; j < lenp; j++)
				{
					r = IdeLetCom(sp[j]);
					mp_ca[k--][r]++;
				}

			}
		}
		/*
		else
		{
			nce_p++;
			if (strcmp(name, namep) != 0)
			{
				npeak_p++;
				strcpy(namep, name);
			}
			for (j = 0; j < lena; j++)
			{
				r = IdeLet(sa[j]);
				ma_cp[j][r]++;
			}
			if (ori_best[0] == 'D')//DirectPA, DirectAP
			{
				for (j = 0; j < lenp; j++)
				{
					r = IdeLet(sp[j]);
					mp_cp[j][r]++;
				}
			}
			else
			{
				k = lenp - 1;
				for (j = 0; j < lenp; j++)
				{
					r = IdeLetCom(sp[j]);
					mp_cp[k--][r]++;
				}
			}
		}*/
	}
	fclose(in);
	out = outa = NULL;//outp = 
	if (strcmp(partner_num, "1") == 0)
	{
		if ((out = fopen(fileo, "wt")) == NULL)
		{
			printf("Output file %s can't be opened!", fileo);
			return -1;
		}
		if ((outa = fopen(fileo_mata, "wt")) == NULL)
		{
			printf("Output file %s can't be opened!", fileo_mata);
			return -1;
		}
		/*if ((outp = fopen(fileo_matp, "wt")) == NULL)
		{
			printf("Output file %s can't be opened!", fileo_matp);
			return -1;
		}*/
		//fprintf(out, "Partner_Number\tParner_Name\tOri\tLoc\tFr_max\tAnchorSta\tAnchorEnd\tPartnerSta\tPartnerEnd\tLenTot\tOverlap\tPeaks\tCEs\tAnchorPeaks\tCEs\tSeqCE\tSeqAnchor\tPartnerPeaks\tCEs\tSeqCE\tSeqAnchor\n");
	}
	else
	{
		if ((out = fopen(fileo, "at")) == NULL)
		{
			printf("Output file %s can't be opened!", fileo);
			return -1;
		}
		if ((outa = fopen(fileo_mata, "at")) == NULL)
		{
			printf("Output file %s can't be opened!", fileo_mata);
			return -1;
		}
		/*if ((outp = fopen(fileo_matp, "at")) == NULL)
		{
			printf("Output file %s can't be opened!", fileo_matp);
			return -1;
		}*/
	}
	fprintf(out, "Partner_Number\tParner_Name\tOri\tLoc\tFr_max\tAnchorSta\tAnchorEnd\tPartnerSta\tPartnerEnd\tLenTot\tOverlap\tPeaks\tCEs\tAnchorPeaks\tCEs\tSeqCE\tSeqAnchor\tPartnerPeaks\tCEs\tSeqCE\tSeqAnchor\n");
	int len_tot1 = len_tot - 1;
	char alfavit1[5] = "ACGT", alfavit1_small[5] = "acgt";
	char alfavit2[7] = "WRMKYS", alfavit2_small[7] = "wrmkys";
	char alfavit3[5] = "BDHV", alfavit3_small[5] = "bdhv";
	int pair[6][2] = { {0,3},{0,2},{0,1},{3,2},{3,1},{2,1} };
	char consena[PWMLEN], consenp[PWMLEN], consena_a[PWMLEN], consenp_a[PWMLEN];
	for (j = 0; j < len_tot; j++)consenp[j] = consena[j] = 'n';
	consena[len_tot] = consenp[len_tot] = '\0';
	//char letter[] = "ACGT";
	for (j = 0; j < 4; j++)
	{
		k = 0;
		if (rev == 0)
		{
			for (i = 0; i < lena; i++)m_ca[k++][j] = ma_ca[i][j];
			for (i = overlap; i < lenp; i++)m_ca[k++][j] = mp_ca[i][j];
		}
		else
		{
			for (i = 0; i < lenp; i++)m_ca[k++][j] = mp_ca[i][j];
			for (i = overlap; i < lena; i++)m_ca[k++][j] = ma_ca[i][j];
		}
	}
	for (j = 0; j < 4; j++)
	{
		k = 0;
		if (rev == 0)
		{
			for (i = 0; i < lena; i++)m_cp[k++][j] = ma_cp[i][j];
			for (i = overlap; i < lenp; i++)m_cp[k++][j] = mp_cp[i][j];
		}
		else
		{
			for (i = 0; i < lenp; i++)m_cp[k++][j] = mp_cp[i][j];
			for (i = overlap; i < lena; i++)m_cp[k++][j] = ma_cp[i][j];
		}
	}
	int nseq, fthr;
	nseq = nce_a;
	fthr = (int)(nseq*thr);
	if (nseq > 0)
	{
		for (i = 0; i < len_tot; i++)
		{
			for (j = 0; j < 4; j++)
			{
				int fr = m_ca[i][j];
				if (fr >= fthr)
				{
					if (fr == nseq)consena[i] = alfavit1[j];
					else consena[i] = alfavit1_small[j];
					break;
				}
			}
			if (consena[i] == 'n')
			{
				for (k = 0; k < 6; k++)
				{
					int fr = m_ca[i][pair[k][0]] + m_ca[i][pair[k][1]];
					if (fr >= fthr)
					{
						if (fr == nseq)consena[i] = alfavit2[k];
						else consena[i] = alfavit2_small[k];
						break;
					}
				}
			}
			if (consena[i] == 'n')
			{
				for (j = 0; j < 4; j++)
				{
					int fr = nseq - m_ca[i][j];
					if (fr >= fthr)
					{
						if (fr == nseq)consena[i] = alfavit3[j];
						else consena[i] = alfavit3_small[j];
						break;
					}
				}
			}
		}
	}
	nseq = nce_p;
	if (nseq > 0)
	{
		fthr = (int)(nseq*thr);
		for (i = 0; i < len_tot; i++)
		{
			for (j = 0; j < 4; j++)
			{
				int fr = m_cp[i][j];
				if (fr >= fthr)
				{
					if (fr == nseq)consenp[i] = alfavit1[j];
					else consenp[i] = alfavit1_small[j];
					break;
				}
			}
			if (consenp[i] == 'n')
			{
				for (k = 0; k < 6; k++)
				{
					int fr = m_cp[i][pair[k][0]] + m_ca[i][pair[k][1]];
					if (fr >= fthr)
					{
						if (fr == nseq)consenp[i] = alfavit2[k];
						else consenp[i] = alfavit2_small[k];
						break;
					}
				}
			}
			if (consenp[i] == 'n')
			{
				for (j = 0; j < 4; j++)
				{
					int fr = nseq - m_cp[i][j];
					if (fr >= fthr)
					{
						if (fr == nseq)consenp[i] = alfavit3[j];
						else consenp[i] = alfavit3_small[j];
						break;
					}
				}
			}
		}
	}
	{
		int ast = asta - 1;
		for (i = 0; i < lena; i++)
		{
			consena_a[i] = consena[i + ast];
			consenp_a[i] = consenp[i + ast];
		}
		consena_a[lena] = '\0';
		consenp_a[lena] = '\0';
	}
	fprintf(out, "%s\t%s\t%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t", partner_num, partner_name, ori_best, loc_best, fr_max, asta, aend, psta, pend, len_tot, overlap);
	fprintf(out, "%d\t%d\t", npeak, nce_a + nce_p);
	fprintf(out, "%d\t%d\t%s\t%s\t", npeak_a, nce_a, consena, consena_a);
	fprintf(out, "%d\t%d\t%s\t%s\n", npeak_p, nce_p, consenp, consenp_a);
	fprintf(outa, ">%s\t%s\t%s\t%s\tLen %d\tOverlap %d\tNpeak\t%d\t%s\tAnchor [%d,%d]\tPartner [%d,%d]\n", partner_num, partner_name, ori_best, loc_best, len_tot, overlap, npeak_a, consena, asta, aend, psta, pend);
	//fprintf(outp, ">%s\t%s\t%s\t%s\tLen %d\tOverlap %d\tNpeak\t%d\t%s\tAnchor [%d,%d]\tPartner [%d,%d]\n", partner_num, partner_name, ori_best, loc_best, len_tot, overlap, npeak_p, consenp, asta, aend, psta, pend);
/*	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < len_tot; i++)
		{
			fprintf(outp, "%d", m_cp[i][j]);
			if (i == len_tot1)fprintf(outp, "\n");
			else fprintf(outp, "\t");
		}
	}*/
	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < len_tot; i++)
		{
			fprintf(outa, "%d", m_ca[i][j]);
			if (i == len_tot1)fprintf(outa, "\n");
			else fprintf(outa, "\t");
		}
	}
	fclose(out);
	fclose(outa);
//	fclose(outp);
	return 1;
}
