void PWMScore(double &min, double &raz, int len1, double(*pwm)[OLIGNUM])
{
	int i, j;
	for (i = 0; i < len1; i++)
	{
		double pwmmin = 100;
		double pwmmax = -100;
		for (j = 0; j < OLIGNUM; j++)
		{
			if (pwm[i][j] < pwmmin)pwmmin = pwm[i][j];
			if (pwm[i][j] > pwmmax)pwmmax = pwm[i][j];
		}
		raz += pwmmax;
		min += pwmmin;
	}
	raz -= min;
}
// ras4et 4asot oligonukleotidov po stroke (zdes' - nukleotidov)
int GetSostPro(char* d, int word, int* sost)
{
	int i, j, k, i_sost, let;
	int ten[6] = { 1, 4, 16, 64, 256, 1024 };
	char letter[5];
	strcpy(letter, "acgt");//atgc
	letter[4] = '\0';
	int lens = (int)strlen(d);
	int size = 1;
	for (k = 0; k < word; k++)size *= 4;
	for (i = 0; i < size; i++)sost[i] = 0;
	for (i = 0; i < lens - word + 1; i++)
	{
		i_sost = 0;
		let = -1;
		for (j = word - 1; j >= 0; j--)
		{
			for (k = 0; k < 4; k++)
			{
				if (d[i + j] == letter[k]) { let = k; break; }
			}
			if (let == -1)return -1;
			i_sost += ten[word - 1 - j] * let;
		}
		sost[i] = i_sost;
	}
	return 0;
}
void Mix(double *a, double *b)
{
	double buf = *a;
	*a = *b;
	*b = buf;
}
int pwm_iz_pwm_thr_dist0(double pwm_source[][OLIGNUM], int lenp, char *file_pro, int nthr, int &nthr_dist, double *thr, double *fpr, char *species, int nseq_pro, int len_pro, double pvalue, double dpvalue)
{
	int i, j, k, n;
	char head[1000], dp[2][SEQLEN];
	FILE *in;

	int nseq = 0;
	int len1 = 0;
	int word = 1;
	len1 = lenp + word - 1;//dlina vyborki obu4eniya
	int ten[6] = { 1, 4, 16, 64, 256, 1024 };
	double min = 0, raz = 0;
	PWMScore(min, raz, lenp, pwm_source);
	double min0 = min, raz0 = raz;
	int cod[SEQLEN];
	if ((in = fopen(file_pro, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", file_pro);
		return -1;
	}
	int all_pos = 0;
	int count_val = 0;
	int nthr_max = nthr - 1;
	double thr_bot = 0;
	for (i = 0; i < nthr; i++)thr[i] = thr_bot;	
	for (n = 0; n < nseq_pro; n++)
	{
		fgets(head, sizeof(head), in);
		fgets(dp[0], len_pro + 2, in);
		DelChar(dp[0], '\n');
		int len_pro1 = (int)strlen(dp[0]);
		int len21 = len_pro1 - len1;
		TransStr(dp[0]);
		strcpy(dp[1], dp[0]);
		ComplStr(dp[1]);				
	//	if (n % 50 == 0)printf("%5d %f\t%d\n", n, thr[nthr_max], count_val);
		/*if (n % 1000 == 0)
		{
			int di = nthr_max / 10;
			printf("%d\t", n + 1);
			for (i = 0; i < nthr_max; i += di)printf("%d %f ", i + 1, thr[i]);
			printf("\n");			
		}*/
		for (i = 0; i <= len21; i++)
		{			
			char d2[SEQLEN];
			int gom = 0;
			double sco2 = -1000;
			int compl1;
			for (compl1 = 0; compl1 < 2; compl1++)
			{
				int kpos;
				if (compl1 == 0)kpos = i;
				else kpos = len21 - i;
				strncpy(d2, &dp[compl1][i], len1);
				d2[len1] = '\0';
				if (strstr(d2, "n") != NULL) { gom = -1; break; }
				gom = GetSostPro(d2, word, cod);
				if (gom == -1)break;
				double score = 0;
				for (j = 0; j < lenp; j++)
				{
					score += pwm_source[j][cod[j]];
				}
				score -= min0;
				score /= raz0;
				if (score > sco2)sco2 = score;
			}
			if (gom == 0)
			{
				all_pos++;
				double thr_check = Max(thr_bot, thr[nthr_max]);
				if (sco2 >= thr_check)
				{
					int gomc = 0;
					for (j = 0; j < nthr; j++)
					{
						if (sco2 >= thr[j])
						{
							//if (thr[j] != 0)
							{
								int ksta = Min(nthr_max, count_val);
								for (k = ksta; k > j; k--)
								{
									Mix(&thr[k - 1], &thr[k]);
								}
							}
							thr[j] = sco2;
							gomc = 1;
							break;
						}
						if (gomc == 1)break;
					}
					count_val++;
				}
			}						
		}
	}
	fclose(in);
	nthr_dist = 0;
	int nthr_final = nthr - 1;
//	all_prom_pos = (int)all_pos;
	double fpr_pred_step = (double)1 / all_pos;	
	double fpr_pred = fpr_pred_step;
	for (j = 1; j < nthr; j++)
	{
		double fpr = (double)(j + 1) / all_pos;
		double thr_pred = thr[j-1];
		if ((thr[j] != thr_pred && fpr - fpr_pred_step > dpvalue) || j == nthr_final)
		{
			nthr_dist++;		
			if (fpr_pred >= pvalue)
			{
				break;
			}			
			fpr_pred_step = fpr;
		}
		fpr_pred = fpr;
	//	if (j == nthr_final && thr[j] == thr_pred)nthr_dist++;
	}
	double *thr_dist, *fpr_dist;
	thr_dist = new double[nthr_dist];
	if (thr_dist == NULL) { puts("Out of memory..."); return -1; }
	fpr_dist = new double[nthr_dist];
	if (fpr_dist == NULL) { puts("Out of memory..."); return -1; }
	int count = 0;
	fpr_pred = (double)1 / all_pos;
	fpr_pred_step = fpr_pred;	
	for (j = 1; j < nthr; j++)
	{
		double fpr = (double)(j + 1) / all_pos;
		double thr_pred = thr[j-1];
		if ((thr[j] != thr_pred && fpr - fpr_pred_step > dpvalue) || j == nthr_final)
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
			count++;
			if (fpr_pred >= pvalue)
			{
				break;
			}
			fpr_pred_step = fpr;										
		}
		fpr_pred = fpr;	   
	}
	nthr_dist = count;
	for (j = 0; j < count; j++)
	{
		thr[j] = thr_dist[j];
		fpr[j] = fpr_dist[j];
	}
	delete[] thr_dist;
	delete[] fpr_dist;
	return 1;
}
