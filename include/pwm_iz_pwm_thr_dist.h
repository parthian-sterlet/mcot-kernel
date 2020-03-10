void PWMScore(double &min,double &raz, int len1, double (*pwm)[OLIGNUM])
{
	 int i,j;
	  for(i=0;i<len1;i++)
	  {
		  double pwmmin=100;
		  double pwmmax=-100;
		  for(j=0;j<OLIGNUM;j++)
		  {							 
			  if(pwm[i][j]<pwmmin)pwmmin=pwm[i][j];
			  if(pwm[i][j]>pwmmax)pwmmax=pwm[i][j];
		  }
		  raz+=pwmmax;
		  min+=pwmmin;
	 }	
	raz-=min;
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
	char dp[SEQLEN];
	double p;
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
	double all_pos = 0;	
	int count_val = 0;
	int nthr_max = nthr - 1;
	double score_min = 0.75;
	double pvalue2 = pvalue * 2;
	for (n = 0; n<nseq_pro; n++)
	{
		//if (n % 100 == 0)printf("%5d %f\t%d\n", n, thr[nthr_max],count_val);
		fgets(dp, len_pro + 2, in);
		DelChar(dp, '\n');
		int len_pro1 = strlen(dp);
		int len21 = len_pro1 - len1;
		double thresh_min;
		if(n==0)thresh_min = score_min;//double thresh_min = Max(thr[nthr_max], score_min);
		else
		{
			if (count_val >= nthr)thresh_min = thr[nthr_max];
			else
			{
				int n_check = (int)(all_pos*pvalue2)-1;
				n_check = (int)((n*n_check + count_val) / (n+1));
				thresh_min = thr[n_check];
			}
		}
		int compl1;
		for (compl1 = 0; compl1<2; compl1++)
		{
			if (compl1 == 1) if (ComplStr(dp) != 1)  { puts("Out of memory..."); return -1; }
			char d2[SEQLEN];
			p = -1000;
			for (i = 0; i <= len21; i++)
			{
				strncpy(d2, &dp[i], len1);
				d2[len1] = '\0';
				if (strstr(d2, "n") != NULL) { continue; }
				GetSostPro(d2, word, cod);				
				all_pos++;
				double score = 0;
				for (j = 0; j<lenp; j++)
				{
					score += pwm_source[j][cod[j]];
				}
				score -= min0;
				score /= raz0;				
				if (score >= thresh_min)
				{
					int gom = 0;
					for (j = 0; j < nthr; j++)
					{
						if (score >= thr[j])
						{
							//if (thr[j] != 0)
							{
								int ksta = Min(nthr_max, count_val);
								for (k = ksta; k > j; k--)
								{
									Mix(&thr[k - 1], &thr[k]);
								}
							}
							thr[j] = score;
							gom = 1;
							break;
						}
						if (gom == 1)break;
					}
					count_val++;
				}
			}
		}
	}
	fclose(in);
	double fpr_pred=1/(double)all_pos;
	double thr0 = thr[0];
	nthr_dist = 0;
	int nthr_final = nthr;
	for (j = 1; j < nthr; j++)
	{
		double fprx = (double)(j + 1) / all_pos;
		if ((thr[j] != thr0 || j == nthr_max) && fprx - fpr_pred > dpvalue)
		{
			nthr_dist++;
			thr0 = thr[j];			
			if (fprx >= pvalue)
			{
				nthr_final = j;
				break;
			}
			fpr_pred = fprx;
		}
	}
	double *thr_dist, *fpr_dist;
	thr_dist = new double[nthr_dist];
	if (thr_dist == NULL) { puts("Out of memory..."); return -1; }
	fpr_dist = new double[nthr_dist];
	if (fpr_dist == NULL) { puts("Out of memory..."); return -1; }
	int count = 0;
	thr0 = thr[0];
	nthr_max = nthr_final;
	fpr_pred = 1 / (double)all_pos;
	for (j = 1; j <= nthr_final; j++)
	{
		double fprx = (double)(j + 1) / all_pos;
		if ((thr[j] != thr0 || j == nthr_max) && fprx - fpr_pred > dpvalue)
		{			
			thr_dist[count] = thr0;
			fpr_dist[count] = fprx;
			thr0 = thr[j];			
			fpr_pred = fprx;
			count++;
		}
	}
	for (j = 0; j < nthr_dist; j++)
	{
		thr[j] = thr_dist[j];
		fpr[j] = fpr_dist[j];
	}
	delete[] thr_dist;
	delete[] fpr_dist;
	return 1;
}
