int select_thresholds_from_pvalues(int n_thr_touzet, double *thr_touzet, double *fp_rate, double pvalue_large, double ratio, 
								   double *fpr_select, double *thr_select, int index[NUM_THR])
{
    int i, j;

	double ratio_cur=ratio;
	fpr_select[NUM_THR-1]=pvalue_large;
	for(i=NUM_THR-2;i>=0;i--)
	{
		fpr_select[i]=fpr_select[i+1]/ratio_cur;
		ratio_cur=1+ratio_cur/2;
	}
	if (n_thr_touzet <= NUM_THR) return -1;
	if (fp_rate[0] > fpr_select[0])
	{		
		for (i = 0; i < NUM_THR; i++)index[i] = i;
	}
	else
	{
		int jsta=n_thr_touzet-1;
		for(i=NUM_THR-1;i>=0;i--)
		{
			int cat=0;
			double fpr=fpr_select[i];
			for(j=jsta;j>0;j--)
			{
				int j1=j-1;
				if(fp_rate[j]>=fpr && fp_rate[j1]<fpr)
				{
					index[i]=j1;
					jsta=j;
					cat=1;
				}
			}
			if(cat==0)index[i]=0;
		}
	}	
	for(i=0;i<NUM_THR;i++)
	{
		fpr_select[i]=fp_rate[index[i]];
		printf("\t%g",fpr_select[i]);
	}
	printf("\n");
	for(i=0;i<NUM_THR;i++)
	{
		thr_select[i]=thr_touzet[index[i]];
		printf("\t%g",thr_select[i]);
	}
	printf("\n");	
	return 1;
}
