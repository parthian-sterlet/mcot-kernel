int select_thresholds_from_pvalues(int n_thr_touzet, double *thr_touzet, double *fp_rate, double *fpr_select_i,
								   double *fpr_select_o, double *thr_select, int index[NUM_THR])
{
    int i, j;

	if (n_thr_touzet <= NUM_THR) return -1;
	if (fp_rate[0] < fpr_select_i[0])
	{		
		for (i = 0; i < NUM_THR; i++)index[i] = i;
	}
	else
	{
		int jsta=n_thr_touzet-1;
		for(i=NUM_THR-1;i>=0;i--)
		{
			int cat=0;
			double fpr=fpr_select_i[i];
			for(j=jsta;j>0;j--)
			{
				int j1=j-1;
				if(fp_rate[j1]> fpr && fpr>=fp_rate[j])
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
		fpr_select_o[i]=fp_rate[index[i]];
		printf("\t%g",fpr_select_o[i]);
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
