void measure(int pfm1[][OLIGNUM], int pfm2[][OLIGNUM],  int len1, int len2,int alphabet, int size, double over_scores[][MATLEN], int measure_type)
{	
	int i,j,k;
	switch(measure_type) 
	{
		case 0://SSD
		{
			for(i=0;i<len1;i++)
			{
				for(j=0;j<len2;j++)
				{
					double ret=0;
					for(k=0;k<alphabet;k++)
					{
						double dif=(double)(pfm1[i][k]-pfm2[j][k])/size;
						ret+=dif*dif;
					}
					ret = 2-sqrt(ret);
					over_scores[i][j]=ret;
				}
			}
			break;
		}
		default://PCC
		{
			double av=1/(double)alphabet;
			double x2[MATLEN], y2[MATLEN];
			for(i=0;i<len1;i++)
			{
				x2[i]=0;
				for(k=0;k<alphabet;k++)
				{
					x2[i]+=((double)pfm1[i][k]/size-av)*((double)pfm1[i][k]/size-av);				
				}
				x2[i]=sqrt(x2[i]);
			}
			for(i=0;i<len2;i++)
			{
				y2[i]=0;
				for(k=0;k<alphabet;k++)
				{
					y2[i]+=((double)pfm2[i][k]/size-av)*((double)pfm2[i][k]/size-av);				
				}
				y2[i]=sqrt(y2[i]);
			}
			for(i=0;i<len1;i++)
			{
				for(j=0;j<len2;j++)
				{
					double xy=0;
					for(k=0;k<alphabet;k++)
					{
						xy+=((double)pfm1[i][k]/size-av)*((double)pfm2[j][k]/size-av);				
					}
					over_scores[i][j]=xy/(x2[i]*y2[j]);
				}
			}	
			break;
		}
	}
}
double Min_trackReal(double mat[][MATLEN],int size_short,int size_long, int track_len, int sta_sh, int sta_lo)//size1 <=size2
{
	double ret=-1000;
	int i_short, i_long, j;
	for(i_long=0;i_long<=size_long-track_len;i_long++)
	{
		for(i_short=0;i_short<=size_short-track_len;i_short++)
		{
			double sum=0;
			for(j=0;j<track_len;j++)
			{
				sum+=mat[i_short+j][i_long+j];				
			}			
			if(sum>ret)
			{
				ret=sum;
				sta_sh=i_short;
				sta_lo=i_long;
			}
		}
	}
	return ret;
}
double Min_trackRand(double mat[][MATLEN],int size_short,int size_long, int track_len)//size1 <=size2
{
	double ret=-1000;
	int i_short, i_long, j;
	for(i_long=0;i_long<=size_long-track_len;i_long++)
	{
		for(i_short=0;i_short<=size_short-track_len;i_short++)
		{
			double sum=0;
			for(j=0;j<track_len;j++)
			{
				double val=mat[i_short+j][i_long+j];
				sum+=val;				
			}
			//sum=sqrt(sum);
			if(sum>ret)ret=sum;
		}
	}
	return ret;
}
void CellCount(double mat[][MATLEN],int size_short,int size_long, double av[][MATLEN], double sd[][MATLEN])//size1 <=size2
{	
	int i_short, i_long;
	for(i_long=0;i_long<size_long;i_long++)
	{
		for(i_short=0;i_short<size_short;i_short++)
		{			
			double val=mat[i_short][i_long];
			av[i_short][i_long]+=val;
			sd[i_short][i_long]+=val*val;
		}
	}	
}
void Permute_columns(int pfm[][OLIGNUM], int len, int alphabet, int size)
{
	int i1, i2, j1, j2;
	int len1=len-1;
	int alphabet1=alphabet-1;	
	for(i1=0;i1<len;i1++)
	{
		i2=rand()%len1;
		if(i2>=i1)i2++;
		j1=rand()%alphabet;
		j2=rand()%alphabet1;
		if(j2==j1)j2++;
		j1c = rand() % alphabet;
		j2c = rand() % alphabet1;
		if (j2c == j1c)j2c++;
		int go_up=1, go_do=1;
		if(pfm[i1][j1]==size && pfm[i2][j1c] == size)go_up=0;
		if(pfm[i1][j2] == size && pfm[i2][j2c]==size)go_up=0;//
		if(pfm[i1][j1] == 0 && pfm[i2][j1c]==0)go_do=0;//
		if(pfm[i1][j2]==0 && pfm[i2][j2c] == 0)go_do=0;//
		if(go_up==0 && go_do==0)continue;		
		int up1,up2,up12=0,do1,do2,do12=0;
		if(go_up==1)
		{
			up1=Min(size-pfm[i1][j1],size-pfm[i2][j2c]);
			up2=Min(pfm[i2][j1c],pfm[i1][j2]);
			up12=Min(up1,up2);		
		}
		if(go_do==1)
		{
			do1=Min(pfm[i1][j1],pfm[i2][j2c]);
			do2=Min(size-pfm[i2][j1c],size-pfm[i1][j2]);
			do12=Min(do1,do2);
		}
		int sum_up_do=up12+do12;
		if(sum_up_do==0)
		{
			continue;
		}
		int delta;
		int dir=rand()%sum_up_do;
		if(dir<up12)delta=1+rand()%up12;
		else delta=-1-rand()%do12;
		//printf("Do- First col %d,%d Second col %d,%d\n",pfm[i1][j1],pfm[i1][j2],pfm[i2][j1],pfm[i2][j2]);
		pfm[i1][j1]+=delta;
		pfm[i1][j2]-=delta;
		pfm[i2][j1c]-=delta;
		pfm[i2][j2c]+=delta;
		/*if((pfm[i1][j1]>size || pfm[i1][j2]<0) || (pfm[i2][j1c]<0 || pfm[i2][j2c]>size))
		{
			int yy=0;
		}*/
		//printf("Po- First col %d,%d Second col %d,%d\n",pfm[i1][j1],pfm[i1][j2],pfm[i2][j1],pfm[i2][j2]);
	}
}
double Laplace(double z)
{
    const double pi=3.141592653589793; 
    const double step=0.0001;
	const double step2=step/2;
    double x, y=0, dy;    
    double prec = 1E-20;	
	x=z;	    
	double mnoj=sqrt(2*pi);
	prec/=mnoj;
	do
    {
		x+=step2;
		dy=step*exp(-x*x/2);
		//dy/=mnoj;
		y+=dy;
		x+=step2;
    }		
	while(dy>prec);
	y/=mnoj;
	y*=2;
    return(y);
}
void PFM_compl(int pfm[][OLIGNUM], int pfmc[][OLIGNUM], int olen)
{
	int i,j;
	int rel_mono[4]={3,2,1,0};
	for(i=0;i<olen;i++)
	{
		int ii=olen-1-i;		
		for(j=0;j<OLIGNUM;j++)
		{
			pfmc[i][j]=pfm[ii][rel_mono[j]];
		}
	}	
}
double pfm_similarity(matrices *mat1, matrices *mat2, double granul, int overlap_min, int n_cycle_small, int n_cycle_large, double pval[2], int reini_permut)//int measure_type
{	
	//int ***pfm, ***pfmr, ***pfmc;	
	int pfm[2][MATLEN][OLIGNUM];
	int pfmc[2][MATLEN][OLIGNUM];
	int pfmr[2][MATLEN][OLIGNUM];
	int i, j, k, n; 
	//int overlap_min=6;

	//measure_type=0   euclid
	//measure_type=1   pearson

	//int n_cycle_large=20000;
	//int n_cycle_small=1000;
	//double granul=0.001;
	int granul_back=(int)(1/granul);

	//start
	//srand((unsigned)time(NULL));
	int olen[2];
	olen[0]=mat1->len;
	olen[1]=mat2->len;
	/*pfm= new int**[2];
	if(pfm==NULL){puts("Out of memory...");return -1;}
	pfmc= new int**[2];
	if(pfmc==NULL){puts("Out of memory...");return -1;}
	pfmr= new int**[2];
	if(pfmr==NULL){puts("Out of memory...");return -1;}
	for(k=0;k<2;k++)
	{
		pfm[k] = new int*[olen[k]];				
		if(pfm[k]==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<olen[k];i++)
		{			
			pfm[k][i] = new int[OLIGNUM];
			if(pfm[k][i]==NULL){puts("Out of memory...");return -1;}
		}
	}*/				
	//okruglenue i normirovka
	for(k=0;k<2;k++)
	{
		for(i=0;i<olen[k];i++)
		{						
			double sum=0;		
			for(j=0;j<OLIGNUM;j++)
			{			
				double val;
				if(k==0)val=mat1->fre[i][j];
				else val=mat2->fre[i][j];
				sum+=val;
			}
			for(j=0;j<OLIGNUM;j++)
			{				
				double val;
				if(k==0)val=mat1->fre[i][j];
				else val=mat2->fre[i][j];	
				pfm[k][i][j]=(int)(granul_back*val/sum);
			}
		}
	}
	int gom=0;
	for(k=0;k<2;k++)
	{
		for(i=0;i<olen[k];i++)
		{						
			int sum=0;		
			for(j=0;j<OLIGNUM;j++)
			{				
				sum+=pfm[k][i][j];
			}
			//for(j=0;j<OLIGNUM;j++)printf("%d\t",pfm[k][i][j]);
			//printf("%d\n",sum);
			if(sum!=granul_back)gom=1;
		}
		//printf("\n");
	}
	//popravka okrugleniya
	if(gom==1)
	{		
		for(k=0;k<2;k++)
		{
			for(i=0;i<olen[k];i++)
			{						
				int sum=0;
				for(j=0;j<OLIGNUM;j++)sum+=pfm[k][i][j];				
				int delta=granul_back-sum;
				if(delta!=0)
				{			
					for(j=0;j<abs(delta);j++)
					{
						int rr=rand()%sum;
						int v, sum2=0;
						for(v=0;v<OLIGNUM;v++)
						{
							sum2+=pfm[k][i][v];
							if(rr<=sum2)
							{
								break;
							}
						}
						if(sum>0){pfm[k][i][v]++;sum++;}
						else {pfm[k][i][v]--;sum--;}
					}
				}				
			}
		}
	}	
	/*
	for(k=0;k<2;k++)
	{
		printf("Matrix %d\n",k+1);
		for(i=0;i<olen[k];i++)
		{						
			int sum=0;		
			for(j=0;j<OLIGNUM;j++)
			{				
				sum+=pfm[k][i][j];
			}
			for(j=0;j<OLIGNUM;j++)printf("%d\t",pfm[k][i][j]);
			printf("%d\n",sum);
		}
		printf("\n");
	}*/
	/*
	for(k=0;k<2;k++)
	{
		pfmr[k] = new int*[olen[k]];				
		if(pfmr[k]==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<olen[k];i++)
		{			
			pfmr[k][i] = new int[OLIGNUM];
			if(pfmr[k][i]==NULL){puts("Out of memory...");return -1;}
		}
	}	
	for(k=0;k<2;k++)
	{
		pfmc[k] = new int*[olen[k]];				
		if(pfmc[k]==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<olen[k];i++)
		{			
			pfmc[k][i] = new int[OLIGNUM];
			if(pfmc[k][i]==NULL){puts("Out of memory...");return -1;}
		}
	}*/	
	for(k=0;k<2;k++)
	{
		for(i=0;i<olen[k];i++)
		{						
			for(j=0;j<OLIGNUM;j++)
			{				
				pfmc[k][i][j]=pfm[k][i][j];
			}
		}
		PFM_compl(pfm[k],pfmc[k],olen[k]);
	}
	int kk[2];
	if(olen[0]<=olen[1])
	{
		kk[0]=0;kk[1]=1;
	}
	else 
	{
		kk[0]=1,kk[1]=0;
	}
	int k0=kk[0], k1=kk[1];
	int olen_min=olen[k0], olen_max=olen[k1];
	int overlap_max=olen_min;
	int n_track=olen_min-overlap_min+1;
	int ori;
	double over_scores_real[2][MATLEN][MATLEN];
	double over_scores_rand_one[2][MATLEN][MATLEN];
	double over_scores_rand_two[2][MATLEN][MATLEN];
/*
	double ***over_scores_real;		//indices 1st short 2nd long
	over_scores_real=new double**[2];
	if(over_scores_real==NULL){puts("Out of memory...");return -1;} 
	for(k=0;k<2;k++)// strands
	{
		over_scores_real[k] = new double*[olen_min];				
		if(over_scores_real[k]==NULL){puts("Out of memory...");return -1;} 
		for(i=0;i<olen_min;i++)
		{			
			over_scores_real[k][i] = new double[olen_max];
			if(over_scores_real[k][i]==NULL){puts("Out of memory...");return -1;}
		}		
	}	
	double ***over_scores_rand_one;		//indices 1st short 2nd long	
	over_scores_rand_one=new double**[2];
	if(over_scores_rand_one==NULL){puts("Out of memory...");return -1;} 
	for(k=0;k<2;k++)// strands
	{
		over_scores_rand_one[k] = new double*[olen_min];				
		if(over_scores_rand_one[k]==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<olen_min;i++)
		{			
			over_scores_rand_one[k][i] = new double[olen_max];
			if(over_scores_rand_one[k][i]==NULL){puts("Out of memory...");return -1;}
		}		
	}
	double ***over_scores_rand_two;		//indices 1st short 2nd long
	over_scores_rand_two=new double**[2];
	if(over_scores_rand_two==NULL){puts("Out of memory...");return -1;} 
	for(k=0;k<2;k++)// strands
	{
		over_scores_rand_two[k] = new double*[olen_min];				
		if(over_scores_rand_two[k]==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<olen_min;i++)
		{			
			over_scores_rand_two[k][i] = new double[olen_max];
			if(over_scores_rand_two[k][i]==NULL){puts("Out of memory...");return -1;}
		}		
	}*/
	double av_cell[2][MATLEN][MATLEN];
	double sd_cell[2][MATLEN][MATLEN];
	/*
	double ***av_cell;		//indices 1st short 2nd long
	av_cell=new double**[2];
	if(av_cell==NULL){puts("Out of memory...");return -1;} 
	for(k=0;k<2;k++)// strands
	{
		av_cell[k] = new double*[olen_min];				
		if(av_cell[k]==NULL){puts("Out of memory...");return -1;} 
		for(i=0;i<olen_min;i++)
		{			
			av_cell[k][i] = new double[olen_max];
			if(av_cell[k][i]==NULL){puts("Out of memory...");return -1;}
		}		
	}
	double ***sd_cell;		//indices 1st short 2nd long
	sd_cell=new double**[2];
	if(sd_cell==NULL){puts("Out of memory...");return -1;} 
	for(k=0;k<2;k++)// strands
	{
		sd_cell[k] = new double*[olen_min];				
		if(sd_cell[k]==NULL){puts("Out of memory...");return -1;} 
		for(i=0;i<olen_min;i++)
		{			
			sd_cell[k][i] = new double[olen_max];
			if(sd_cell[k][i]==NULL){puts("Out of memory...");return -1;}
		}		
	}*/
	for(k=0;k<2;k++)for(i=0;i<olen_min;i++)for(j=0;j<olen_max;j++)av_cell[k][i][j]=sd_cell[k][i][j]=0;	
	
	double sim_real[MATLEN], sim_rand[MATLEN];
	/*
	double *sim_real;
	sim_real=new double[n_track];
	if(sim_real==NULL){puts("Out of memory...");return -1;}	
	double *sim_rand;
	sim_rand=new double[n_track];
	if(sim_rand==NULL){puts("Out of memory...");return -1;}	
	*/
	for(i=0;i<n_track;i++)sim_real[i]=10000;
	
	double av[MATLEN], sd[MATLEN];
	/*
	double *av, *sd, *zs, *pv;
	av=new double[n_track];
	if(av==NULL){puts("Out of memory...");return -1;}		
	for(i=0;i<n_track;i++)av[i]=0;
	sd=new double[n_track];
	if(sd==NULL){puts("Out of memory...");return -1;}	
	for(i=0;i<n_track;i++)sd[i]=0;
	zs=new double[n_track];
	if(zs==NULL){puts("Out of memory...");return -1;}		
	pv=new double[n_track];
	if(pv==NULL){puts("Out of memory...");return -1;}		
	*/
	int sta_sh[MATLEN], sta_lo[MATLEN], ori_select[MATLEN], best_not_ali[MATLEN];
	/*int *sta_sh, *sta_lo, *ori_select, *best_not_ali;
	sta_sh=new int[n_track];
	if(sta_sh==NULL){puts("Out of memory...");return -1;}		
	sta_lo=new int[n_track];
	if(sta_lo==NULL){puts("Out of memory...");return -1;}			
	ori_select=new int[n_track];
	if(ori_select==NULL){puts("Out of memory...");return -1;}		
	best_not_ali=new int[n_track];
	if(best_not_ali==NULL){puts("Out of memory...");return -1;}	
	for(i=0;i<n_track;i++)sta_sh[i]=sta_lo[i]=ori_select[i]=best_not_ali[i]=0;
	*/
	int test, measure_type;	
	double pv_min_thr[4]={0.05, 0.05, 0.05, 0.05};
	double pv_best=1;
	for(test=0;test<4;test++)
	{
		measure_type=test%2;
		int n_cycle;
		if(test<2)n_cycle=n_cycle_small;
		else n_cycle=n_cycle_large;
		//real two orientations
		measure(pfm[k0],pfm[k1],olen_min,olen_max,OLIGNUM,granul_back,over_scores_real[0],measure_type);
		measure(pfmc[k0],pfm[k1],olen_min,olen_max,OLIGNUM,granul_back,over_scores_real[1],measure_type);
		for(i=0;i<n_track;i++)sim_real[i]=-10000;	
		for(j=0;j<n_track;j++)
		{
			for(ori=0;ori<2;ori++)
			{
				int sta_sh_cur=0, sta_lo_cur=0;
				double sim_score=Min_trackReal(over_scores_real[ori],olen_min,olen_max,j+overlap_min,sta_sh_cur,sta_lo_cur);
				if(sim_score>sim_real[j])
				{
					sim_real[j]=sim_score;
					sta_sh[j]=sta_sh_cur;
					sta_lo[j]=sta_lo_cur;
					ori_select[j]=ori;
				}
			}
		}
		for(k=0;k<2;k++)
		{
			for(i=0;i<olen[k];i++)
			{						
				for(j=0;j<OLIGNUM;j++)
				{				
					pfmr[k][i][j]=pfm[k][i][j];
				}
			}
			for(n=1;n<=olen[k];n++)Permute_columns(pfmr[k],olen[k],OLIGNUM,granul_back);
		}
		for(i=0;i<n_track;i++)av[i]=sd[i]=0;
		for(n=1;n<=n_cycle;n++)
		{
			if (n % reini_permut == 0)
			{
				for (k = 0; k < 2; k++)
				{
					for (i = 0; i < olen[k]; i++)
					{
						for (j = 0; j < OLIGNUM; j++)
						{
							pfmr[k][i][j] = pfm[k][i][j];
						}
					}					
				}
			}
			for(k=0;k<2;k++)Permute_columns(pfmr[k],olen[k],OLIGNUM,granul_back);		
			if(k1==1) //vtoraya dlinnee k0=0 k1=1
			{						
				measure(pfm[0],pfmr[1],olen_min,olen_max,OLIGNUM,granul_back,over_scores_rand_one[0],measure_type);//2nd matrix permute		
				measure(pfmc[0],pfmr[1],olen_min,olen_max,OLIGNUM,granul_back,over_scores_rand_one[1],measure_type);//2nd matrix permute		
				measure(pfmr[0],pfm[1],olen_min,olen_max,OLIGNUM,granul_back,over_scores_rand_two[0],measure_type);//1st matrix permute
				measure(pfmr[0],pfmc[1],olen_min,olen_max,OLIGNUM,granul_back,over_scores_rand_two[1],measure_type);//1st matrix permute	
			}		
			else//pervaya dlinnee k0=1 k1=0
			{			
				measure(pfmr[1],pfm[0],olen_min,olen_max,OLIGNUM,granul_back,over_scores_rand_one[0],measure_type);//2nd matrix permute		
				measure(pfmr[1],pfmc[0],olen_min,olen_max,OLIGNUM,granul_back,over_scores_rand_one[1],measure_type);//2nd matrix permute		
				measure(pfm[1],pfmr[0],olen_min,olen_max,OLIGNUM,granul_back,over_scores_rand_two[0],measure_type);//1st matrix permute
				measure(pfmc[1],pfmr[0],olen_min,olen_max,OLIGNUM,granul_back,over_scores_rand_two[1],measure_type);//1st matrix permute
			}		
			for(ori=0;ori<2;ori++)
			{
				for(j=0;j<n_track;j++)sim_rand[j]=-10000;		
				for(j=0;j<n_track;j++)
				{
					double sim_score=Min_trackRand(over_scores_rand_one[ori],olen_min,olen_max,j+overlap_min);
					if(sim_score>sim_rand[j])sim_rand[j]=sim_score;
					sim_score=Min_trackRand(over_scores_rand_two[ori],olen_min,olen_max,j+overlap_min);
					if(sim_score>sim_rand[j])sim_rand[j]=sim_score;
				}
				CellCount(over_scores_rand_one[ori],olen_min,olen_max,av_cell[ori],sd_cell[ori]);
				CellCount(over_scores_rand_two[ori],olen_min,olen_max,av_cell[ori],sd_cell[ori]);
				for(j=0;j<n_track;j++)av[j]+=sim_rand[j];
				for(j=0;j<n_track;j++)sd[j]+=sim_rand[j]*sim_rand[j];
			}
		}
		int n_cycle2=2*n_cycle;
		for(ori=0;ori<2;ori++)
		{
			for(i=0;i<olen_min;i++)
			{
				for(j=0;j<olen_max;j++)
				{
					av_cell[ori][i][j]/=n_cycle2;
					sd_cell[ori][i][j]=(sd_cell[ori][i][j]-n_cycle2*av_cell[ori][i][j]*av_cell[ori][i][j])/(n_cycle2-1);//dispersion out of alignment
				}
			}
		}		
		for(j=0;j<n_track;j++)av[j]/=n_cycle2;
		for(j=0;j<n_track;j++)
		{
			sd[j]=(sd[j]-n_cycle2*av[j]*av[j])/(n_cycle2-1);//dispersion in alignment
			//dobavka ot rayona v kototom net alignment v 5'
			int sta_sh_cur=sta_sh[j], sta_lo_cur=sta_lo[j];
			while(sta_sh_cur>0 && sta_lo_cur>0)
			{
				sta_sh_cur--;
				sta_lo_cur--;
				sd[j]+=sd_cell[ori_select[j]][sta_sh_cur][sta_lo_cur];
				best_not_ali[j]++;
			}
			sta_sh_cur=sta_sh[j]+j+overlap_min, sta_lo_cur=sta_lo[j]+j+overlap_min;
			while(sta_sh_cur<olen_min && sta_lo_cur<olen_max)
			{
				sd[j]+=sd_cell[ori_select[j]][sta_sh_cur][sta_lo_cur];
				sta_sh_cur++;
				sta_lo_cur++;
				best_not_ali[j]++;
			}
		}
		double zmax=0;
		for(j=0;j<n_track;j++)
		{
			sd[j]=sqrt(sd[j]);
			double zs=(sim_real[j]-av[j])/sd[j];
			if(zs>zmax)zmax=zs;			
		}
		double pv, z_thr=1.644805;//pv=0.1;
		if(zmax<z_thr)pv=1;
		else pv=Laplace(zmax);
		if(test<=1)pval[test]=pv;
		if(test>=2)pval[test-2]=pv;
		if(test==1)
		{
			pv_best=Max(pval[0],pval[1]);
			if(pval[0]>pv_min_thr[0] && pval[1]>pv_min_thr[1])
			{
				break;		
			}
		}
		if(test==3)
		{
			pv_best=Max(pval[0],pval[1]);							
		}		
	}
/*
	fprintf(out,"RealAv");
	for(j=0;j<n_track;j++)fprintf(out,"\t%f",sim_real[j]);	
	fprintf(out,"\nRandAv");
	for(j=0;j<n_track;j++)fprintf(out,"\t%f",av[j]);
	fprintf(out,"\nStDev");
	for(j=0;j<n_track;j++)fprintf(out,"\t%f",sd[j]);
	fprintf(out,"\nZscore");
	for(j=0;j<n_track;j++)fprintf(out,"\t%f",zs[j]);
	fprintf(out,"\nP-value");
	for(j=0;j<n_track;j++)fprintf(out,"\t%g",pv[j]);
	fprintf(out,"\n");
	*/
	//fprintf(out,"Best p-value\t%g\tMoifConst\t%d\tMoifPermut\t%d\tOffset\t%d\tStrand\t%s\n",pv_best,const_motif,permuted_motif,offset_best,cepi[ori_best]);	
	//char metrika[2][10]={"Euclid","Pearson"};
//	fprintf(out,"%s\t%s\t%s\t%d\t%d\t%d\t%d\t%g\t\t",metrika[measure_type],head[0],head[1],overlap_min,overlap_max,best_ali,best_not_ali_select,pv_min);	
//	fprintf(out,"\t");
//	for(j=0;j<n_track;j++)fprintf(out,"\t%g",pv[j]);
//	fprintf(out,"\n");
//	fclose(out);
	/*
	for(k=0;k<2;k++)
	{
		for(i=0;i<olen[k];i++)delete [] pfm[k][i];
		delete [] pfm[k];
		for(i=0;i<olen[k];i++)delete [] pfmr[k][i];
		delete [] pfmr[k];
		for(i=0;i<olen[k];i++)delete [] pfmc[k][i];
		delete [] pfmc[k];
	}
	delete [] pfm;
	delete [] pfmc;
	delete [] pfmr;
	for(k=0;k<2;k++)
	{
		for(i=0;i<olen_min;i++)
		{			
			delete [] over_scores_real[k][i];			
		}	
		delete [] over_scores_real[k];
	}
	delete [] over_scores_real;
	for(k=0;k<2;k++)
	{
		for(i=0;i<olen_min;i++)
		{			
			delete [] over_scores_rand_one[k][i];			
		}	
		delete [] over_scores_rand_one[k];
	}
	delete [] over_scores_rand_one;
	for(k=0;k<2;k++)
	{
		for(i=0;i<olen_min;i++)
		{			
			delete [] over_scores_rand_two[k][i];			
		}	
		delete [] over_scores_rand_two[k];
	}
	delete [] over_scores_rand_two;	
	for(k=0;k<2;k++)
	{
		for(i=0;i<olen_min;i++)
		{			
			delete [] sd_cell[k][i];			
		}	
		delete [] sd_cell[k];
	}
	delete [] sd_cell;
	for(k=0;k<2;k++)
	{
		for(i=0;i<olen_min;i++)
		{			
			delete [] av_cell[k][i];			
		}	
		delete [] av_cell[k];
	}
	delete [] av_cell;*/
	/*
	delete [] sta_sh;
	delete [] sta_lo;
	delete [] best_not_ali;
	delete[] av;
	delete[] sd;
	delete [] sim_real;
	delete[] pv;
	delete[] zs;
	delete [] sim_rand;*/	
	return pv_best;
}
