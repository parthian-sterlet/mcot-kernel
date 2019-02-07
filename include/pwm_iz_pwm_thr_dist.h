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
int pwm_iz_pwm_thr_dist(double pwm_source[][OLIGNUM], int lenp, char *file_pro, int n_thr_touzet, double *thr_touzet, double *fp_rate, char *partner_db)
{
	int i, j, n, nseq_pro, len_pro; 
	char dp[SEQLEN];	
	double p;
	FILE *in;

	if(strstr(partner_db,"hs")!=NULL){len_pro=2000;nseq_pro=19795;}
	else 
	{
		if(strstr(partner_db,"mm")!=NULL){len_pro=2000;nseq_pro=19991;}
		else 
		{
			if(strstr(partner_db,"dapseq")!=NULL){len_pro=1500;nseq_pro=27202;}
			else
			{
				if(strstr(partner_db,"moss")!=NULL){len_pro=1500;nseq_pro=18672;}
				else
				{
					printf("Partner databasee %s is wrong\n",partner_db);
					return -1;
				}
			}
		}
	}	            	
   int compl2[2]={0,1};
   int compl1;
   //if(compl_==0)compl2[1]=-1;
   //if(compl_==1)compl2[0]=-1;	

	int nseq=0;
	int len1=0;	
	int word=1;               
	len1=lenp+word-1;//dlina vyborki obu4eniya
	int ten[6]={1,4,16,64,256,1024};	
	double min=0, raz=0;
	PWMScore(min,raz,lenp,pwm_source);
	double min0=min, raz0=raz;
    int cod[SEQLEN];	
	for(i=0;i<n_thr_touzet;i++)fp_rate[i]=0;
	double thr0=thr_touzet[n_thr_touzet-1];//min				
	if((in=fopen(file_pro,"rt"))==NULL)
	{
		   printf("Input file %s can't be opened!",file_pro);
		   return -1;
	}
	double all_pos=0;
	int n_thr_touzet1=n_thr_touzet-1;
	for(n=0;n<nseq_pro;n++)
	{	    		
	   	fgets(dp,len_pro+2,in);	   	
	    DelChar(dp,'\n');					
		int len_pro1=strlen(dp);		
		int len21=len_pro1-len1;				
		for(compl1=0;compl1<2;compl1++)
		{		
			if(compl2[compl1]==-1)continue;
			if(compl2[compl1]==1) if(ComplStr(dp)!=1)  {puts("Out of memory...");return -1;}		
			char d2[SEQLEN];							
			p=-1000;
			for(i=0;i<=len21;i++)
			{	
				strncpy(d2,&dp[i],len1);
				d2[len1]='\0';
				GetSostPro(d2,word,cod);
				if(strstr(d2,"n")!=NULL){continue;}
				all_pos++;
				  double score=0;
    			  for(j=0;j<lenp;j++)
				  {
					 score+=pwm_source[j][cod[j]];
				  }
				  score-=min0;
				  score/=raz0;			
				  if(score>=thr0)			  
				  {
						for(j=n_thr_touzet1;j>=0;j--)
						{
							if(score>=thr_touzet[j])fp_rate[j]++;
							else break;
						}
				  }			 
			}						
		}							
	}
	fclose(in);	
	for(j=0;j<n_thr_touzet;j++)fp_rate[j]/=all_pos;	
	return 1;
}