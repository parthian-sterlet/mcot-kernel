int pwm_rec0(matrices *mat, double thr, int len_pro, int nseq_pro, char ***seq, profile *real, int &all_pos)
{
	int i, j, n; 			
    int compl1;	
	int word=1;	
	int len1=mat->len+word-1;//dlina vyborki obu4eniya	
	int cod[SEQLEN];
	//int all_pos=0;
//	int rec_pos=0;
	int nseq_rec=0;
	for(n=0;n<nseq_pro;n++)
	{		
//		if ((n + 1) % 100 == 0)printf("%d ", n + 1);
		real->nsit[n]=0;
		int len_pro1=strlen(seq[0][n]);				
		int len21=len_pro1-len1;				
		for(compl1=0;compl1<2;compl1++)
		{		
			//if(compl1==1) if(ComplStr(d)!=1)  {puts("Complem error...");return -1;}		
			char d2[MATLEN];										
			for(i=0;i<=len21;i++)
			{	
				strncpy(d2,&seq[compl1][n][i],len1);
				d2[len1]='\0';
				if(strstr(d2,"n")!=NULL){continue;}
				GetSostPro(d2,word,cod);								
				all_pos++;
				  double score=0;
    			  for(j=0;j<len1;j++)
				  {
					 score+=mat->wei[j][cod[j]];
				  }
				  score-=mat->min;
				  score/=mat->raz;			
				  if(score>=thr)			  
				  {
	//				  rec_pos++;
					  real->nsit[n]++;
		//			  printf("Seq %d, StaPos %d_Compl %d Sco %f ->Interval %d Site %d\n",n+1,i+1,compl1,score,k+1,real_sites[k][n]);						
				  }			 
			}						
		}					
	}	
	/*for(j=0;j<num_thr;j++)
	{
		printf("Internal Num_thr %d\t",j);
		for(n=0;n<nseq_pro;n++)
		{
			printf("\t%d",real_sites[j][n]);
			//for(m=0;m<real_sites[j][n];m++)
		}
		printf("\n");
	}*/
	all_pos/=2;
	for(n=0;n<nseq_pro;n++)if(real->nsit[n]>0)nseq_rec++;
	/*FILE *out_stat;
	if(mot==0)
	{
		if((out_stat=fopen("rec_pos.txt","wt"))==NULL)    
		{
			printf("Input file can't be opened!\n");
			return -1;
		}
	}
	else
	{
		if((out_stat=fopen("rec_pos.txt","at"))==NULL)    
		{
			printf("Input file can't be opened!\n");
			return -1;
		}

	}
	fprintf(out_stat,"Motif %d\tThr %f\t%3f\t%d\t%d\t%.3e\t%d\t%d\t\n",mot,thr,(double)nseq_rec/nseq_pro,nseq_rec,nseq_pro,(double)rec_pos/all_pos,rec_pos,all_pos);
	fclose(out_stat);*/
	return 1;
}
int pwm_rec1(matrices *mat,double thr,int len_pro,int nseq_pro,char ***seq, profile *real)
{
	int i, j, n; 			 
    int compl1;	
	int word=1;	
	int len1=mat->len+word-1;//dlina vyborki obu4eniya	
	int cod[SEQLEN];		
	char cepx[3]="+-";			
	for(n=0;n<nseq_pro;n++)
	{		
	//	if ((n + 1) % 100 == 0)printf("%d ", n + 1);
		int len_pro1=strlen(seq[0][n]);				
		int len21=len_pro1-len1;							
		int x=0;
		for(i=0;i<=len21;i++)
		{												
			for(compl1=0;compl1<2;compl1++)			
			{	
				char d[MATLEN];
				int ista;
				if(compl1==0)ista=i;
				else ista=len21-i;
				strncpy(d,&seq[compl1][n][ista],len1);
				d[len1]='\0';
				if(strstr(d,"n")!=NULL){continue;}
				GetSostPro(d,word,cod);				
				double score=0;
    			for(j=0;j<len1;j++)
				{
					score+=mat->wei[j][cod[j]];
				}
				score-=mat->min;
				score/=mat->raz;			
				if(score>=thr)			  
				{
					real->sta[n][x]=i;							  
					real->cep[n][x]=cepx[compl1];
					real->sco[n][x]=score;
					x++;
				}								
			}							
		}		
	}
	//for(n=0;n<nseq_pro;n++)printf("%d ",nsit[n]);
	//printf("\n");
	return 1;
}
