int compare_qq( const void *X1, const void *X2 )
{
	int X=( *(int*)X1 - *(int*)X2 );
	if(X>0)return 1;
	if(X<0)return -1;		
	return 0;
}
int throw_predictions(int *peak_len, profile *anc, profile *par, int len_a, int len_p, int rand_prof, int *throw_err_real, 
					  int nseq_real, int nseq_rand, char **seq, int height, char *file_err)
{	
	char peak[SEQLEN], spacer[SEQLEN];// dna sequence of peak
	char peak_anc[SEQLEN], peak_par[SEQLEN];
    int  i, j, k, n;
	int ret_value=0;	
	
//	puts("Syntax: 1fasta_length 2ipre1 3ipre2 4insit1 5insit2 6int mat1 7int mat2 8opre1 9opre2 10int height ");
//	puts("11int size_min 12int size_max 13int rand_profile(0 0.5/0.5 1 = 1st 2 = 2nd 14file error_list_real 15file error_list_rand");
			
	FILE *out_err;		

	int test=0;
	for(n=0;n<nseq_rand;n++)
	{
		if(anc->nsit[n]>0 && par->nsit[n]>0)test++;
	}
	if(test==0)return 0;
	//int over=0, tot=0;
	int rand_ini=rand()%2;
	for(n=0;n<nseq_rand;n++)
	{
		if(anc->nsit[n]==0 || par->nsit[n]==0)continue;
		int n_clust[2], *nsit_clust[2], *clust_len[2], clust_tot_len[2];
		int n_real=n/height;// nomer realnoy posl-ti
		int len[2], nsit[2];
		memset(peak_anc,'n',peak_len[n]);
		peak_anc[peak_len[n]]='\0';							
		memset(peak_par,'n',peak_len[n]);
		peak_par[peak_len[n]]='\0';							
		len[0]=len_a, len[1]=len_p;
		nsit[0]=anc->nsit[n], nsit[1]=par->nsit[n];		
//		if((n+1)%height==0){over=0;tot=0;}
		for(i=0;i<2;i++)
		{			
			if(i==0)
			{
				//printf("Anchor sites:\n");
				for(j=0;j<anc->nsit[n];j++)
				{
	//				printf("%d%c ",anc->sta[n][j],anc->cep[n][j]);
					for(k=anc->sta[n][j];k<anc->sta[n][j]+len_a;k++)
					{
						peak_anc[k]=seq[n_real][k];
		//				printf("%c",peak_anc[k]);						
					}					
			//		printf("\n");
				}			
			}
			else
			{
				//printf("Partner sites:\n");
				for(j=0;j<par->nsit[n];j++)
				{
					//printf("%d%c ",par->sta[n][j],par->cep[n][j]);
					for(k=par->sta[n][j];k<par->sta[n][j]+len_p;k++)
					{
						peak_par[k]=seq[n_real][k];
						//printf("%c",peak_par[k]);						
					}
				}	
				//printf("\n");
			}
		}
		/*for(i=0;i<2;i++)
		{
			printf("Peak anchor\n");
			for(k=0;k<peak_len[n];k++)
			{
				printf("%c",peak_anc[k]);						
				if((k+1)%120==0)printf("\n");
			}
			printf("\n");
			printf("Peak partner\n");
			for(k=0;k<peak_len[n];k++)
			{
				printf("%c",peak_par[k]);						
				if((k+1)%120==0)printf("\n");
			}
			printf("\n");
		}*/
		for(i=0;i<2;i++)
		{		
			n_clust[i]=1;
			for(j=1;j<nsit[i];j++)
			{
				int dif;
				if(i==0)dif=anc->sta[n][j]-anc->sta[n][j-1];
				else dif=par->sta[n][j]-par->sta[n][j-1];
				if(dif>=len[i])n_clust[i]++;
			}
			clust_len[i] = new int[n_clust[i]];
			if(clust_len[i]==NULL ){puts("Out of memory...");return -1;}
			nsit_clust[i] = new int[n_clust[i]];
			if(nsit_clust[i]==NULL ){puts("Out of memory...");return -1;}
			int clust_cur=0;
			nsit_clust[i][0]=1;
			for(j=1;j<n_clust[i];j++)nsit_clust[i][j]=0;			
			for(j=1;j<nsit[i];j++)
			{
				int dif;
				if(i==0)dif=anc->sta[n][j]-anc->sta[n][j-1];
				else dif=par->sta[n][j]-par->sta[n][j-1];
				if(dif<len[i])nsit_clust[i][clust_cur]++;
				else 
				{
					clust_cur++;
					nsit_clust[i][clust_cur]++;
				}
			}			
			int nsit_cur=0;
			clust_tot_len[i]=0;
			for(j=0;j<n_clust[i];j++)
			{
				if(i==0)clust_len[i][j]=len[i]+anc->sta[n][nsit_cur+nsit_clust[i][j]-1]-anc->sta[n][nsit_cur];
				else clust_len[i][j]=len[i]+par->sta[n][nsit_cur+nsit_clust[i][j]-1]-par->sta[n][nsit_cur];
				clust_tot_len[i]+=clust_len[i][j];
				nsit_cur+=nsit_clust[i][j];
			}
		}		
		char ***clust_seq;
		clust_seq = new char**[2];				
		if(clust_seq==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<2;i++)
		{
			clust_seq[i] = new char*[n_clust[i]];
			if(clust_seq[i]==NULL ){puts("Out of memory...");return -1;}
			int nsit_cur=0;
			for(j=0;j<n_clust[i];j++)
			{
				clust_seq[i][j] = new char[clust_len[i][j]+1];
				if(clust_seq[i][j]==NULL ){puts("Out of memory...");return -1;}
				if(i==0)
				{
					for(k=0;k<clust_len[i][j];k++)clust_seq[i][j][k]=peak_anc[anc->sta[n][nsit_cur]+k];
				}
				else
				{
					for(k=0;k<clust_len[i][j];k++)clust_seq[i][j][k]=peak_par[par->sta[n][nsit_cur]+k];
				}
				clust_seq[i][j][clust_len[i][j]]='\0';
				nsit_cur+=nsit_clust[i][j];
			}
		}
		char ***clust_cep;
		clust_cep = new char**[2];				
		if(clust_cep==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<2;i++)
		{
			clust_cep[i] = new char*[n_clust[i]];
			if(clust_cep[i]==NULL ){puts("Out of memory...");return -1;}
			int nsit_cur=0;
			for(j=0;j<n_clust[i];j++)
			{
				clust_cep[i][j] = new char[nsit_clust[i][j]+1];
				if(clust_cep[i][j]==NULL ){puts("Out of memory...");return -1;}
				if(i==0)
				{
					for(k=0;k<nsit_clust[i][j];k++)clust_cep[i][j][k]=anc->cep[n][nsit_cur++];
				}
				else
				{
					for(k=0;k<nsit_clust[i][j];k++)clust_cep[i][j][k]=par->cep[n][nsit_cur++];
				}
				clust_cep[i][j][nsit_clust[i][j]]='\0';
				//nsit_cur+=nsit_clust[i][j];
			}
		}
		int ***cel_pv;
		cel_pv = new int**[2];				
		if(cel_pv==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<2;i++)
		{
			cel_pv[i] = new int*[n_clust[i]];
			if(cel_pv[i]==NULL ){puts("Out of memory...");return -1;}
			int nsit_cur=0;
			for(j=0;j<n_clust[i];j++)
			{
				cel_pv[i][j] = new int[nsit_clust[i][j]+1];
				if(cel_pv[i][j]==NULL ){puts("Out of memory...");return -1;}
				if(i==0)
				{
					for(k=0;k<nsit_clust[i][j];k++)
					{
						cel_pv[i][j][k]=anc->cel[n][nsit_cur++];
						//printf("I %d J %d K %d Cel %d\n",i,j,k,pv[i][j][k]);
					}
				}
				else
				{
					for(k=0;k<nsit_clust[i][j];k++)
					{
						cel_pv[i][j][k]=par->cel[n][nsit_cur++];
						//printf("I %d J %d K %d Cel %d\n",i,j,k,pv[i][j][k]);
					}
				}
				//nsit_cur+=nsit_clust[i][j];
			}
		}
		double ***pv;
		pv = new double**[2];				
		if(pv==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<2;i++)
		{
			pv[i] = new double*[n_clust[i]];
			if(pv[i]==NULL ){puts("Out of memory...");return -1;}
			int nsit_cur=0;
			for(j=0;j<n_clust[i];j++)
			{
				pv[i][j] = new double[nsit_clust[i][j]+1];
				if(pv[i][j]==NULL ){puts("Out of memory...");return -1;}
				if(i==0)
				{
					for(k=0;k<nsit_clust[i][j];k++)
					{
						pv[i][j][k]=anc->pv[n][nsit_cur++];
						//printf("I %d J %d K %d Cel %d\n",i,j,k,pv[i][j][k]);
					}
				}
				else
				{
					for(k=0;k<nsit_clust[i][j];k++)
					{
						pv[i][j][k]=par->pv[n][nsit_cur++];
						//printf("I %d J %d K %d Cel %d\n",i,j,k,pv[i][j][k]);
					}
				}
				//nsit_cur+=nsit_clust[i][j];
			}
		}		
		int ***shift;
		shift = new int**[2];				
		if(shift==NULL){puts("Out of memory...");return -1;}
		for(i=0;i<2;i++)
		{
			shift[i] = new int*[n_clust[i]];
			if(shift[i]==NULL ){puts("Out of memory...");return -1;}
			int nsit_cur=0;
			for(j=0;j<n_clust[i];j++)
			{
				shift[i][j] = new int[nsit_clust[i][j]];
				if(shift[i][j]==NULL ){puts("Out of memory...");return -1;}
				if(i==0)
				{
					for(k=0;k<nsit_clust[i][j];k++)shift[i][j][k]=anc->sta[n][nsit_cur+k]-anc->sta[n][nsit_cur];
				}
				else
				{
					for(k=0;k<nsit_clust[i][j];k++)shift[i][j][k]=par->sta[n][nsit_cur+k]-par->sta[n][nsit_cur];
				}
				nsit_cur+=nsit_clust[i][j];
			}
		}
		if(clust_tot_len[0]==peak_len[n] || clust_tot_len[1]==peak_len[n])
		{						
			if ((out_err = fopen(file_err, "at")) == NULL)
			{
				printf("Input file %s can't be opened!\n", file_err);
				return -1;
			}							
			fprintf(out_err,"Too many sites! Mot(Thr) %d(%d) or Mot %d(%d) motifs Seq(Nsites1,2) %d (%d,%d)\n", anc->mot, anc->nam,par->mot, par->nam,n+1,anc->nsit[n], par->nsit[n]);
			throw_err_real[n_real]=1;
			ret_value++;	
			fclose(out_err);
			for(i=0;i<2;i++)
			{
				for(j=0;j<n_clust[i];j++)delete[] clust_seq[i][j]; 
				for(j=0;j<n_clust[i];j++)delete[] clust_cep[i][j]; 
				for(j=0;j<n_clust[i];j++)delete[] cel_pv[i][j]; 
				for(j=0;j<n_clust[i];j++)delete[] pv[i][j]; 
				for(j=0;j<n_clust[i];j++)delete[] shift[i][j]; 
				delete[] clust_seq[i];
				delete[] clust_cep[i];
				delete[] cel_pv[i];
				delete[] pv[i];
				delete[] shift[i];
			}
			delete [] clust_seq;
			delete [] clust_cep;
			delete [] cel_pv;
			delete [] pv;
			delete [] shift;
			for(i=0;i<2;i++)delete [] nsit_clust[i];
			for(i=0;i<2;i++)delete [] clust_len[i];			
			continue;
		}
		//int exit_h=0;		
		int spacer_lenr[SEQLEN];		
		int n_clust_max=Max(n_clust[0],n_clust[1]);
		n_clust_max++;
		//spacer_lenr = new int[n_clust_max];
		//if(spacer_lenr==NULL ){puts("Out of memory...");return -1;}
		int order[2][SEQLEN];
		//int **order;		
		//order = new int*[2];				
		//if(order==NULL){puts("Out of memory...");return -1;}
		/*
		for(j=0;j<2;j++)
		{			
			order[j] = new int[n_clust[j]];
			if(order[j]==NULL){puts("Out of memory...");return -1;}
		}*/
		//for(h=0;h<height;h++)					
		//printf("%d,%d\t",nsit[0],nsit[1]);
		//if(exit_h==1)break;
		int implant, permut;
		if(rand_prof==0)
		{
			//implant = rand()%2;//randomly chosen
			int rand_here=n%2;
			if(rand_here==rand_ini)implant=0;
			else implant=1;
			permut=1-implant;
		}
		else
		{
			if(rand_prof==1)implant=1, permut=0;//permutation for anchor
			if(rand_prof==2)implant=0, permut=1;//permutation for partner
		}
		//printf("Impl %d Perm %d\t",implant,permut);
		int cycle=0;
		int cycle_max=10000;
		int gom=0;
		for(j=0;j<n_clust[permut];j++)order[permut][j]=j;
		do
		{
		//	if(cycle!=0)printf("Cycle %d\n",cycle+1);
			BigMix1(order[permut],n_clust[permut]);
			int spacer_len_tot=peak_len[n]-clust_tot_len[permut];
			//strcpy(peak,peak0[implant]);
		//	for(j=0;j<nsit[permut];j++)pre[permut][j].get_copy(&pra[permut][order[permut][j]]);							
			//for(j=0;j<nsit[permut];j++)pre[permut][j].get_copy(&pra[permut][j]);										
			for(j=0;j<n_clust[permut];j++)
			{
				spacer_lenr[j]=rand()%spacer_len_tot;
			}				
			qsort(spacer_lenr,n_clust[permut],sizeof(int),compare_qq);	
			spacer_lenr[n_clust[permut]]=spacer_len_tot;
			for(j=n_clust[permut];j>0;j--)
			{
				spacer_lenr[j]=spacer_lenr[j]-spacer_lenr[j-1];
			}			
			memset(peak,'n',spacer_lenr[0]);			
			peak[spacer_lenr[0]]='\0';
			int n_cur=spacer_lenr[0];
			//printf("SpacLen perm: ");
			//for(j=0;j<n_clust[permut];j++)printf("%d ",spacer_lenr[j]);
			//printf("\t");
			for(j=0;j<n_clust[permut];j++)
			{
				int j1=order[permut][j];					
			//	for(k=0;k<clust_len[permut][j1];k++)peak[m++]=clust_seq[permut][j1][k];
				strcat(peak,clust_seq[permut][j1]);
				n_cur+=clust_len[permut][j1];
				memset(spacer,'n',spacer_lenr[j+1]);
				spacer[spacer_lenr[j+1]]='\0';
				strcat(peak,spacer);
				n_cur+=spacer_lenr[j+1];
			}
			peak[peak_len[n]]='\0';
			/*printf("Permuted peak:\n");
			for(k=0;k<peak_len[n];k++)
			{
				printf("%c",peak[k]);						
				if((k+1)%120==0)printf("\n");
			}
			printf("\n");*/
			gom=0;
			//for(j=0;j<peak_len;j++)
			int nsit_cur=0;
			for(j=0;j<n_clust[implant];j++)
			{		
				for(k=0;k<clust_len[implant][j];k++)
				{
					int pos;
					char y;
					if(implant==1)
					{
						pos=par->sta[n][nsit_cur]+k;
						y = peak_par[pos];					
					}
					else 
					{
						pos=anc->sta[n][nsit_cur]+k;
						y = peak_anc[pos];					
					}
					if(y == 'n')continue;
					char x = peak[pos];
					if(x == 'n')continue;						
					if (x != y)
					{
						gom=1;
						break;
					}
				}
				nsit_cur+=nsit_clust[implant][j];				
				if(gom==1)break;
			}				
			if(gom==1)
			{
				cycle++;
//				printf("No homology! Nseq %d Cycle %d\n",n, j,cycle);
				if(cycle>=cycle_max)
				{						
					if ((out_err = fopen(file_err, "at")) == NULL)
					{
						printf("Input file %s can't be opened!\n", file_err);
						return -1;
					}							
					if(permut==0)fprintf(out_err,"Unable to throw permute Mot(Thr) %d(%d) motif Seq(Nsites permute,implant) %d (%d,%d)\n", anc->mot, anc->nam,n/height+1,anc->nsit[n], par->nsit[n]);
					else fprintf(out_err,"Unable to throw permute Mot(Thr) %d(%d) motif Seq(Nsites permute,implant) %d (%d,%d)\n", par->mot, par->nam,n/height+1,par->nsit[n], anc->nsit[n]);
					throw_err_real[n_real]=1;
				//	exit_h=1;
					ret_value++;							
					fclose(out_err);
					break;
				//	return -1;
				}
			//	printf ("Nseq %d Cycle %d\n",n+1,cycle);
				continue;
			}
	/*		else 
			{
				//printf("Success!\n");
				if(h==height-1)
				{
					int yy=1;
				}
			}*/
			nsit_cur=0;
			int m=spacer_lenr[0];
			for(j=0;j<n_clust[permut];j++)
			{
				int j1=order[permut][j];
				//printf("Perm ");
				//if(permut==0)printf("anc ");
				//else printf("part ");
				//printf("sites:\t");
				for(k=0;k<nsit_clust[permut][j1];k++)
				{
					if(permut==0)
					{
						anc->sta[n][nsit_cur]=m+shift[permut][j1][k];//prf[implant]->sta[n][nsit_cur]
						anc->cep[n][nsit_cur]=clust_cep[permut][j1][k];
						anc->cel[n][nsit_cur]=cel_pv[permut][j1][k];
						anc->pv[n][nsit_cur]=pv[permut][j1][k];
						//printf("I %d J %d K %d Cel %d\n",permut,j1,k,pv[permut][j1][k]);
					//	printf("%d%c ",anc->sta[n][nsit_cur],anc->cep[n][nsit_cur]);
					}
					else
					{
						par->sta[n][nsit_cur]=m+shift[permut][j1][k];//prf[implant]->sta[n][nsit_cur]
						par->cep[n][nsit_cur]=clust_cep[permut][j1][k];
						par->cel[n][nsit_cur]=cel_pv[permut][j1][k];
						par->pv[n][nsit_cur]=pv[permut][j1][k];
					//	printf("I %d J %d K %d Cel %d\n",permut,j1,k,pv[permut][j1][k]);
						//printf("%d%c ",par->sta[n][nsit_cur],par->cep[n][nsit_cur]);
					}
			//		pra.sco=1;
					//int t;
				//	for(t=0;t<motif_len[permut];t++)pra.seq[t]=clust_seq[permut][j1][t+shift[permut][j1][k]];
				//	pra.seq[motif_len[permut]]='\0';
				//	if(clust_cep[permut][j1][k]=='-')ComplStr(pra.seq);
					nsit_cur++;
				}
				//printf("\n");
				//pra.position();
			//	fprintf(out[permut], "%d\t%f\t%c\t%s\n", pra.pos, pra.sco, pra.cep, pra.seq);
				m+=clust_len[permut][j];
				m+=spacer_lenr[j+1];
			}
		/*	if(gom==0)
			{
				int x,y;
				tot++;
				int n_over=0;
				for(x=0;x<anc->nsit[n];x++)
				{
					int anc1=anc->sta[n][x];
					for(y=0;y<par->nsit[n];y++)
					{						
						int par1=par->sta[n][y];						
						if((par1<=anc1 && anc1-par1<len_p) || (anc1<par1 && par1-anc1<len_a)) 
						{
							n_over++;
							printf("Nseq %d Overlap %d out of %d\tPar %d Anc %d\n",n+1,n_over,tot,par1,anc1);
						}
					}
				}
				if(n_over>0)
				{					
					over++;
					printf("Nseq %d Tot %d Overlap %d\n",n+1,tot,over);
				}
			}*/			
		}
		while(gom==1);						
		//printf("%s\n",peak);		
		for(i=0;i<2;i++)
		{
			for(j=0;j<n_clust[i];j++)delete[] clust_seq[i][j]; 
			for(j=0;j<n_clust[i];j++)delete[] clust_cep[i][j]; 
			for(j=0;j<n_clust[i];j++)delete[] cel_pv[i][j]; 
			for(j=0;j<n_clust[i];j++)delete[] pv[i][j]; 
			for(j=0;j<n_clust[i];j++)delete[] shift[i][j]; 
			delete[] clust_seq[i];
			delete[] clust_cep[i];
			delete[] cel_pv[i];
			delete[] pv[i];
			delete[] shift[i];
			//delete [] order[i];//printf("Order\n");
		}
		delete [] clust_seq;
		delete [] clust_cep;
		delete [] cel_pv;
		delete [] pv;
		delete [] shift;
		//delete [] spacer_lenr;
		//delete [] order;
		for(i=0;i<2;i++)delete [] nsit_clust[i];//printf("Nsit_clust\n");		
		for(i=0;i<2;i++)delete [] clust_len[i];//printf("Clust_len\n");				
	}
	//if(ret_value>0)
	return ret_value;
}