int projoin(char *rera, char *motif,profile prf_a, profile prf_p, int shift_min, int shift_max, int len_a, int len_p, int *thr_pre_err,
			int nseq, char ***seq, result *sam, combi *hist, int *peak_len, asy_plot *plot, int &nseq_two_sites)
{	
	int n, k,j, x,y;
	char filebest[120], file_nsit_par[120], file_nsit_anc[120]; 

	FILE *outbest, *out_nsit_par, *out_nsit_anc;// *outhead,*out,  
//	FILE *out_cepi_seq;//*out_cepi_sit, *out_over_spac;//*outnsite_any,*out_hist,

	memset(filebest,'\0',sizeof(filebest));
	strcpy(filebest,rera);//real or random
	strcat(filebest,"_");
	// strcat(filebest,motif);//hocomoco or dapseq
	char buf[10];
	sprintf(buf,"%d",prf_a.mot);
	strcat(filebest,buf);
	sprintf(buf,"%d",prf_p.mot);
	strcat(filebest,buf);
	strcat(filebest,"_thr");
	sprintf(buf,"%d",prf_a.nam);
	strcat(filebest,buf);
	sprintf(buf,"%d",prf_p.nam);
	strcat(filebest,buf);
	strcpy(file_nsit_anc,filebest);
	strcpy(file_nsit_par,filebest);
	strcat(filebest,".best");
	strcat(file_nsit_anc,".anchor");
	strcat(file_nsit_par,".partner");	
	nseq_two_sites = 0;
		
	int *cepi_seq[4], *cepi_seq_dir, *cepi_seqa, *cepi_seqc;
	//int *cepi_sit[4];
	int *cepi_sit_one[4];	
	outbest=out_nsit_anc=out_nsit_par=NULL;
	int real=0;
	if(strstr(rera,"real")!=NULL)real=1;
	if(real==1)
	{
		if((outbest=fopen(filebest,"wt"))==NULL)
		{
 			printf("Input file %s can't be opened!\n",filebest);
			return -1;
		}
		fprintf(outbest,"#Seq\tA Start\tA End\tP Start\tP End\tMutual Loc\tLoc Type\tStrands\tMutual Ori\tA Score\tP Score\tA Seq\tP Seq\n");
		if((out_nsit_anc=fopen(file_nsit_anc,"wt"))==NULL)
		{
 			printf("Input file %s can't be opened!\n",file_nsit_anc);
			return -1;
		}		
		if((out_nsit_par=fopen(file_nsit_par,"wt"))==NULL)
		{
 			printf("Input file %s can't be opened!\n",file_nsit_par);
			return -1;
		}
	}
	int noverp=Min(len_p,len_a);//partial overlap
	noverp--;
	int noveri=1+abs(len_p-len_a)/2;// full overlap	
	int cepi_len=shift_max+noveri+noverp+1;		
	int nover=noveri+noverp;		
	for(j=0;j<4;j++)
	{
		cepi_seq[j]=new int[cepi_len];
		if(cepi_seq[j]==NULL){printf("Not  enough memory!");return -1;}	
		cepi_sit_one[j]=new int[cepi_len];
		if(cepi_sit_one[j]==NULL){printf("Not  enough memory!");return -1;}	
		for(k=0;k<cepi_len;k++)cepi_seq[j][k]=0;
	}
	cepi_seq_dir=new int[cepi_len];
	if(cepi_seq_dir==NULL){printf("Not  enough memory!");return -1;}	
	for(k=0;k<cepi_len;k++)cepi_seq_dir[k]=0;		
	cepi_seqa = new int[cepi_len]; //any ori
	if (cepi_seqa == NULL) { printf("Not  enough memory!"); return -1; }
	for (k = 0; k < cepi_len; k++)cepi_seqa[k] = 0;
	cepi_seqc = new int[cepi_len]; //cumulative for any
	if (cepi_seqc == NULL) { printf("Not  enough memory!"); return -1; }
	for (k = 0; k < cepi_len; k++)cepi_seqc[k] = 0;

	int nseq_rec_err[2]={0,0};//sequences that should be ignored in comparison of two profiles
	char oris[4][10]={"DirectAP","DirectPA","Inverted","Everted"};//direct_anchor_partner, direct partner_anchor, invert_anchor_partner, evert_anchor_partner
	char locs[3][10]={"Full","Partial","Spacer"};
	int dlen=0;
	{
		int dif_len=len_p-len_a;
		if(dif_len>=0)dlen=dif_len-(dif_len)/2;
		else dlen=dif_len/2;
	}
	int nseq_both = 0;//both motifs are present
	for(n=0;n<nseq;n++)
	{						
		int join_sites[NUM_THR][NUM_THR], join_sites_anc=0, join_sites_par=0;		
		int asym_anc=0,asym_par=0;
		if(prf_a.nsit[n]==0 || prf_p.nsit[n]==0)
		{		
			if(real==1)
			{
				fprintf(out_nsit_anc,"%d\n",asym_anc);
				fprintf(out_nsit_par,"%d\n",asym_par);
			}
			continue;
		}
	//	printf("%8d",n+1);
		int lenp=peak_len[n];
		int ce_full[NUM_THR][NUM_THR], ce_part[NUM_THR][NUM_THR], ce_spac[NUM_THR][NUM_THR];// no. of CE 
		for(j=0;j<NUM_THR;j++)for(k=0;k<NUM_THR;k++)join_sites[j][k]=ce_full[j][k]=ce_part[j][k]=ce_spac[j][k]=0;		
		int ce_full_anc = 0, ce_part_anc = 0, ce_spac_anc = 0;// no. of CE , anc
		int ce_full_par = 0, ce_part_par = 0, ce_spac_par = 0;// no. of CE , par
		int ce_full_eq = 0, ce_part_eq = 0, ce_spac_eq = 0;// no. of CE , equal
		if(thr_pre_err[n]==0)
		{				
			if (prf_a.nsit[n] >= 1 && prf_p.nsit[n] >= 1)nseq_both++;
//			printf("Nseq %d Anc %d Par %d\n",n+1,prf_a.nsit[n],prf_p.nsit[n]);
	//		for(y=0;y<prf_a.nsit[n];y++)printf("%d%c\t",prf_a.sta[n][y],prf_a.cep[n][y]);
		//	printf("\n");
		//	for(x=0;x<prf_p.nsit[n];x++)printf("%d%c\t",prf_p.sta[n][x],prf_p.cep[n][x]);
		//	printf("\n");
			for(k=0;k<cepi_len;k++)for(j=0;j<4;j++)cepi_sit_one[j][k]=0;																	
			for(x=0;x<prf_p.nsit[n];x++)
			{
				//printf("\b\b\b\b\b\b\b\b\b%8d %8d",x);
				int xsta,xend;// start & end position of motif				
				xsta=prf_p.sta[n][x];	
				xend=xsta+len_p-1;
			//	printf("x=%d\t%d\t%d\t\n",x,xsta,xend);
				int xsum=xsta+xend;
				int y0;			
				y0=prf_a.nsit[n]-1;												
				int ori_ce;
				double pv_par=prf_p.pv[n][x];
				int r_par=prf_p.cel[n][x];
				int p_karman;
				if(pv_par<plot->min)p_karman=0;
				else
				{
					if(pv_par>=plot->max)p_karman=plot->n_karman;
					else p_karman=1+(int)((pv_par-plot->min)/plot->inter);
				}
				for(y=0;y<prf_a.nsit[n];y++)
				{
					if(prf_a.mot==prf_p.mot)
					{
						if(x<=y)continue;//homodimer
					}																		
					double pv_anc=prf_a.pv[n][y];					
					int r_anc=prf_a.cel[n][y];
					join_sites[r_anc][r_par]=1;
					int a_karman;
					if(pv_anc<plot->min)a_karman=0;
					else
					{
						if(pv_anc>=plot->max)a_karman=plot->n_karman;
						else a_karman=1+(int)((pv_anc-plot->min)/plot->inter);
					}
					double r_dif=pv_anc-pv_par;
					if (r_dif > 0)join_sites_anc = 1;
					else join_sites_par=1;					
					int ysta,yend;// start & end position of motif
					ysta=prf_a.sta[n][y];												
					yend=ysta+len_a-1;												
	//				printf("1x=%d\t%d\t%d\t\t",x,xsta,xend);
		//			printf("1y=%d\t%d\t%d\t\n",y,ysta,yend);
					int ysum=ysta+yend;						
					int take_distance=0;	
					int partlial_len=-1;
					int full_len=-1;
					int spacer_len=-1;
					int cepi_pos;
					int cat[4]={0,0,0,0};
					if(xsta<=ysta && ysta<=xend)cat[0]=1;
					if(xsta<=yend && yend<=xend)cat[1]=1;
					if(ysta<=xsta && xsta<=yend)cat[2]=1;
					if(ysta<=xend && xend<=yend)cat[3]=1;
				//	if(sta_max-sta_min<shift_max || end_max-end_min<shift_max)take_distance=1;
					if(cat[0]==1 && cat[1]==1)
					{
						take_distance=1;//full overlap
						full_len=Min(ysta-xsta,xend-yend);
						sam->sit.full++;
						sam->sit.overlap++;
						plot->full[a_karman][p_karman]++;
						plot->overlap[a_karman][p_karman]++;
						if (r_dif > 0)
						{
							sam->anc_sit.full++;
							sam->anc_sit.overlap++;
						}
						else
						{
							sam->par_sit.full++;
							sam->par_sit.overlap++;
						}
					}
					else
					{
						if(cat[2]==1 && cat[3]==1)
						{
							take_distance=1;//full overlap
							full_len=Min(xsta-ysta,yend-xend);
							sam->sit.full++;
							sam->sit.overlap++;
							plot->full[a_karman][p_karman]++;
							plot->overlap[a_karman][p_karman]++;
							if (r_dif > 0)
							{
								sam->anc_sit.full++;
								sam->anc_sit.overlap++;
							}
							else
							{
								sam->par_sit.full++;
								sam->par_sit.overlap++;
							}
						}
					}
					if(take_distance==0)
					{
						if(cat[0]==1 && cat[3]==1)
						{
							take_distance=1;//partial overlap
							partlial_len=xend-ysta+1;
							sam->sit.overlap++;
							sam->sit.partial++;
							plot->partial[a_karman][p_karman]++;
							plot->overlap[a_karman][p_karman]++;
							if (r_dif > 0)
							{
								sam->anc_sit.partial++;
								sam->anc_sit.overlap++;
							}
							else
							{
								sam->par_sit.partial++;
								sam->par_sit.overlap++;
							}
						}
						else 
						{
							if(cat[1]==1 && cat[2]==1)
							{
								take_distance=1;//partial overlap
								partlial_len=yend-xsta+1;
								sam->sit.overlap++;
								sam->sit.partial++;
								plot->partial[a_karman][p_karman]++;
								plot->overlap[a_karman][p_karman]++;
								if (r_dif > 0)
								{
									sam->anc_sit.partial++;
									sam->anc_sit.overlap++;
								}
								else
								{
									sam->par_sit.partial++;
									sam->par_sit.overlap++;
								}
							}
						}
					}
					if(take_distance==0)//check for spacer
					{
						if(xend<ysta)spacer_len=ysta-xend-1;
						if(yend<xsta)spacer_len=xsta-yend-1;
						if (spacer_len >= shift_min && spacer_len <= shift_max)
						{
							take_distance = 1;//spacer
							sam->sit.spacer++;
							if (r_dif > 0)sam->anc_sit.spacer++;
							else sam->par_sit.spacer++;
							plot->spacer[a_karman][p_karman]++;
						}
					}					
					if(take_distance==1)
					{		
						sam->sit.any++;
						if (r_dif > 0)
						{
							sam->anc_sit.any++;
							asym_anc=1;
						}
						else 
						{
							sam->par_sit.any++;
							asym_par=1;
						}
						plot->any[a_karman][p_karman]++;
//						printf("2x=%d\t%d\t%d\t\t",x,xsta,xend);
//						printf("2y=%d\t%d\t%d\t\n",y,ysta,yend);
						int xsta0=ysta-dlen;
						int xend0=yend+dlen;
						int d0=0;							
						if(cepi_len%2==1)
						{
							xsta0--;xend0++;d0=1;
						}
						if(xsum==ysum)cepi_pos=0;
						else
						{
							if(ysum<xsum)cepi_pos=d0+xend-xend0;
							else cepi_pos=d0+xsta0-xsta;	
						}
						if(cepi_pos<0 || cepi_pos>cepi_len-1)
						{
							printf("Cepi pos error! %d\n", cepi_pos);
							exit(1);
						}
						if(prf_a.cep[n][y]=='+')
						{										
							if(prf_p.cep[n][x]=='+')
							{									
								if(ysum<=xsum)ori_ce=0;
								else ori_ce=1;																			
							}
							else 
							{
								if(ysum<=xsum)ori_ce=2;										
								else ori_ce=3;				 
							}
						}
						else
						{
							if(prf_p.cep[n][x]=='-')
							{									
								if(ysum<xsum)ori_ce=1;
								else ori_ce=0; 
							}
							else 
							{
								if(ysum<xsum)ori_ce=3;
								else ori_ce=2; 
							}
						}							
						cepi_sit_one[ori_ce][cepi_pos]++;
						if(r_dif==0)
						{
							if (cepi_pos<noveri)ce_full_eq=1;
							else
							{
								if (cepi_pos<nover)ce_part_eq=1;
								else ce_spac_eq++;
							}
						}
						else
						{
							if(r_dif<0)
							{
								if (cepi_pos<noveri)ce_full_anc=1;
								else
								{
									if (cepi_pos<nover)ce_part_anc=1;
									else ce_spac_anc++;
								}
							}
							else 
							{
								if (cepi_pos<noveri)ce_full_par=1;
								else
								{
									if (cepi_pos<nover)ce_part_par=1;
									else ce_spac_par++;
								}

							}
						}
						if(cepi_pos<noveri)ce_full[r_anc][r_par]++;
						else
						{
							if(cepi_pos<nover)ce_part[r_anc][r_par]++;	
							else ce_spac[r_anc][r_par]++;
						}																							
						int dir;
						if(prf_p.cep[n][x]=='+')dir=1;
						else dir=-1;																									
						if(real==1)
						{
							fprintf(outbest,"Seq %d\t",n+1);							
							fprintf(outbest,"%d\t%d\t",prf_a.sta[n][y],prf_a.sta[n][y]+len_a-1);
							fprintf(outbest,"%d\t%d\t",prf_p.sta[n][x],prf_p.sta[n][x]+len_p-1);							
							if(full_len!=-1)fprintf(outbest,"%d%c\t%s",full_len,locs[0][0],locs[0]);
							else
							{
								if(partlial_len!=-1)fprintf(outbest,"%d%c\t%s",partlial_len,locs[1][0],locs[1]);	
								else fprintf(outbest,"%d%c\t%s",spacer_len,locs[2][0],locs[2]);	
							}
							fprintf(outbest,"\t%c%c\t%s\t",prf_a.cep[n][y],prf_p.cep[n][x],oris[ori_ce]);																								
							fprintf(outbest,"%f\t%f\t",prf_a.pv[n][y],prf_p.pv[n][x]);							
						//	if (strstr(rera, "real") != NULL)
							{
								char dseq[2][MATLEN], compl1;											
								int pos;
								if(prf_a.cep[n][y]=='+')
								{
									pos=prf_a.sta[n][y];
									compl1=0;
								}
								else
								{
									pos=lenp-len_a-prf_a.sta[n][y];
									compl1=1;
								}
								strncpy(dseq[0],&seq[compl1][n][pos],len_a);	
								dseq[0][len_a]='\0';
								if(prf_p.cep[n][x]=='+')
								{
									pos=prf_p.sta[n][x];
									compl1=0;
								}
								else
								{
									pos=lenp-len_p-prf_p.sta[n][x];
									compl1=1;
								}
								strncpy(dseq[1],&seq[compl1][n][pos],len_p);	
								dseq[1][len_p]='\0';
								fprintf(outbest,"%s\t%s",dseq[0],dseq[1]);							
							}	
							fprintf(outbest, "\n");
							//printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%8d %8d",n_site,tot);													
						}
					}									
				}
			}
			if(join_sites_anc==1 || join_sites_par==1)
			{
				for(j=0;j<NUM_THR;j++)
				{
					for(k=0;k<NUM_THR;k++)
					{
						if(join_sites[j][k]==1)sam->cell[j][k].two_sites++;
					}
				}
			}
			if(join_sites_anc==1)
			{
				sam->anc.two_sites++;
				if(ce_full_anc+ce_part_anc+ce_spac_anc>0)
				{
					sam->anc.any++;
					if(ce_full_anc+ce_part_anc>0)sam->anc.overlap++;
					if(ce_spac_anc>0)sam->anc.spacer++;
					if(ce_part_anc>0)sam->anc.partial++;
					if(ce_full_anc>0)sam->anc.full++;
				}
			}
			if(join_sites_par==1)
			{
				sam->par.two_sites++;
			//if(join_sites_eq==1)sam->eq.two_sites++;			
				if(ce_full_par+ce_part_par+ce_spac_par>0)
				{
					sam->par.any++;
					if(ce_full_par+ce_part_par>0)sam->par.overlap++;
					if(ce_spac_par>0)sam->par.spacer++;
					if(ce_part_par>0)sam->par.partial++;
					if(ce_full_par>0)sam->par.full++;
				}
			}
			if (join_sites_anc == 1 || join_sites_par == 1)nseq_two_sites++;
			/*if(ce_full_eq+ce_part_eq+ce_spac_eq>0)
			{
				sam->eq.any++;
				if(ce_full_eq+ce_part_eq>0)sam->eq.overlap++;
				if(ce_spac_eq>0)sam->eq.spacer++;
				if(ce_part_eq>0)sam->eq.partial++;
				if(ce_full_eq>0)sam->eq.full++;
			}*/
			for(j=0;j<NUM_THR;j++)
			{
				for(k=0;k<NUM_THR;k++)
				{
					if(ce_full[j][k]+ce_part[j][k]+ce_spac[j][k]>0)
					{
						sam->cell[j][k].any++;
						if(ce_full[j][k]+ce_part[j][k]>0)sam->cell[j][k].overlap++;
						if(ce_spac[j][k]>0)sam->cell[j][k].spacer++;
						if(ce_part[j][k]>0)sam->cell[j][k].partial++;
						if(ce_full[j][k]>0)sam->cell[j][k].full++;
					}
				}
			}
			int kmin = -1;
			for(k=0;k<cepi_len;k++)
			{
				int take = 0;
				for(x=0;x<4;x++)
				{
					if(cepi_sit_one[x][k]>0)
					{
						take = 1;
						cepi_seq[x][k]++;
						//cepi_sit[x][k]+=cepi_sit_one[x][k];
					}
				}
				if(cepi_sit_one[0][k]+cepi_sit_one[1][k]>0)cepi_seq_dir[k]++;
				if (take == 1)
				{
					cepi_seqa[k]++;
					if(kmin==-1)kmin = k;
				}
			}
			if(kmin!=-1)for (k = kmin; k < cepi_len; k++)cepi_seqc[k]++;
		}
		else
		{				
			if(prf_p.nsit[n]>0)nseq_rec_err[0]++;
			if(prf_a.nsit[n]>0)nseq_rec_err[1]++;				
		}
		if(real==1)
		{
			fprintf(out_nsit_anc,"%d\n",asym_anc);
			fprintf(out_nsit_par,"%d\n",asym_par);
		}
	}	
	//strcpy(filerec,"projoin.txt");
//	printf("Debug - Start print!\n");
	if(real==1)
	{
		fclose(outbest);
		fclose(out_nsit_anc);
		fclose(out_nsit_par);
	}
	if(prf_a.mot==prf_p.mot)
	{
		for(j=0;j<cepi_len;j++)
		{			
			for(k=0;k<2;k++)
			{	
				hist->freq[k][j]=(double)cepi_seq_dir[j]/nseq;
			}
			for(k=2;k<4;k++)hist->freq[k][j]=(double)cepi_seq[k][j]/nseq;
		}
	}
	else
	{
		for(k=0;k<4;k++)
		{
			for(j=0;j<cepi_len;j++)
			{
				hist->freq[k][j]=(double)cepi_seq[k][j]/nseq;				
			}
		}		
	}
	for (j = 0; j < cepi_len; j++)hist->freqa[j] = (double)cepi_seqa[j] / nseq_both;
	for (j = 0; j < cepi_len; j++)hist->freqc[j] = (double)cepi_seqc[j] / nseq_both;
	//	printf("Debug - At the end! delete\n");
	for(j=0;j<4;j++)delete [] cepi_seq[j];
//	for(j=0;j<4;j++)delete [] cepi_sit[j];
	for(j=0;j<4;j++)delete [] cepi_sit_one[j];
	delete [] cepi_seq_dir;
	delete [] cepi_seqa;
	delete [] cepi_seqc;
	//printf("Debug - At the end! return\n");
	return 1;	
}
