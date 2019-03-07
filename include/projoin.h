int projoin_one(char *rera, char *motif_a, char *motif_p,profile prf_a, profile prf_p, int shift_min, int shift_max, int len_a, int len_p, int *thr_pre_err,
			int nseq, char ***seq, count *sample, combi *hist, int *peak_len)
{	
	int n, k,j, x,y;
	char filebest[120], filerec[80]; 

	FILE *out2, *outbest;// *outhead,*out,  
//	FILE *out_cepi_seq;//*out_cepi_sit, *out_over_spac;//*outnsite_any,*out_hist,

	memset(filebest,'\0',sizeof(filebest));
	memset(filerec,'\0',sizeof(filerec));
	strcpy(filebest,rera);//real or random
	strcat(filebest,"_");
	strcat(filebest,motif_a);//hocomoco or dapseq
	strcat(filebest,"_");
	strcat(filebest,motif_p);//hocomoco or dapseq
	char buf[10];
	sprintf(buf,"%d",prf_a.mot);
	strcat(filebest,buf);
	sprintf(buf,"%d",prf_p.mot);
	strcat(filebest,buf);
	strcpy(filerec,filebest);
	strcat(filerec,".txt");
	strcat(filebest,"_thr");
	sprintf(buf,"%d",prf_a.nam);
	strcat(filebest,buf);
	sprintf(buf,"%d",prf_p.nam);
	strcat(filebest,buf);
	strcat(filebest,".best");
		
	int *cepi_seq[4], *cepi_seq_dir;
	//int *cepi_sit[4];
	int *cepi_sit_one[4];
	//if(strncmp(rera,"rand",4)!=0)
	if(strstr(rera,"real_one")!=NULL)
	{
		if((outbest=fopen(filebest,"wt"))==NULL)
		{
 			printf("Input file %s can't be opened!\n",filebest);
			return -1;
		}
		fprintf(outbest,"#Seq\tA Start\tA End\tP Start\tP End\tMutual Loc\tLoc Type\tStrands\tMutual Ori\tA Score\tP Score\tA Seq\tP Seq\n");
	}
	int noverp=Min(len_p,len_a);//partial overlap
	noverp--;
	int noveri=1+abs(len_p-len_a)/2;// full overlap	
	int cepi_len=shift_max+noveri+noverp+1;	
	//hist->n_tot=cepi_len;
	int nover=noveri+noverp;
	int n_peak_over=0, n_peak_spacer=0, n_peak_full_over=0, n_peak_part_over=0;// no. of peak CE overlap, spacer only CE, full overlap, partial overlap
	int n_ce_part_over=0, n_ce_full_over=0, n_ce_spacer=0;// no. of CE overlap, spacer CE
	for(j=0;j<4;j++)
	{
		cepi_seq[j]=new int[cepi_len];
		if(cepi_seq[j]==NULL){printf("Not  enough memory!");return -1;}	
		//cepi_sit[j]=new int[cepi_len];
		//if(cepi_sit[j]==NULL){printf("Not  enough memory!");return -1;}	
		cepi_sit_one[j]=new int[cepi_len];
		if(cepi_sit_one[j]==NULL){printf("Not  enough memory!");return -1;}	
		for(k=0;k<cepi_len;k++)cepi_seq[j][k]=0;
		//for(k=0;k<cepi_len;k++)cepi_sit[j][k]=0;		
	}
	cepi_seq_dir=new int[cepi_len];
	if(cepi_seq_dir==NULL){printf("Not  enough memory!");return -1;}	
	for(k=0;k<cepi_len;k++)cepi_seq_dir[k]=0;		
	
	int tot=0;	
	int rec_pos_seq=0;
	int rec_seq=0;
	int nseq_rec[2]={0,0};//sequences recognized by 1st(2nd) profile but not recognized by 2nd(1st) profile
	int nseq_join_sites=0;// both sites are present
	int nseq_rec_tot[2]={0,0};//sequences recognized by 1st(2nd) profile	
	int nseq_rec_err[2]={0,0};//sequences that should be ignored in comparison of two profiles
	char oris[4][10]={"DirectAP","DirectPA","Inverted","Everted"};//direct_anchor_partner, direct partner_anchor, invert_anchor_partner, evert_anchor_partner
	char locs[3][10]={"Full","Partial","Spacer"};

	for(n=0;n<nseq;n++)
	{						
		if(prf_a.nsit[n]==0 || prf_p.nsit[n]==0)continue;
		int lenp=peak_len[n];
		int n_ce_full_over_here=0, n_ce_part_over_here=0, n_ce_spacer_here=0;// no. of CE overlap per peak, no. of spacer CE per peak
		if(thr_pre_err[n]==0)
		{								
//			printf("Nseq %d Anc %d Par %d\n",n+1,prf_a.nsit[n],prf_p.nsit[n]);
	//		for(y=0;y<prf_a.nsit[n];y++)printf("%d%c\t",prf_a.sta[n][y],prf_a.cep[n][y]);
		//	printf("\n");
		//	for(x=0;x<prf_p.nsit[n];x++)printf("%d%c\t",prf_p.sta[n][x],prf_p.cep[n][x]);
		//	printf("\n");
			for(k=0;k<cepi_len;k++)for(j=0;j<4;j++)cepi_sit_one[j][k]=0;														
			//double scobest[2]={-1,-1};
			//double max_score1=0, max_score2=0;
			if(prf_a.nsit[n]>0 && prf_p.nsit[n]>0)nseq_join_sites++; 
			int dlen=0;
			{
				int dif_len=len_p-len_a;
				if(dif_len>=0)dlen=dif_len-(dif_len)/2;
				else dlen=dif_len/2;
			}
			rec_pos_seq=0;
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
				for(y=0;y<prf_a.nsit[n];y++)
				{
					if(prf_a.mot==prf_p.mot)
					{
						if(x<=y)continue;//homodimer
					}
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
					if(cat[0]==1 && cat[1]==1)
					{
						take_distance=1;//full overlap
						full_len=Min(ysta-xsta,xend-yend);
					}
					else
					{
						if(cat[2]==1 && cat[3]==1)
						{
							take_distance=1;//full overlap
							full_len=Min(xsta-ysta,yend-xend);
						}
					}
					if(take_distance==0)
					{
						if(cat[0]==1 && cat[3]==1)
						{
							take_distance=1;//partial overlap
							partlial_len=xend-ysta+1;
						}
						else 
						{
							if(cat[1]==1 && cat[2]==1)
							{
								take_distance=1;//partial overlap
								partlial_len=yend-xsta+1;
							}
						}
					}
					if(take_distance==0)//check for spacer
					{
						if(xend<ysta)spacer_len=ysta-xend-1;
						if(yend<xsta)spacer_len=xsta-yend-1;
						if(spacer_len>=shift_min && spacer_len<=shift_max)take_distance=1;//spacer
					}					
					if(take_distance==1)
					{															
//							printf("2x=%d\t%d\t%d\t\t",x,xsta,xend);
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
						if(full_len!=-1)n_ce_full_over_here++;
						else
						{
							if(partlial_len!=-1)n_ce_part_over_here++;	
							else n_ce_spacer_here++;
						}																							
						int dir;
						if(prf_p.cep[n][x]=='+')dir=1;
						else dir=-1;														
						if(rec_pos_seq==0)rec_seq++;								
						tot++;							
						rec_pos_seq++;							
	//					if(strncmp(rera,"rand",4)!=0)
						if(strstr(rera,"real_one")!=NULL)
						{
							fprintf(outbest,"Seq %d\t",n+1);	
						//	if(pre[1][y].cep=='-'){Mix(&x1,&y1);Mix(&zero1,&one1);}
							fprintf(outbest,"%d\t%d\t",prf_a.sta[n][y],prf_a.sta[n][y]+len_a-1);
							fprintf(outbest,"%d\t%d\t",prf_p.sta[n][x],prf_p.sta[n][x]+len_p-1);
							if(full_len!=-1)fprintf(outbest,"%d%c\t%s",full_len,locs[0][0],locs[0]);
							else
							{
								if(partlial_len!=-1)fprintf(outbest,"%d%c\t%s",partlial_len,locs[1][0],locs[1]);	
								else fprintf(outbest,"%d%c\t%s",spacer_len,locs[2][0],locs[2]);	
							}
							fprintf(outbest,"\t%c%c\t%s\t",prf_a.cep[n][y],prf_p.cep[n][x],oris[ori_ce]);												
							fprintf(outbest,"%f\t%f\t",prf_a.sco[n][y],prf_p.sco[n][x]);							
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
								fprintf(outbest,"%s\t%s\n",dseq[0],dseq[1]);							
							}																						
							//printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%8d %8d",n_site,tot);													
						}
					}									
				}
			}
			if(n_ce_full_over_here+n_ce_part_over_here>0)
			{
				n_peak_over++;					
				if(n_ce_full_over_here>0)
				{
					n_peak_full_over++;
					n_ce_full_over+=n_ce_full_over_here;
				}
				if(n_ce_part_over_here>0)
				{
					n_peak_part_over++;
					n_ce_part_over+=n_ce_part_over_here;
				} 
			}
			if(n_ce_spacer_here>0)
			{
				n_peak_spacer++;				
				n_ce_spacer+=n_ce_spacer_here;
			}																
			for(k=0;k<cepi_len;k++)
			{
				for(x=0;x<4;x++)
				{
					if(cepi_sit_one[x][k]>0)
					{
						cepi_seq[x][k]++;
						//cepi_sit[x][k]+=cepi_sit_one[x][k];
					}
				}
				if(cepi_sit_one[0][k]+cepi_sit_one[1][k]>0)cepi_seq_dir[k]++;
			}							
		}
		else
		{				
			if(prf_p.nsit[n]>0)nseq_rec_err[0]++;
			if(prf_a.nsit[n]>0)nseq_rec_err[1]++;				
		}
	}	
	//strcpy(filerec,"projoin.txt");
//	printf("Debug - Start print!\n");
	if(prf_a.nam==1 && prf_p.nam==1)
	{
		if((out2=fopen(filerec,"wt"))==NULL)
		{
 			printf("Input file %s can't be opened!\n",filerec);
			return -1;
		}
	}
	else
	{
		if((out2=fopen(filerec,"at"))==NULL)
		{
 			printf("Input file %s can't be opened!\n",filerec);
			return -1;
		}
	}
	//                                      join_file seq(join_sites) 4islo par,4islo par (na kajdy sayt po 1 best pare  n_site[0,1]=1st,2nd profile_total_sites  sequences with both sites(no overlap)   
	fprintf(out2,"%s Anc_%d(%d) Par %d(%d)\tSpacer %d..%d\t",rera,prf_a.mot,prf_a.nam,prf_p.mot,prf_p.nam,shift_min,shift_max);
	fprintf(out2,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",rec_seq,tot,tot,prf_a.nsit_all,prf_p.nsit_all,nseq_rec[0], nseq_rec[1],nseq_join_sites,prf_a.nseq_rec-nseq_rec_err[0], prf_p.nseq_rec-nseq_rec_err[1]);
	fprintf(out2,"%d\t%d\t%d\t%d\n",n_peak_full_over,n_peak_part_over,n_peak_over,n_peak_spacer);
	fclose(out2);
//  fclose(outhead);
	//if(strncmp(rera,"rand",4)!=0)
	if(strstr(rera,"real_one")!=NULL)fclose(outbest);
	//if(strcmp(motif_a,motif_p)==0)
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
				
	//printf("Debug - End print1!\n");

//	printf("Debug - At the end! cepi\n");
/*	if((out_cepi_seq=fopen(file_cepi_seq,"at"))==NULL)
	{
 		printf("Input file %s can't be opened!\n",file_cepi_seq);
		return -1;
	}
	if((out_cepi_sit=fopen(file_cepi_sit,"at"))==NULL)
	{
 		printf("Input file %s can't be opened!\n",file_cepi_sit);
		return -1;
	}
	if((out_over_spac=fopen(file_over_spac,"at"))==NULL)
	{
 		printf("Input file %s can't be opened!\n",file_over_spac);
		return -1;
	}*/	
	/*
	int print_overlap=1;	//0 ne pe4araem 2 pe4ataem polovinu 1 pe4ataem vse
	{
		char spacer[]="S";
		char border[]="P";//partial
		char internal_loc[]="F";	//full			
		{
			fprintf(out_cepi_seq,"\t");
		//	fprintf(out_cepi_sit,"\t");
			if(print_overlap==1)
			{
				if(len_p==len_a)
				{
					fprintf(out_cepi_seq,"\t0%s",internal_loc);
			//		fprintf(out_cepi_sit,"\t0%s",internal_loc);
				}
				else
				{				
					for(j=noveri-1;j>=0;j--)
					{					
				//		fprintf(out_cepi_sit,"\t%d%s",j,internal_loc);
						fprintf(out_cepi_seq,"\t%d%s",j,internal_loc);
					}
				}		
			}
		}
		if(print_overlap==1)
		{
			//for(j=noverp;j>=1;j--)fprintf(out_cepi_sit,"\t%d%s",j,border);
			for(j=noverp;j>=1;j--)fprintf(out_cepi_seq,"\t%d%s",j,border);
		}
		if(print_overlap==2)
		{
			int noverp2=noverp/2;
			//for(j=noverp2;j>=1;j--)fprintf(out_cepi_sit,"\t%d%s",j,border);
			for(j=noverp2;j>=1;j--)fprintf(out_cepi_seq,"\t%d%s",j,border);
		}
		//for(j=0;j<shift;j++)fprintf(out_cepi_sit,"\t%d%s",j,spacer);
		for(j=0;j<shift;j++)fprintf(out_cepi_seq,"\t%d%s",j,spacer);
	}
	fprintf(out_cepi_seq,"\n");
	//fprintf(out_cepi_sit,"\n");
	if(prf_a.mot==prf_p.mot)
	{		
//		fprintf(out_over_spac,"%s Anc_%d Par_%d\t\t%d\t%d\t%d\t\t%d\t%d\t%d\t%d\t\t",rera,prf_a.mot,prf_p.mot,nseq,nseq_join_sites,rec_seq,n_peak_full_over,n_peak_part_over,n_peak_over,n_peak_spacer);
	//	fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/nseq,(double)n_peak_full_over/nseq,(double)n_peak_part_over/nseq,(double)n_peak_spacer/nseq);
//		fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/nseq_join_sites,(double)n_peak_full_over/nseq_join_sites,(double)n_peak_part_over/nseq_join_sites,(double)n_peak_spacer/nseq_join_sites);
	//	fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/rec_seq,(double)n_peak_full_over/rec_seq,(double)n_peak_part_over/rec_seq,(double)n_peak_spacer/rec_seq);		
		//int n_ce_any=n_ce_full_over+n_ce_part_over+n_ce_spacer;
//		fprintf(out_over_spac,"%d\t\t",n_ce_any);		
	//	fprintf(out_over_spac,"%d\t%d\t%d\t\t",n_ce_full_over,n_ce_part_over,n_ce_spacer);		
//		fprintf(out_over_spac,"%f\t%f\t%f\n",(double)n_ce_full_over/n_ce_any,(double)n_ce_part_over/n_ce_any,(double)n_ce_spacer/n_ce_any);		
		for(k=3;k>=0;k--)
		{
			switch (k)
			{
			case 0: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tDirectAP");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tDirectAP");
				break;}
			case 1: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tDirectPA");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tDirectPA");
				break;}
			case 2: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tInvert");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tInvert");
				break;}
			default: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tEvert");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tEvert");
				break;}			
			}			
			int umk=1-k;	
			int j0;
			if(print_overlap==1)j0=0;
			if(print_overlap==0)j0=noveri+noverp;
			if(print_overlap==2)j0=noveri+noverp-noverp/2;
			if(k>1)
			{
				for(j=j0;j<cepi_len;j++)
				{
					double freq=(double)cepi_seq[k][j]/nseq;
					hist->freq[k][j-j0]=freq;
					fprintf(out_cepi_seq,"\t%f",freq);
					//fprintf(out_cepi_sit,"\t%f",(double)cepi_sit[k][j]/nseq);
				}
			}
			else 
			{
				for(j=j0;j<cepi_len;j++)
				{								
					double freq=(double)cepi_seq_dir[j]/nseq;
					hist->freq[k][j-j0]=freq;
					fprintf(out_cepi_seq,"\t%f",freq);					
					//fprintf(out_cepi_sit,"\t%f",(double)(cepi_sit[k][j]+cepi_sit[umk][j])/nseq);
				}			
			}
			fprintf(out_cepi_seq,"\n");
			//fprintf(out_cepi_sit,"\n");
		}
	}
	else
	{			
		//fprintf(out_over_spac,"%s Anc_%d Par_%d\t\t%d\t%d\t%d\t\t%d\t%d\t%d\t%d\t\t",rera, prf_a.mot,prf_p.mot,nseq,nseq_join_sites,rec_seq,n_peak_full_over,n_peak_part_over,n_peak_over,n_peak_spacer);
//		fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/nseq,(double)n_peak_full_over/nseq,(double)n_peak_part_over/nseq,(double)n_peak_spacer/nseq);
	//	fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/nseq_join_sites,(double)n_peak_full_over/nseq_join_sites,(double)n_peak_part_over/nseq_join_sites,(double)n_peak_spacer/nseq_join_sites);
//		fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/rec_seq,(double)n_peak_full_over/rec_seq,(double)n_peak_part_over/rec_seq,(double)n_peak_spacer/rec_seq);		
		//int n_ce_any=n_ce_full_over+n_ce_part_over+n_ce_spacer;
	//	fprintf(out_over_spac,"%d\t\t",n_ce_any);		
		//fprintf(out_over_spac,"%d\t%d\t%d\t\t",n_ce_full_over,n_ce_part_over,n_ce_spacer);		
		//fprintf(out_over_spac,"%f\t%f\t%f\n",(double)n_ce_full_over/n_ce_any,(double)n_ce_part_over/n_ce_any,(double)n_ce_spacer/n_ce_any);		
		for(k=3;k>=0;k--)
		{			
			switch (k)
			{
			case 0: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tDirectAP");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tDirectAP");
				break;}
			case 1: {
				fprintf(out_cepi_seq,"Par_%d Anc_%d",prf_p.mot,prf_a.mot);fprintf(out_cepi_seq,"\tDirectPA");
				//fprintf(out_cepi_sit,"Par_%d Anc_%d",prf_p.mot,prf_a.mot);fprintf(out_cepi_sit,"\tDirectPA");
				break;}
			case 2: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tInvert");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tInvert");
				break;}
			default: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tEvert");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tEvert");
				break;}			
			}			
			int j0;
			if(print_overlap==0)j0=noveri+noverp;
			if(print_overlap==1)j0=0;
			if(print_overlap==2)j0=noveri+noverp-noverp/2;
			for(j=j0;j<cepi_len;j++)
			{
				double freq=(double)cepi_seq[k][j]/nseq;
				hist->freq[k][j-j0]=freq;
				fprintf(out_cepi_seq,"\t%f",freq);
				//fprintf(out_cepi_sit,"\t%f",(double)cepi_sit[k][j]/nseq);
			}
			fprintf(out_cepi_seq,"\n");
			//fprintf(out_cepi_sit,"\n");
		}
	}
	fclose(out_cepi_seq);
	*/
	//fclose(out_cepi_sit);
	//fclose(out_over_spac);
	//	printf("Debug - At the end! delete\n");
	for(j=0;j<4;j++)delete [] cepi_seq[j];
//	for(j=0;j<4;j++)delete [] cepi_sit[j];
	for(j=0;j<4;j++)delete [] cepi_sit_one[j];
	delete [] cepi_seq_dir;
	//printf("Debug - At the end! return\n");
	sample->two_sites=nseq_join_sites;
	sample->any=rec_seq;
	sample->full=n_peak_full_over;
	sample->partial=n_peak_part_over;
	sample->over=n_peak_over;
	sample->spacer=n_peak_spacer;
	return 1;	
}
int projoin(char *rera, char *motif,profile prf_a, profile prf_p, int shift_min, int shift_max, int len_a, int len_p, int *thr_pre_err,
			int nseq, char ***seq, count *sample, combi *hist, int *peak_len)
{	
	int n, k,j, x,y;
	char filebest[120], filerec[30]; 

	FILE *out2, *outbest;// *outhead,*out,  
//	FILE *out_cepi_seq;//*out_cepi_sit, *out_over_spac;//*outnsite_any,*out_hist,

	memset(filebest,'\0',sizeof(filebest));
	memset(filerec,'\0',sizeof(filerec));
	strcpy(filebest,rera);//real or random
	strcat(filebest,"_");
	strcat(filebest,motif);//hocomoco or dapseq
	char buf[10];
	sprintf(buf,"%d",prf_a.mot);
	strcat(filebest,buf);
	sprintf(buf,"%d",prf_p.mot);
	strcat(filebest,buf);
	strcpy(filerec,filebest);
	strcat(filerec,".txt");
	strcat(filebest,"_thr");
	sprintf(buf,"%d",prf_a.nam);
	strcat(filebest,buf);
	sprintf(buf,"%d",prf_p.nam);
	strcat(filebest,buf);
	strcat(filebest,".best");
		
	int *cepi_seq[4], *cepi_seq_dir;
	//int *cepi_sit[4];
	int *cepi_sit_one[4];
	//if(strncmp(rera,"rand",4)!=0)
	if(strstr(rera,"real_one")!=NULL)
	{
		if((outbest=fopen(filebest,"wt"))==NULL)
		{
 			printf("Input file %s can't be opened!\n",filebest);
			return -1;
		}
		fprintf(outbest,"#Seq\tA Start\tA End\tP Start\tP End\tMutual Loc\tLoc Type\tStrands\tMutual Ori\tA Score\tP Score\tA Seq\tP Seq\n");
	}
	int noverp=Min(len_p,len_a);//partial overlap
	noverp--;
	int noveri=1+abs(len_p-len_a)/2;// full overlap	
	int cepi_len=shift_max+noveri+noverp+1;	
	//hist->n_tot=cepi_len;
	int nover=noveri+noverp;
	int n_peak_over=0, n_peak_spacer=0, n_peak_full_over=0, n_peak_part_over=0;// no. of peak CE overlap, spacer only CE, full overlap, partial overlap
	int n_ce_part_over=0, n_ce_full_over=0, n_ce_spacer=0;// no. of CE overlap, spacer CE
	for(j=0;j<4;j++)
	{
		cepi_seq[j]=new int[cepi_len];
		if(cepi_seq[j]==NULL){printf("Not  enough memory!");return -1;}	
		//cepi_sit[j]=new int[cepi_len];
		//if(cepi_sit[j]==NULL){printf("Not  enough memory!");return -1;}	
		cepi_sit_one[j]=new int[cepi_len];
		if(cepi_sit_one[j]==NULL){printf("Not  enough memory!");return -1;}	
		for(k=0;k<cepi_len;k++)cepi_seq[j][k]=0;
		//for(k=0;k<cepi_len;k++)cepi_sit[j][k]=0;		
	}
	cepi_seq_dir=new int[cepi_len];
	if(cepi_seq_dir==NULL){printf("Not  enough memory!");return -1;}	
	for(k=0;k<cepi_len;k++)cepi_seq_dir[k]=0;		
	
	int tot=0;	
	int rec_pos_seq=0;
	int rec_seq=0;
	int nseq_rec[2]={0,0};//sequences recognized by 1st(2nd) profile but not recognized by 2nd(1st) profile
	int nseq_join_sites=0;// both sites are present
	int nseq_rec_tot[2]={0,0};//sequences recognized by 1st(2nd) profile	
	int nseq_rec_err[2]={0,0};//sequences that should be ignored in comparison of two profiles
	char oris[4][10]={"DirectAP","DirectPA","Inverted","Everted"};//direct_anchor_partner, direct partner_anchor, invert_anchor_partner, evert_anchor_partner
	char locs[3][10]={"Full","Partial","Spacer"};

	for(n=0;n<nseq;n++)
	{						
		if(prf_a.nsit[n]==0 || prf_p.nsit[n]==0)continue;
		int lenp=peak_len[n];
		int n_ce_full_over_here=0, n_ce_part_over_here=0, n_ce_spacer_here=0;// no. of CE overlap per peak, no. of spacer CE per peak
		if(thr_pre_err[n]==0)
		{								
//			printf("Nseq %d Anc %d Par %d\n",n+1,prf_a.nsit[n],prf_p.nsit[n]);
	//		for(y=0;y<prf_a.nsit[n];y++)printf("%d%c\t",prf_a.sta[n][y],prf_a.cep[n][y]);
		//	printf("\n");
		//	for(x=0;x<prf_p.nsit[n];x++)printf("%d%c\t",prf_p.sta[n][x],prf_p.cep[n][x]);
		//	printf("\n");
			for(k=0;k<cepi_len;k++)for(j=0;j<4;j++)cepi_sit_one[j][k]=0;														
			//double scobest[2]={-1,-1};
			//double max_score1=0, max_score2=0;
			if(prf_a.nsit[n]>0 && prf_p.nsit[n]>0)nseq_join_sites++; 
			int dlen=0;
			{
				int dif_len=len_p-len_a;
				if(dif_len>=0)dlen=dif_len-(dif_len)/2;
				else dlen=dif_len/2;
			}
			rec_pos_seq=0;
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
				for(y=0;y<prf_a.nsit[n];y++)
				{
					if(prf_a.mot==prf_p.mot)
					{
						if(x<=y)continue;//homodimer
					}
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
					}
					else
					{
						if(cat[2]==1 && cat[3]==1)
						{
							take_distance=1;//full overlap
							full_len=Min(xsta-ysta,yend-xend);
						}
					}
					if(take_distance==0)
					{
						if(cat[0]==1 && cat[3]==1)
						{
							take_distance=1;//partial overlap
							partlial_len=xend-ysta+1;
						}
						else 
						{
							if(cat[1]==1 && cat[2]==1)
							{
								take_distance=1;//partial overlap
								partlial_len=yend-xsta+1;
							}
						}
					}
					if(take_distance==0)//check for spacer
					{
						if(xend<ysta)spacer_len=ysta-xend-1;
						if(yend<xsta)spacer_len=xsta-yend-1;
						if(spacer_len>=shift_min && spacer_len<=shift_max)take_distance=1;//spacer
					}					
					if(take_distance==1)
					{															
//							printf("2x=%d\t%d\t%d\t\t",x,xsta,xend);
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
						if(cepi_pos<noveri)n_ce_full_over_here++;
						else
						{
							if(cepi_pos<nover)n_ce_part_over_here++;	
							else n_ce_spacer_here++;
						}																							
						int dir;
						if(prf_p.cep[n][x]=='+')dir=1;
						else dir=-1;														
						if(rec_pos_seq==0)rec_seq++;								
						tot++;							
						rec_pos_seq++;							
	//					if(strncmp(rera,"rand",4)!=0)
						if(strstr(rera,"real_one")!=NULL)
						{
							fprintf(outbest,"Seq %d\t",n+1);	
						//	if(pre[1][y].cep=='-'){Mix(&x1,&y1);Mix(&zero1,&one1);}
							fprintf(outbest,"%d\t%d\t",prf_a.sta[n][y],prf_a.sta[n][y]+len_a-1);
							fprintf(outbest,"%d\t%d\t",prf_p.sta[n][x],prf_p.sta[n][x]+len_p-1);							
							if(full_len!=-1)fprintf(outbest,"%d%c\t%s",full_len,locs[0][0],locs[0]);
							else
							{
								if(partlial_len!=-1)fprintf(outbest,"%d%c\t%s",partlial_len,locs[1][0],locs[1]);	
								else fprintf(outbest,"%d%c\t%s",spacer_len,locs[2][0],locs[2]);	
							}
							fprintf(outbest,"\t%c%c\t%s\t",prf_a.cep[n][y],prf_p.cep[n][x],oris[ori_ce]);																								
							fprintf(outbest,"%f\t%f\t",prf_a.sco[n][y],prf_p.sco[n][x]);							
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
								fprintf(outbest,"%s\t%s\n",dseq[0],dseq[1]);							
							}																						
							//printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%8d %8d",n_site,tot);													
						}
					}									
				}
			}
			if(n_ce_full_over_here+n_ce_part_over_here>0)
			{
				n_peak_over++;					
				if(n_ce_full_over_here>0)
				{
					n_peak_full_over++;
					n_ce_full_over+=n_ce_full_over_here;
				}
				if(n_ce_part_over_here>0)
				{
					n_peak_part_over++;
					n_ce_part_over+=n_ce_part_over_here;
				} 
			}
			if(n_ce_spacer_here>0)
			{
				n_peak_spacer++;				
				n_ce_spacer+=n_ce_spacer_here;
			}																
			for(k=0;k<cepi_len;k++)
			{
				for(x=0;x<4;x++)
				{
					if(cepi_sit_one[x][k]>0)
					{
						cepi_seq[x][k]++;
						//cepi_sit[x][k]+=cepi_sit_one[x][k];
					}
				}
				if(cepi_sit_one[0][k]+cepi_sit_one[1][k]>0)cepi_seq_dir[k]++;
			}							
		}
		else
		{				
			if(prf_p.nsit[n]>0)nseq_rec_err[0]++;
			if(prf_a.nsit[n]>0)nseq_rec_err[1]++;				
		}
	}	
	//strcpy(filerec,"projoin.txt");
//	printf("Debug - Start print!\n");
	if(prf_a.nam==1 && prf_p.nam==1)
	{
		if((out2=fopen(filerec,"wt"))==NULL)
		{
 			printf("Input file %s can't be opened!\n",filerec);
			return -1;
		}
	}
	else
	{
		if((out2=fopen(filerec,"at"))==NULL)
		{
 			printf("Input file %s can't be opened!\n",filerec);
			return -1;
		}
	}
	//                                      join_file seq(join_sites) 4islo par,4islo par (na kajdy sayt po 1 best pare  n_site[0,1]=1st,2nd profile_total_sites  sequences with both sites(no overlap)   
	fprintf(out2,"%s Anc_%d(%d) Par %d(%d)\tSpacer %d..%d\t",rera,prf_a.mot,prf_a.nam,prf_p.mot,prf_p.nam,shift_min,shift_max);
	fprintf(out2,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",rec_seq,tot,tot,prf_a.nsit_all,prf_p.nsit_all,nseq_rec[0], nseq_rec[1],nseq_join_sites,prf_a.nseq_rec-nseq_rec_err[0], prf_p.nseq_rec-nseq_rec_err[1]);
	fprintf(out2,"%d\t%d\t%d\t%d\n",n_peak_full_over,n_peak_part_over,n_peak_over,n_peak_spacer);
	fclose(out2);
//  fclose(outhead);
	//if(strncmp(rera,"rand",4)!=0)
	if(strstr(rera,"real_one")!=NULL)fclose(outbest);
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
				
	//printf("Debug - End print1!\n");

//	printf("Debug - At the end! cepi\n");
/*	if((out_cepi_seq=fopen(file_cepi_seq,"at"))==NULL)
	{
 		printf("Input file %s can't be opened!\n",file_cepi_seq);
		return -1;
	}
	if((out_cepi_sit=fopen(file_cepi_sit,"at"))==NULL)
	{
 		printf("Input file %s can't be opened!\n",file_cepi_sit);
		return -1;
	}
	if((out_over_spac=fopen(file_over_spac,"at"))==NULL)
	{
 		printf("Input file %s can't be opened!\n",file_over_spac);
		return -1;
	}*/	
	/*
	int print_overlap=1;	//0 ne pe4araem 2 pe4ataem polovinu 1 pe4ataem vse
	{
		char spacer[]="S";
		char border[]="P";//partial
		char internal_loc[]="F";	//full			
		{
			fprintf(out_cepi_seq,"\t");
		//	fprintf(out_cepi_sit,"\t");
			if(print_overlap==1)
			{
				if(len_p==len_a)
				{
					fprintf(out_cepi_seq,"\t0%s",internal_loc);
			//		fprintf(out_cepi_sit,"\t0%s",internal_loc);
				}
				else
				{				
					for(j=noveri-1;j>=0;j--)
					{					
				//		fprintf(out_cepi_sit,"\t%d%s",j,internal_loc);
						fprintf(out_cepi_seq,"\t%d%s",j,internal_loc);
					}
				}		
			}
		}
		if(print_overlap==1)
		{
			//for(j=noverp;j>=1;j--)fprintf(out_cepi_sit,"\t%d%s",j,border);
			for(j=noverp;j>=1;j--)fprintf(out_cepi_seq,"\t%d%s",j,border);
		}
		if(print_overlap==2)
		{
			int noverp2=noverp/2;
			//for(j=noverp2;j>=1;j--)fprintf(out_cepi_sit,"\t%d%s",j,border);
			for(j=noverp2;j>=1;j--)fprintf(out_cepi_seq,"\t%d%s",j,border);
		}
		//for(j=0;j<shift;j++)fprintf(out_cepi_sit,"\t%d%s",j,spacer);
		for(j=0;j<shift;j++)fprintf(out_cepi_seq,"\t%d%s",j,spacer);
	}
	fprintf(out_cepi_seq,"\n");
	//fprintf(out_cepi_sit,"\n");
	if(prf_a.mot==prf_p.mot)
	{		
//		fprintf(out_over_spac,"%s Anc_%d Par_%d\t\t%d\t%d\t%d\t\t%d\t%d\t%d\t%d\t\t",rera,prf_a.mot,prf_p.mot,nseq,nseq_join_sites,rec_seq,n_peak_full_over,n_peak_part_over,n_peak_over,n_peak_spacer);
	//	fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/nseq,(double)n_peak_full_over/nseq,(double)n_peak_part_over/nseq,(double)n_peak_spacer/nseq);
//		fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/nseq_join_sites,(double)n_peak_full_over/nseq_join_sites,(double)n_peak_part_over/nseq_join_sites,(double)n_peak_spacer/nseq_join_sites);
	//	fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/rec_seq,(double)n_peak_full_over/rec_seq,(double)n_peak_part_over/rec_seq,(double)n_peak_spacer/rec_seq);		
		//int n_ce_any=n_ce_full_over+n_ce_part_over+n_ce_spacer;
//		fprintf(out_over_spac,"%d\t\t",n_ce_any);		
	//	fprintf(out_over_spac,"%d\t%d\t%d\t\t",n_ce_full_over,n_ce_part_over,n_ce_spacer);		
//		fprintf(out_over_spac,"%f\t%f\t%f\n",(double)n_ce_full_over/n_ce_any,(double)n_ce_part_over/n_ce_any,(double)n_ce_spacer/n_ce_any);		
		for(k=3;k>=0;k--)
		{
			switch (k)
			{
			case 0: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tDirectAP");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tDirectAP");
				break;}
			case 1: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tDirectPA");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tDirectPA");
				break;}
			case 2: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tInvert");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tInvert");
				break;}
			default: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tEvert");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tEvert");
				break;}			
			}			
			int umk=1-k;	
			int j0;
			if(print_overlap==1)j0=0;
			if(print_overlap==0)j0=noveri+noverp;
			if(print_overlap==2)j0=noveri+noverp-noverp/2;
			if(k>1)
			{
				for(j=j0;j<cepi_len;j++)
				{
					double freq=(double)cepi_seq[k][j]/nseq;
					hist->freq[k][j-j0]=freq;
					fprintf(out_cepi_seq,"\t%f",freq);
					//fprintf(out_cepi_sit,"\t%f",(double)cepi_sit[k][j]/nseq);
				}
			}
			else 
			{
				for(j=j0;j<cepi_len;j++)
				{								
					double freq=(double)cepi_seq_dir[j]/nseq;
					hist->freq[k][j-j0]=freq;
					fprintf(out_cepi_seq,"\t%f",freq);					
					//fprintf(out_cepi_sit,"\t%f",(double)(cepi_sit[k][j]+cepi_sit[umk][j])/nseq);
				}			
			}
			fprintf(out_cepi_seq,"\n");
			//fprintf(out_cepi_sit,"\n");
		}
	}
	else
	{			
		//fprintf(out_over_spac,"%s Anc_%d Par_%d\t\t%d\t%d\t%d\t\t%d\t%d\t%d\t%d\t\t",rera, prf_a.mot,prf_p.mot,nseq,nseq_join_sites,rec_seq,n_peak_full_over,n_peak_part_over,n_peak_over,n_peak_spacer);
//		fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/nseq,(double)n_peak_full_over/nseq,(double)n_peak_part_over/nseq,(double)n_peak_spacer/nseq);
	//	fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/nseq_join_sites,(double)n_peak_full_over/nseq_join_sites,(double)n_peak_part_over/nseq_join_sites,(double)n_peak_spacer/nseq_join_sites);
//		fprintf(out_over_spac,"%f\t%f\t%f\t%f\t\t",(double)n_peak_over/rec_seq,(double)n_peak_full_over/rec_seq,(double)n_peak_part_over/rec_seq,(double)n_peak_spacer/rec_seq);		
		//int n_ce_any=n_ce_full_over+n_ce_part_over+n_ce_spacer;
	//	fprintf(out_over_spac,"%d\t\t",n_ce_any);		
		//fprintf(out_over_spac,"%d\t%d\t%d\t\t",n_ce_full_over,n_ce_part_over,n_ce_spacer);		
		//fprintf(out_over_spac,"%f\t%f\t%f\n",(double)n_ce_full_over/n_ce_any,(double)n_ce_part_over/n_ce_any,(double)n_ce_spacer/n_ce_any);		
		for(k=3;k>=0;k--)
		{			
			switch (k)
			{
			case 0: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tDirectAP");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tDirectAP");
				break;}
			case 1: {
				fprintf(out_cepi_seq,"Par_%d Anc_%d",prf_p.mot,prf_a.mot);fprintf(out_cepi_seq,"\tDirectPA");
				//fprintf(out_cepi_sit,"Par_%d Anc_%d",prf_p.mot,prf_a.mot);fprintf(out_cepi_sit,"\tDirectPA");
				break;}
			case 2: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tInvert");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tInvert");
				break;}
			default: {
				fprintf(out_cepi_seq,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_seq,"\tEvert");
				//fprintf(out_cepi_sit,"Anc_%d Par_%d",prf_a.mot,prf_p.mot);fprintf(out_cepi_sit,"\tEvert");
				break;}			
			}			
			int j0;
			if(print_overlap==0)j0=noveri+noverp;
			if(print_overlap==1)j0=0;
			if(print_overlap==2)j0=noveri+noverp-noverp/2;
			for(j=j0;j<cepi_len;j++)
			{
				double freq=(double)cepi_seq[k][j]/nseq;
				hist->freq[k][j-j0]=freq;
				fprintf(out_cepi_seq,"\t%f",freq);
				//fprintf(out_cepi_sit,"\t%f",(double)cepi_sit[k][j]/nseq);
			}
			fprintf(out_cepi_seq,"\n");
			//fprintf(out_cepi_sit,"\n");
		}
	}
	fclose(out_cepi_seq);
	*/
	//fclose(out_cepi_sit);
	//fclose(out_over_spac);
	//	printf("Debug - At the end! delete\n");
	for(j=0;j<4;j++)delete [] cepi_seq[j];
//	for(j=0;j<4;j++)delete [] cepi_sit[j];
	for(j=0;j<4;j++)delete [] cepi_sit_one[j];
	delete [] cepi_seq_dir;
	//printf("Debug - At the end! return\n");
	sample->two_sites=nseq_join_sites;
	sample->any=rec_seq;
	sample->full=n_peak_full_over;
	sample->partial=n_peak_part_over;
	sample->over=n_peak_over;
	sample->spacer=n_peak_spacer;
	return 1;	
}

