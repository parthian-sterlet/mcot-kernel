char TransStr(char x)
{
  int c=int(x);
  if(c<97) x=char(c+32);	
   return x;
}
int fasta_to_plain0(char *file_in_fasta, int &length_fasta_max, int &nseq_fasta)
{
	char head[200];		
	FILE *in;

	if((in=fopen(file_in_fasta,"rt"))==NULL)
	{
		printf("Input file %s can't be opened\n",file_in_fasta);
		return -1;
	}
	char c, symbol = '>';	
	nseq_fasta=0;
	int len=0;
	//double sum_len=0;
	int fl=1;
	length_fasta_max=0;
	do  
	{
		c=getc(in);
		if(c==EOF)fl=-1;
		if(c==symbol || fl==-1)
		{
			if(nseq_fasta>0)
			{
				if(len>length_fasta_max)length_fasta_max=len;
				if(len>SEQLEN)
				{
					printf("Sequence N %d too long... %d nt\n",nseq_fasta,len);
					return -1;
				}
			}
			if(fl!=-1)
			{
				nseq_fasta++;
				len=0;
			}
			if(fl==1)
			{
				fgets(head,sizeof(head),in);			
				continue;
			}
		}
		if(c=='\n')continue;	
		if(c=='\t')continue;	
		if(c==' ')continue;	
		if(c=='\r')continue;	
		len++;
	}
	while(fl==1);
	fclose(in);
	return 1;
}
int fasta_to_plain1(char *file_in_fasta, int length_fasta_max, int nseq_fasta, char ***seq, int *peak_len)
{
	int i, fl=1, n=0, len=0;
	char head[200];		
	char c, symbol = '>';
	char alfavit4[]="acgtnACGTN";
	FILE *in;
	if((in=fopen(file_in_fasta,"rt"))==NULL)
	{
		printf("Input file %s can't be opened\n",file_in_fasta);
		return -1;
	}	
	do  
	{
		c=getc(in);
		if(c==EOF)
		{
			fl=-1;
		}
		if(c==symbol || fl==-1)
		{		
			if(len>0)
			{	
				peak_len[n]=len;
				seq[0][n][len]='\0';
				strncpy(seq[1][n],seq[0][n],len);
				seq[1][n][len]='\0';
				ComplStr(seq[1][n]);				
				{
					if(fl!=-1)
					{
						n++;						
						len=0;
					}
				}
			}
			else
			{
				if(n>0 && fl!=-1)
				{
					printf("Peak length error! peak %d\n",n+1);
					return -1;
				}
			}
			if(fl==-1)break;
			if(fl==1)
			{
				fgets(head,sizeof(head),in);			
				continue;
			}
		}
		if(c=='\n')continue;	
		if(c=='\t')continue;	
		if(c==' ')continue;	
		if(c=='\r')continue;	
		if(strchr(alfavit4,c)!=NULL)
		{
			c=TransStr(c);
			seq[0][n][len++]=c;
		}
		else seq[0][n][len++]='n';
	}
	while(fl==1);
	fclose(in);	
	return 1;
}