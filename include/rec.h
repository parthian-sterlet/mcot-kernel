
int rec0(char *file_pro, int nseq_pro, profile *real)
{
	int n; 			
	FILE *in;
	char str[SEQLEN];
	int nseq_rec=0;
	if ((in = fopen(file_pro, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", file_pro);
		return -1;
	}
	char sym = '>';	
	fgets(str, sizeof(str), in);
	for(n=0;n<nseq_pro;n++)real->nsit[n] = 0;
	n = 0;
	while (fgets(str, sizeof(str), in) != NULL)
	{		
		char c = str[0];
		if (c == sym) 
		{ 
			n++; 			
			continue;
		}		
		if(isdigit(c))real->nsit[n]++;															
	}	
	for (n = 0; n < nseq_pro; n++)
	{
		if (real->nsit[n] > 0)nseq_rec++;
	}
	fclose(in);
	return 1;
}
int rec1(char *file_pro,int nseq_pro, profile *real)
{
	int n, x;
	FILE *in;
	char str[SEQLEN];
	int nseq_rec = 0;
	if ((in = fopen(file_pro, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", file_pro);
		return -1;
	}
	char sym = '>';
	fgets(str, sizeof(str), in);	
	n = x = 0;
	char sep = '\t';
	while (fgets(str, sizeof(str), in) != NULL)
	{
		char c = str[0];
		if (c == sym)
		{
			n++;
			x = 0;
			continue;
		}
		if (isdigit(c))
		{
			real->sta[n][x] = atoi(str);
			char s[30];
			int test = UnderStol(str, 1, s, sep);
			if (test == -1) { printf("Wrong format %s\n", str); return -1;}
			real->sco[n][x] = atof(s);
			test = UnderStol(str, 2, s, sep);
			if (test == -1) { printf("Wrong format %s\n", str); return -1;}
			real->cep[n][x] = s[0];						
			x++;
		}
	}
	fclose(in);
	return 1;

}
