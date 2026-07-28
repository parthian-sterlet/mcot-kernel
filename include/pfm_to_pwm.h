double pfm[MATLEN][OLIGNUM];

int UnderStol(char* str, int nstol, char* ret, size_t size, char sep)
{
	memset(ret, 0, size);
	int p1, p2, len;
	if (nstol == 0)
	{
		p2 = StrNStr(str, sep, 1);
		if (p2 == -1)p2 = (int)strlen(str);
		if (p2 == 0) return -1;
		strncpy(ret, str, p2);
		ret[p2] = '\0';
		return 1;
	}
	else
	{
		p1 = StrNStr(str, sep, nstol);
		p2 = StrNStr(str, sep, nstol + 1);
		if (p2 == -1)
		{
			p2 = (int)strlen(str);
			if (p2 == 0) return -1;
		}
		if (p1 == -1 || p2 == -1) return -1;
		len = p2 - p1 - 1;
		strncpy(ret, &str[p1 + 1], len);
		ret[len] = '\0';
		return 1;
	}
}
int pfm_to_pwm_one(char *file_pfm, matrices *mat)
{
	char d[500], head[500], s[500];
	int i, j, len = 0;
	int shift_col;//=atoi(argv[4]);
	int mtype;

	FILE* in_pfm;
	if ((in_pfm = fopen(file_pfm, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file_pfm);
		return -1;
	}
	int alfabet = 0;
	fgets(head, sizeof(head), in_pfm);
	fgets(head, sizeof(head), in_pfm);
	if (strchr("ATGC", head[0]) != NULL)
	{// jaspar
		mtype = 1;
		shift_col = 1;
		alfabet++;
		while (fgets(head, sizeof(head), in_pfm) != NULL)
		{
			if (strchr("ATGC", head[0]) != NULL)alfabet++;
		}
		if (alfabet != 4 && alfabet != 16)
		{
			printf("Reading error in Jaspar matrix file %s !\n", file_pfm);
			return -1;
		}
	}
	else
	{//homer, cis-bp
		mtype = 0;
		//shift_col=0;
	}
	if (mtype == 0)
	{
		DelChar(head, '\n');
		DelChar(head, '\r');
		int headlen = (int)strlen(head);
		if (head[headlen - 1] == '\t')
		{
			head[headlen - 1] = '\0';
			headlen--;
		}
		int counttab = 0;
		for (i = 0; i < headlen; i++)
		{
			if (head[i] == '\t')counttab++;
		}
		if (counttab == 4 || counttab == 16)
		{
			shift_col = 1;
			alfabet = counttab;
		}
		else
		{
			if (counttab == 3 || counttab == 15)
			{
				shift_col = 0;
				alfabet = counttab + 1;
			}
			else
			{
				printf("Reading errorin Cisbp/Homer matrix file %s !\n", file_pfm);
				return -1;
			}
		}
	}
	rewind(in_pfm);
	int test;
	char sep = '\t';
	int olen = 0;
	double olen_perf = 0;
	fgets(head, sizeof(head), in_pfm);
	if (mtype == 0)//homer, cis-bp
	{
		for (i = 0; i < MATLEN; i++)
		{
			if (fgets(d, sizeof(d), in_pfm) != NULL)
			{
				char c = d[0];
				if (isdigit(c) || (strchr("-ATGC", c) != 0))olen++;
			}
			else break;
		}
		rewind(in_pfm);
		mat->mem_in(olen);
		fgets(head, sizeof(head), in_pfm);
		for (i = 0; i < MATLEN; i++)
		{
			if (fgets(d, sizeof(d), in_pfm) != NULL)
			{
				DelChar(d, '\n');
				DelChar(d, '\r');
				char c = d[0];
				if (isdigit(c) || (strchr("-ATGC", c) != 0))
				{
					//	 printf("%s",d[i]);						
					for (j = 0; j < alfabet; j++)
					{
						test = UnderStol(d, j + shift_col, s, sizeof(s), sep);
						if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
						mat->fre[i][j] = atof(s);						 
					}
				}
				else break;
			}
			else break;
		}
	}
	else
	{// jaspar
		fgets(d, sizeof(d), in_pfm);
		int dlen = (int)strlen(d);
		olen = 0;
		for (i = 1; i < dlen; i++)
		{
			if (d[i - 1] == '\t' && isdigit(d[i]))olen++;
		}
		rewind(in_pfm);
		mat->mem_in(olen);
		fgets(head, sizeof(head), in_pfm);
		for (i = 0; i < alfabet; i++)
		{
			if (fgets(d, sizeof(d), in_pfm) != NULL)
			{
				DelChar(d, '\n');
				char c = d[0];
				if (isdigit(c) || (strchr("-ATGC", c) != 0))
				{
					//	 printf("%s",d[i]);	
					for (j = 0; j < olen; j++)
					{
						test = UnderStol(d, j + shift_col, s, sizeof(s), sep);
						if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
						mat->fre[j][i] = atof(s);						
					}
				}
				else break;
			}
			else break;
		}
	}
	fclose(in_pfm);
	//logodds score
	double pse = 0.25;
	double pse4 = 1;
	double vych = log10(pse);
	int olen1 = olen - 1;
	for (i = 0; i < olen; i++)
	{
		double sum = 0;
		for (j = 0; j < alfabet; j++)sum +=  mat->fre[i][j];
		int alfabet1 = alfabet - 1;
		for (j = 0; j < alfabet; j++)
		{
			double count = mat->fre[i][j];
			double ves = (count + pse) / (sum + pse4);
			mat->wei[i][j] = log10(ves) - vych;
		}
	}
	return olen;
}
int pfm_to_pwm(char* file_pfm, matrices* mat)
{
	char d[500], head[500], s[500];
	int i, j, len = 0;
	int shift_col;//=atoi(argv[4]);
	int mtype;

	FILE* in_pfm;
	if ((in_pfm = fopen(file_pfm, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file_pfm);
		return -1;
	}
	int alfabet = 0;
	fgets(head, sizeof(head), in_pfm);
	fgets(head, sizeof(head), in_pfm);
	if (strchr("ATGC", head[0]) != NULL)
	{// jaspar
		mtype = 1;
		shift_col = 1;
		alfabet++;
		while (fgets(head, sizeof(head), in_pfm) != NULL)
		{
			if (strchr("ATGC", head[0]) != NULL)alfabet++;
		}
		if (alfabet != 4 && alfabet != 16)
		{
			printf("Reading error in Jaspar matrix file %s !\n", file_pfm);
			return -1;
		}
	}
	else
	{//homer, cis-bp
		mtype = 0;
		//shift_col=0;
	}
	if (mtype == 0)
	{
		DelChar(head, '\n');
		DelChar(head, '\r');
		int headlen = (int)strlen(head);
		if (head[headlen - 1] == '\t')
		{
			head[headlen - 1] = '\0';
			headlen--;
		}
		int counttab = 0;
		for (i = 0; i < headlen; i++)
		{
			if (head[i] == '\t')counttab++;
		}
		if (counttab == 4 || counttab == 16)
		{
			shift_col = 1;
			alfabet = counttab;
		}
		else
		{
			if (counttab == 3 || counttab == 15)
			{
				shift_col = 0;
				alfabet = counttab + 1;
			}
			else
			{
				printf("Reading errorin Cisbp/Homer matrix file %s !\n", file_pfm);
				return -1;
			}
		}
	}
	rewind(in_pfm);
	int test;
	char sep = '\t';
	int olen = 0;
	double olen_perf = 0;
	fgets(head, sizeof(head), in_pfm);
	if (mtype == 0)//homer, cis-bp
	{
		for (i = 0; i < MATLEN; i++)
		{
			if (fgets(d, sizeof(d), in_pfm) != NULL)
			{
				char c = d[0];
				if (isdigit(c) || (strchr("-ATGC", c) != 0))olen++;
			}
			else break;
		}
		rewind(in_pfm);
		mat->mem_in(olen);
		fgets(head, sizeof(head), in_pfm);
		for (i = 0; i < MATLEN; i++)
		{
			if (fgets(d, sizeof(d), in_pfm) != NULL)
			{
				DelChar(d, '\n');
				DelChar(d, '\r');
				char c = d[0];
				if (isdigit(c) || (strchr("-ATGC", c) != 0))
				{
					//	 printf("%s",d[i]);						
					for (j = 0; j < alfabet; j++)
					{
						test = UnderStol(d, j + shift_col, s, sizeof(s), sep);
						if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
						mat->fre[i][j] = atof(s);
					}
				}
				else break;
			}
			else break;
		}
	}
	else
	{// jaspar
		fgets(d, sizeof(d), in_pfm);
		int dlen = (int)strlen(d);
		olen = 0;
		for (i = 1; i < dlen; i++)
		{
			if (d[i - 1] == '\t' && isdigit(d[i]))olen++;
		}
		rewind(in_pfm);
		mat->mem_in(olen);
		fgets(head, sizeof(head), in_pfm);
		for (i = 0; i < alfabet; i++)
		{
			if (fgets(d, sizeof(d), in_pfm) != NULL)
			{
				DelChar(d, '\n');
				char c = d[0];
				if (isdigit(c) || (strchr("-ATGC", c) != 0))
				{
					//	 printf("%s",d[i]);	
					for (j = 0; j < olen; j++)
					{
						test = UnderStol(d, j + shift_col, s, sizeof(s), sep);
						if (test == -1) { printf("Wrong format %s\n", d); exit(1); }
						mat->fre[j][i] = atof(s);
					}
				}
				else break;
			}
			else break;
		}
	}
	fclose(in_pfm);
	//logodds score
	double pse = 0.25;
	double pse4 = 1;
	double vych = log10(pse);
	int olen1 = olen - 1;
	for (i = 0; i < olen; i++)
	{
		double sum = 0;
		for (j = 0; j < alfabet; j++)sum += mat->fre[i][j];
		int alfabet1 = alfabet - 1;
		for (j = 0; j < alfabet; j++)
		{
			double count = mat->fre[i][j];
			double ves = (count + pse) / (sum + pse4);
			mat->wei[i][j] = log10(ves) - vych;
		}
	}
	return olen;
}