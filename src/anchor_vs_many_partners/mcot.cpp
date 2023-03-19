#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 12000 // max length of peak in input fasta
#define MATLEN 50 //max matrix length
#define SPACLEN 100 //max upper bound of spacer length
#define ARGLEN 300 //max argv length
#define OLIGNUM 4// di 16 mono 4
#define NUM_THR 5 //4islo porogov
#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);

//return n-th occurrence of a certain symbol c in a string str
int StrNStr(char *str, char c, int n)
{
	int i, len = strlen(str);
	int k = 1;
	for (i = 0; i<len; i++)
	{
		if (str[i] == c)
		{
			if (k == n)return i;
			k++;
		}
	}
	return -1;
}
// menyaet registr stroki
char *TransStr(char *d)
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i<lens; i++)
	{
		c = int(d[i]);
		if (c<97) d[i] = char(c + 32);
		//else break;
	}
	return(d);
}
char *TransStrBack(char *d)
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i<lens; i++)
	{
		c = int(d[i]);
		if (c >= 97) d[i] = char(c - 32);
		//else break;
	}
	return(d);
}
//udalenie simvola iz stroki
void DelChar(char *str, char c)
{
	int i, lens, size;

	size = 0;
	lens = strlen(str);
	for (i = 0; i<lens; i++)
	{
		if (str[i] != c)str[size++] = str[i];
	}
	str[size] = '\0';
}
// ras4et 4asot oligonukleotidov po stroke (zdes' - nukleotidov)
void GetSostPro(char *d, int word, int *sost)
{
	int i, j, k, i_sost, let;
	char letter[] = "acgt";
	int ten[6] = { 1, 4, 16, 64, 256, 1024 };
	int lens = strlen(d);
	int size = 1;
	for (k = 0; k<word; k++)size *= 4;
	for (i = 0; i<size; i++)sost[i] = 0;
	for (i = 0; i<lens - word + 1; i++)
	{
		i_sost = 0;
		for (j = word - 1; j >= 0; j--)
		{
			for (k = 0; k<4; k++)
			{
				if (d[i + j] == letter[k]){ let = k; break; }
			}
			i_sost += ten[word - 1 - j] * let;
		}
		sost[i] = i_sost;
	}
}
//komplementaciya stroki
int ComplStr(char *d)
{
	char *d1;
	int i, len;
	len = strlen(d);
	d1 = new char[len + 1];
	if (d1 == NULL)
	{
		fprintf(stderr, "Error: Out of memory...");
		return 0;
	}
	strcpy(d1, d);
	//	memset(d,0,sizeof(d));
	for (i = 0; i<len; i++)
	{
		switch (d1[len - i - 1])
		{
		case 'a':{d[i] = 't'; break; }
		case 't':{d[i] = 'a'; break; }
		case 'c':{d[i] = 'g'; break; }
		case 'g':{d[i] = 'c'; break; }
		case 'A':{d[i] = 'T'; break; }
		case 'T':{d[i] = 'A'; break; }
		case 'C':{d[i] = 'G'; break; }
		case 'G':{d[i] = 'C'; break; }
		case 'N':{d[i] = 'N'; break; }
		case 'n':{d[i] = 'n'; break; }
		default: d[i] = 'n';
		}
	}
	delete[] d1;
	return 1;
}
void Mix(int *a, int *b)
{
	int buf = *a;
	*a = *b;
	*b = buf;
}
void Mix(char *a, char *b)
{
	char buf = *a;
	*a = *b;
	*b = buf;
}
void BigMix1(char *d)//me6alka
{
	int r;
	int len = strlen(d);
	for (r = 0; r<len - 1; r++) Mix(&d[r], &d[1 + r + (rand() % (len - 1 - r))]);
}
void BigMix1(int *d1, int len) // pereme6ivanie stroki
{
	int r, s;
	for (r = 0; r<len; r++)
	{
		s = rand() % len;
		if (s != r)Mix(&d1[r], &d1[s]);
	}
}
#include "fasta_to_plain.h" //input = peaks (fasta), output = plain format OR lengths of peaks
#include "pwm_iz_pwm_thr_dist.h" // full list of thresholds for PWM -> FP-rates po genomu
#include "select_thresholds_from_pvalues.h" //fp-rates po genomu -> vybor pyati porogov po pyati fixed FP rates

int bad_matrix_mm_core[] = { 172, 186, 192, 261, 287, -1 };
int bad_matrix_hs_core[] = { 170, 183, 190, 253, 280, 324, -1 };
int bad_matrix_dapseq[] = { 2, 102, 183, 184, 186, 207, 212, 213, 217, 286, 299, 300, 303, 439, -1 };
int bad_matrix_mm_full[] = { 40, 61, 64, 99, 107, 108, 158, 169, 243, 253, 264, 272, 277, 283, 330, 380, 381, 424, 428, 429, 499, 529, -1 };
int bad_matrix_hs_full[] = { 80, 84, 129, 139, 141, 142, 212, 230, 341, 360, 369, 379, 389, 394, 408, 443, 509, 510, 558, 564, 565, 634, 647, 719, -1 };
int hs_core_thr_count[] = { 221, 1984, 1086, 2108, 2206, 177, 2194, 586, 221, 2174, 1131, 861, 1043, 2128, 2300, 57, 1708, 697, 2206, 2015, 1815, 1787, 563, 372, 834, 2046, 18, 1226, 1479, 1782, 757, 1262, 1715, 1524, 2525, 2266, 1557, 1142, 1151, 1705, 3699, 3195, 1792, 901, 272, 926, 1846, 391, 311, 1405, 190, 1248, 1049, 2122, 2390, 2684, 2303, 1855, 2299, 2140, 984, 1848, 223, 1427, 2880, 76, 60, 2407, 2079, 1708, 719, 800, 2755, 356, 1291, 2873, 1956, 184, 954, 2580, 59, 57, 2397, 1770, 1200, 1342, 1789, 3383, 78, 2272, 70, 1220, 146, 947, 835, 244, 36, 41, 68, 181, 2466, 2374, 2486, 437, 298, 2140, 2398, 184, 230, 281, 1445, 3689, 73, 49, 2917, 2430, 2212, 2545, 2553, 2362, 1189, 3001, 1903, 173, 91, 979, 761, 571, 208, 30, 324, 585, 264, 4, 1414, 2192, 2907, 1963, 2760, 120, 2850, 169, 105, 127, 782, 845, 739, 2541, 1684, 255, 2824, 2666, 95, 1494, 2262, 41, 2921, 1195, 1190, 2416, 1375, 1974, 709, 2063, 3748, 3078, 308, 2603, 428, 6, 1305, 1622, 1385, 973, 1429, 3736, 309, 3711, 2003, 1308, 1067, 1621, 4, 2714, 2177, 1155, 2212, 268, 49, 12, 3142, 1599, 40, 47, 168, 2215, 2006, 2277, 13, 620, 2990, 2793, 2301, 186, 9, 30, 1068, 320, 2160, 95, 1956, 2117, 2329, 2712, 3083, 2630, 1451, 1871, 60, 79, 778, 2514, 4246, 2130, 1917, 915, 22, 11474, 3073, 2875, 2397, 2689, 3326, 1261, 286, 2057, 727, 84, 791, 1578, 821, 2273, 696, 1971, 2314, 2492, 2288, 2267, 2195, 1980, 784, 2307, 7, 2406, 2523, 2082, 2102, 1349, 3633, 3150, 3798, 4944, 2307, 1335, 1396, 805, 1914, 314, 2007, 181, 2974, 53, 2243, 1722, 992, 581, 1546, 972, 2061, 9, 239, 2325, 597, 1130, 511, 777, 15, 2504, 4693, 2978, 2470, 2546, 2983, 3383, 523, 1099, 2859, 19, 241, 1234, 2954, 2047, 1360, 583, 633, 626, 514, 1828, 2733, 198, 1051, 330, 2010, 1284, 1699, 2587, 1666, 1887, 222, 1721, 420, 333, 24, 5, 2015, 4111, 2356, 2375, 2614, 1941, 2190, 3523, 2295, 2542, 2938, 2118, 2303, 167, 2151, 684, 1265, 57, 1680, 186, 2504, 1523, 2533, 181, 33, 2539, 3485, 3456, 5291, 2724, 3596, 3488, 3876, 2348, 4384, 1206, 4724, 3072, 2599, 2652, 2190, 2775, 2934, 2831, 4225, 3900, 2746, 2357, 1740, 2981, 1027, 1708, 2299, 5849, 483, 3373, 5375, 4136, 3444, 2635, 2844, 2066, 2105, 3231, 5404, 2120, 3826, 2882, 4359, 3760, 2898, 1545, 2778, 3477, 4568, 2741, 2230, 2391 };
int hs_full_thr_count[] = { 221, 1984, 1086, 737, 2043, 2108, 2330, 947, 2206, 177, 2194, 1025, 3239, 586, 221, 537, 2531, 2174, 56, 1131, 861, 10799, 564, 1043, 2128, 2300, 1727, 57, 1708, 697, 2072, 2434, 889, 183, 2206, 922, 2015, 1815, 1787, 1041, 563, 4850, 971, 372, 3566, 834, 102, 2046, 520, 92, 1735, 734, 18, 1226, 1479, 1782, 757, 1262, 1715, 1120, 3432, 1524, 2525, 2266, 1401, 1557, 2112, 1698, 163, 2744, 1142, 7262, 1703, 1151, 1705, 3699, 3195, 1792, 709, 13, 901, 225, 1850, 13, 272, 2530, 604, 1819, 2686, 2786, 1419, 813, 926, 1057, 1846, 391, 311, 1405, 1810, 190, 1248, 1049, 1040, 113, 2122, 2390, 2083, 278, 177, 2684, 2303, 1855, 2299, 2140, 984, 2189, 1848, 1634, 1558, 1180, 223, 1427, 2880, 76, 60, 2407, 256, 2079, 4, 2944, 1708, 719, 800, 2755, 741, 356, 1291, 39, 3, 2873, 1720, 1956, 184, 954, 594, 2580, 1394, 59, 57, 2397, 1770, 1200, 1342, 1789, 2942, 3383, 1656, 2333, 3900, 1435, 200, 366, 1906, 78, 2272, 70, 1220, 1490, 146, 7132, 947, 835, 244, 36, 2282, 41, 68, 40, 181, 505, 2466, 2374, 1069, 2486, 792, 437, 298, 2124, 2140, 2652, 2341, 2037, 1305, 2398, 1017, 184, 230, 1186, 62, 281, 2029, 3086, 245, 1548, 1176, 1445, 844, 772, 1751, 2122, 2247, 3, 32, 3689, 1544, 2987, 2476, 2104, 1944, 3995, 73, 526, 49, 2917, 2430, 1114, 106, 3224, 576, 2, 3027, 4021, 760, 1160, 2212, 2545, 599, 2553, 2362, 1189, 9343, 3001, 48, 1903, 3906, 2189, 173, 91, 979, 1118, 761, 907, 368, 1750, 571, 70, 208, 3339, 1826, 30, 1224, 324, 585, 1616, 862, 1220, 701, 955, 1918, 264, 86, 354, 281, 414, 884, 12, 592, 13, 578, 4, 1414, 2192, 2907, 1963, 2760, 5504, 120, 2850, 169, 1552, 634, 105, 948, 2819, 127, 2589, 782, 845, 739, 2541, 1167, 1913, 1684, 255, 537, 1932, 2824, 2680, 2666, 95, 1494, 2262, 41, 2921, 1162, 1195, 1190, 2416, 663, 394, 3224, 3020, 3253, 2679, 1375, 1974, 20, 2187, 709, 2063, 68, 3748, 18, 3078, 1859, 308, 2603, 112, 428, 2419, 6, 1305, 1622, 1385, 973, 1429, 76, 3736, 971, 2773, 761, 423, 3186, 309, 654, 1393, 8048, 1590, 3009, 3, 3711, 2003, 842, 1308, 2435, 151, 1067, 1621, 4, 2636, 2714, 634, 2177, 1155, 2212, 936, 268, 49, 12, 3142, 1599, 96, 40, 47, 168, 234, 2215, 2006, 5, 54, 2277, 55, 904, 13, 620, 2990, 2793, 2301, 1601, 186, 13, 835, 9, 30, 1068, 320, 2160, 5, 782, 95, 805, 184, 1956, 5960, 3140, 2117, 167, 2329, 17, 2712, 19, 3083, 19, 2630, 1451, 3375, 1871, 2504, 60, 79, 53, 778, 2514, 4246, 776, 5529, 2130, 194, 1191, 7088, 1182, 1917, 9, 915, 22, 35, 11474, 3073, 2099, 2941, 2875, 576, 2397, 1053, 2689, 14, 2997, 1955, 6649, 1097, 3326, 1261, 977, 1508, 286, 2107, 2057, 727, 141, 84, 94, 791, 731, 748, 1578, 104, 156, 2724, 821, 2197, 143, 2273, 696, 1998, 1971, 2314, 1712, 2865, 2320, 2925, 3053, 2492, 190, 613, 810, 2288, 51, 2064, 2267, 53, 2195, 1980, 6302, 784, 2307, 52, 582, 1657, 6, 7, 2406, 40, 1900, 2523, 116, 2158, 722, 2082, 12, 2379, 1377, 2102, 1349, 3633, 3150, 796, 3798, 1612, 4944, 3186, 2307, 652, 703, 1335, 1396, 3640, 805, 1914, 314, 647, 2007, 182, 181, 2974, 53, 2480, 4233, 3715, 3832, 2243, 1722, 1100, 992, 581, 1546, 972, 2061, 9, 239, 5621, 2325, 991, 3656, 9, 5, 597, 3714, 1130, 1869, 2918, 511, 777, 15, 4236, 2717, 2504, 67, 4693, 1710, 2978, 453, 2470, 2546, 2277, 1884, 2983, 3383, 649, 1713, 523, 1099, 2859, 19, 241, 1234, 2954, 1640, 2047, 1360, 583, 633, 626, 514, 1828, 2733, 1905, 198, 2221, 3143, 31, 1, 7592, 2187, 1051, 330, 2544, 13, 2010, 1284, 1314, 2099, 1699, 1988, 2797, 2587, 1666, 1887, 222, 2039, 1721, 420, 333, 24, 5, 753, 2015, 154, 4111, 2356, 2375, 209, 1629, 2614, 228, 1941, 997, 1, 630, 2190, 3523, 4453, 2366, 2295, 60, 2098, 2542, 235, 2187, 1708, 2938, 225, 750, 2118, 2303, 7612, 167, 2151, 684, 1265, 844, 57, 3561, 7211, 2600, 1680, 186, 838, 1720, 2243, 490, 2504, 1523, 2533, 181, 1746, 33, 2031, 2539, 3081, 3485, 3456, 940, 5291, 2724, 2575, 3596, 3488, 3876, 1785, 2348, 729, 7516, 4384, 1206, 4724, 3072, 385, 2599, 2652, 2190, 1842, 2775, 2934, 2831, 4225, 3900, 71, 2746, 5, 2357, 63, 1740, 3518, 2981, 1027, 1708, 36, 2971, 2299, 1193, 172, 5849, 483, 3373, 5375, 4136, 2567, 3444, 2635, 2844, 2066, 39, 2105, 941, 3231, 5404, 1735, 733, 2120, 3826, 2882, 78, 2104, 1968, 4359, 3760, 85, 836, 2898, 1674, 1545, 2778, 1105, 3477, 4568, 2741, 2114, 4110, 2230, 2391, 786 };
int mm_full_thr_count[] = { 254, 1885, 1084, 2825, 254, 1885, 1084, 2825, 885, 1636, 280, 1910, 1441, 3082, 616, 65, 704, 2824, 59, 338, 887, 62, 2219, 823, 56, 1633, 790, 941, 171, 2050, 756, 2027, 1448, 596, 134, 5997, 117, 2423, 1364, 66, 718, 16, 1170, 1211, 856, 1076, 758, 1330, 2006, 1054, 1491, 2288, 56, 1154, 1575, 1959, 2189, 1124, 1221, 1748, 4150, 3898, 1764, 676, 15, 885, 257, 13, 360, 651, 2904, 2404, 1079, 505, 229, 1894, 2010, 203, 959, 1244, 113, 2358, 1584, 2198, 101, 178, 2604, 2689, 1187, 2240, 2578, 950, 2133, 1442, 239, 1586, 1542, 3098, 227, 2252, 278, 2185, 14, 2126, 618, 2681, 12, 1175, 37, 2827, 1632, 1283, 214, 2164, 2472, 2468, 61, 2227, 768, 188, 1670, 1679, 3773, 1884, 2260, 1276, 116, 332, 2170, 50, 1103, 962, 613, 347, 860, 167, 226, 29, 66, 40, 173, 467, 2585, 2084, 924, 740, 1554, 260, 2043, 301, 2126, 2290, 904, 252, 185, 1297, 61, 1146, 280, 1430, 2195, 3, 32, 3688, 1642, 2161, 6415, 91, 52, 2986, 923, 1095, 2, 3111, 2247, 2370, 1830, 2419, 2291, 1133, 3319, 52, 1908, 144, 90, 938, 741, 347, 1777, 493, 67, 29, 1192, 320, 578, 914, 1921, 77, 614, 12, 12, 9, 1533, 2858, 2824, 1450, 2307, 5439, 112, 2824, 179, 102, 65, 213, 55, 833, 2677, 1344, 1948, 1973, 2733, 2790, 1891, 2012, 2314, 45, 1115, 365, 86, 135, 1347, 2186, 57, 1466, 637, 2858, 246, 4069, 17, 3745, 1962, 1022, 2352, 139, 861, 2410, 10, 1804, 2034, 2197, 1364, 1203, 3409, 265, 8316, 1537, 5, 4123, 4647, 1973, 35, 1214, 2510, 142, 1517, 45, 544, 5, 2741, 713, 2251, 1962, 808, 245, 216, 10, 2038, 2370, 395, 42, 3, 47, 168, 4673, 2920, 2029, 5, 60, 53, 1541, 922, 798, 772, 2926, 2682, 3315, 1559, 274, 669, 276, 32, 1146, 296, 2037, 63, 91, 216, 2009, 6222, 2048, 3215, 2015, 2111, 184, 2248, 20, 2853, 18, 3027, 21, 2498, 1552, 1919, 2488, 42, 65, 53, 682, 1505, 4471, 2081, 175, 6977, 14, 1380, 22, 35, 3022, 2151, 2999, 544, 2502, 1121, 1985, 1473, 2173, 1297, 1537, 866, 1917, 2077, 1945, 2270, 107, 664, 1564, 243, 237, 802, 2165, 2187, 585, 2126, 2217, 3039, 2996, 2589, 70, 669, 2335, 51, 2091, 2190, 179, 2178, 165, 1952, 1880, 40, 2750, 2435, 13, 582, 6, 7, 2353, 40, 1830, 4069, 3050, 2192, 1231, 898, 2116, 17, 2130, 1446, 4556, 3306, 845, 3150, 973, 3973, 303, 541, 2064, 1429, 1318, 3334, 1421, 1255, 1356, 3538, 3184, 235, 214, 4342, 143, 1652, 65, 1944, 1483, 975, 40, 680, 826, 2191, 10, 702, 2304, 341, 8, 3, 606, 3452, 2115, 1138, 484, 624, 14, 2399, 306, 3394, 149, 3125, 646, 2412, 2469, 2257, 3015, 735, 2126, 2616, 3582, 1797, 736, 1186, 2571, 18, 206, 1139, 1602, 2334, 2058, 1104, 175, 2072, 501, 509, 1668, 2686, 163, 7831, 2153, 700, 88, 14, 1176, 921, 1506, 2082, 2287, 1950, 2711, 300, 770, 595, 1752, 2077, 1793, 521, 407, 28, 70, 2457, 52, 4266, 1918, 1731, 2504, 600, 1711, 1, 3314, 3673, 2299, 73, 2038, 2842, 250, 5087, 2740, 704, 84, 1835, 207, 921, 1759, 468, 1433, 430, 3001, 25, 37, 2560, 2652, 2975, 4575, 1497, 2938, 2683, 2723, 7, 169, 3513 };
int mm_core_thr_count[] = { 254, 1885, 1084, 2825, 1636, 280, 1910, 616, 65, 2824, 59, 338, 887, 62, 2219, 56, 1633, 790, 941, 2050, 2027, 1448, 596, 134, 117, 2423, 16, 1170, 1211, 856, 1076, 758, 1330, 2006, 1491, 2288, 56, 1575, 1124, 1221, 1748, 4150, 3898, 1764, 676, 885, 257, 360, 651, 2904, 2404, 1079, 505, 229, 1894, 203, 959, 1244, 2358, 1584, 2604, 2689, 1187, 2240, 2578, 950, 1442, 239, 1586, 1542, 3098, 227, 2252, 2185, 2126, 618, 2681, 12, 37, 2827, 1632, 1283, 214, 2164, 2468, 61, 2227, 768, 188, 1670, 1679, 3773, 1276, 2170, 50, 1103, 613, 347, 860, 167, 226, 29, 66, 173, 2585, 2084, 740, 1554, 260, 301, 2290, 252, 185, 1297, 1430, 2195, 3688, 91, 52, 2986, 923, 2247, 2370, 2419, 2291, 1133, 3319, 1908, 144, 90, 938, 741, 493, 29, 320, 578, 9, 1533, 2858, 2824, 1450, 2307, 112, 2824, 179, 102, 65, 213, 55, 833, 2677, 1973, 2733, 2790, 1891, 2012, 2314, 45, 1115, 365, 86, 135, 1347, 2186, 637, 2858, 4069, 3745, 1022, 2352, 861, 10, 1804, 2034, 2197, 1364, 3409, 265, 1537, 4647, 1973, 1214, 2510, 1517, 45, 5, 2741, 2251, 1962, 245, 216, 10, 2038, 2370, 42, 47, 168, 2920, 2029, 60, 53, 922, 798, 772, 2926, 2682, 3315, 1559, 274, 669, 276, 32, 1146, 296, 2037, 91, 2009, 2048, 2015, 2248, 2853, 3027, 2498, 1552, 1919, 42, 65, 682, 4471, 2081, 1380, 22, 35, 3022, 2999, 2502, 2173, 1297, 866, 1945, 2270, 107, 664, 1564, 243, 802, 2187, 585, 2126, 2217, 2589, 2335, 2190, 2178, 165, 1952, 1880, 2750, 2435, 582, 7, 2353, 4069, 2116, 2130, 1446, 4556, 3306, 3150, 3973, 541, 1429, 1318, 1421, 1255, 1356, 3538, 214, 4342, 65, 1944, 1483, 40, 680, 826, 2191, 10, 702, 2304, 2115, 484, 624, 14, 2399, 3394, 3125, 2412, 2469, 3015, 2126, 2616, 3582, 736, 1186, 2571, 18, 206, 1139, 1602, 2058, 1104, 175, 2072, 501, 509, 1668, 2686, 163, 2153, 700, 88, 1176, 921, 1506, 2287, 2711, 770, 595, 1752, 521, 407, 28, 70, 2457, 4266, 2504, 1711, 3314, 3673, 2299, 2038, 2842, 5087, 2740, 84, 207, 1433, 430, 3001, 37, 2560, 2652, 2975, 4575, 2938, 2683, 2723, 3513 };
int dapseq_thr_count[] = { 2700, 2, 2414, 4348, 2350, 7109, 9379, 3776, 1098, 5097, 3683, 3410, 3255, 4096, 4168, 4089, 3625, 3407, 7145, 5043, 3913, 6753, 2759, 1576, 5522, 3067, 4930, 3648, 4978, 8916, 3446, 3625, 2605, 3455, 6749, 5009, 5213, 3375, 3770, 4611, 4410, 8802, 3641, 5019, 6086, 8439, 4250, 4426, 6778, 4255, 4283, 4542, 3746, 3252, 5207, 3930, 4211, 3739, 4596, 4468, 4612, 3536, 6049, 3655, 4654, 3804, 4038, 4088, 3895, 3964, 2604, 2106, 3894, 5255, 3501, 3359, 4338, 3715, 3951, 5055, 4756, 6261, 271, 367, 2717, 131, 132, 1691, 6946, 12027, 7245, 4719, 5236, 4275, 3302, 875, 2289, 3423, 4241, 3079, 1, 1, 7318, 3783, 6230, 1085, 7112, 4159, 7103, 4206, 7702, 9935, 7053, 7745, 7687, 1833, 3925, 4214, 1448, 436, 5280, 4699, 4352, 3391, 4892, 4158, 1843, 5358, 4817, 4740, 3962, 3424, 4541, 3995, 1060, 4818, 4062, 1012, 2587, 1345, 1304, 1046, 5488, 10207, 4939, 2422, 280, 2181, 3580, 2727, 2208, 308, 2992, 2713, 2530, 3011, 2418, 2570, 246, 1716, 2392, 108, 1665, 2450, 2196, 1336, 2773, 2829, 2427, 2505, 4728, 10529, 976, 33, 2321, 3226, 3868, 4024, 5619, 9, 3236, 888, 1, 2, 2042, 4, 11, 2068, 481, 2537, 2252, 3669, 2781, 2042, 3539, 3480, 3771, 745, 3294, 4738, 4241, 2356, 3028, 704, 3438, 2072, 4, 2342, 2273, 258, 1774, 1, 3, 3947, 2428, 9, 2, 1994, 2158, 2248, 2090, 2225, 8640, 4889, 7027, 2899, 2984, 2100, 6287, 6321, 2268, 523, 624, 1647, 2418, 608, 4950, 309, 1912, 2160, 2763, 37, 18, 497, 3058, 2373, 2399, 3212, 2223, 3237, 2757, 246, 2527, 2097, 2543, 254, 183, 290, 1498, 2449, 304, 3039, 210, 251, 241, 269, 1763, 296, 444, 2493, 138, 154, 264, 2835, 392, 3094, 1861, 1456, 3342, 2900, 5019, 3124, 2693, 4364, 1415, 3, 3774, 2342, 3291, 2391, 2582, 2419, 2509, 4050, 4045, 3435, 3017, 2700, 4, 1, 3558, 2768, 1, 3033, 3934, 2286, 2135, 1976, 2094, 2763, 2867, 1718, 2537, 2169, 519, 2979, 3564, 1165, 2054, 2792, 309, 5502, 2158, 564, 2237, 3421, 3734, 3785, 1738, 2364, 2376, 2522, 2893, 2155, 2393, 1191, 2628, 2061, 5989, 2518, 2571, 2104, 2397, 2969, 1451, 2530, 1750, 618, 2467, 1882, 624, 566, 3330, 276, 2007, 2580, 2417, 2894, 1758, 3205, 2724, 2082, 6227, 621, 3021, 1062, 3413, 2117, 1774, 254, 10679, 2040, 2174, 2490, 285, 2519, 2297, 2372, 1441, 59, 125, 2531, 7778, 10252, 6944, 6100, 4284, 7009, 4645, 3848, 6289, 4949, 3698, 2853, 2957, 3045, 3964, 4408, 2234, 8697, 5168, 2743, 3550, 3510, 4620, 3425, 2874, 2828, 3938, 3512, 4662, 3178, 3890, 4124, 4279, 4504, 4817, 2210, 3890, 3682, 3361, 5864, 2533, 1213, 2725, 3351, 7227, 6186, 3406, 2327, 2773, 1681, 2674, 3608, 2563, 3613, 2752, 2615, 1, 99, 5900, 2572, 6529, 2668, 1530, 3474, 2200, 132, 837, 2032, 3142, 4428, 1583, 1105, 1461, 2309, 4566, 20, 1943, 1150, 8466, 3889, 4732, 4537, 3955, 7334, 7566, 8229, 1223, 1422, 3666, 3199, 5294, 2585, 1783, 2106, 2195, 3059, 1584, 8429, 3760, 3557, 2151, 2823, 1604, 3081, 133, 3420, 3170, 741, 3288, 910, 2765, 751, 3050, 3324, 1005, 3141, 3080, 2153, 725, 878, 1071, 4295, 1159, 2895, 3421, 3063, 1611, 334, 2975, 2489, 2254, 1425, 5579, 660, 3677, 3908, 778, 2625, 5370, 992, 2632, 1726, 1621, 381, 418, 1775 };

struct fpr_thr_t {
	int nthr;
	double *thr;
	double *fpr;
	double thr_selected[5];
	int inx_selected[5];
	int mem_ini(int count, double *t, double *f, double *t_s, int *i_s);
	int hs_core_ini(int *n, int mot);
	int hs_full_ini(int *n, int mot);
	int mm_core_ini(int *n, int mot);
	int mm_full_ini(int *n, int mot);
	int dapseq_ini(int *n, int mot);
	void mem_out(void);
} motifs;

int fpr_thr_t::mem_ini(int count, double *t, double *f, double *t_s, int *i_s)
{
	nthr = count - 1;
	thr = new double[nthr];
	if (thr == NULL) return -1;
	fpr = new double[nthr];
	if (fpr == NULL) return -1;
	int i;
	for (i = 0; i<nthr; i++)thr[i] = t[i];
	for (i = 0; i<nthr; i++)fpr[i] = f[i];
	for (i = 0; i<5; i++)thr_selected[i] = t_s[i];
	for (i = 0; i<5; i++)inx_selected[i] = i_s[i];
	return 1;
}
void fpr_thr_t::mem_out(void)
{
	delete[] thr;
	delete[] fpr;
}

#include "fpr_thr_hs_core.h"
#include "fpr_thr_hs_full.h"
#include "fpr_thr_mm_core.h"
#include "fpr_thr_mm_full.h"
#include "fpr_thr_dapseq.h"

// position frequency mattrix (PFM), position weight matrix (PWM)
struct matrices {
	int len;
	double min;
	double raz;
	double **wei;
	double **fre;
	void init_fre(int i, double a, double b, double c, double d);
	void init_wei(int i, double a, double b, double c, double d);
	void init_dapseq(int num);
	void init_hs_core(int num);
	void init_hs_full(int num);
	void init_mm_core(int num);
	void init_mm_full(int num);
	void get_copy(matrices *a);
	int mem_in(int len);
	void mem_out(int len);
	void norm(void);
};

int matrices::mem_in(int length)
{
	int i;
	len = length;
	wei = new double *[len];
	if (wei == NULL) return -1;
	for (i = 0; i<len; i++)
	{
		wei[i] = new double[OLIGNUM];
		if (wei[i] == NULL) return -1;
	}
	fre = new double *[len];
	if (fre == NULL) return -1;
	for (i = 0; i<len; i++)
	{
		fre[i] = new double[OLIGNUM];
		if (fre[i] == NULL) return -1;
	}
	return 1;
}
void matrices::mem_out(int len)
{
	int i;
	for (i = 0; i<len; i++) delete[] wei[i];
	delete[] wei;
	for (i = 0; i<len; i++) delete[] fre[i];
	delete[] fre;
}
void matrices::norm(void)
{
	int i, j;
	min = raz = 0;
	for (i = 0; i<len; i++)
	{
		double pwmmin = 100;
		double pwmmax = -100;
		for (j = 0; j<OLIGNUM; j++)
		{
			double w = wei[i][j];
			if (w<pwmmin)pwmmin = w;
			if (w>pwmmax)pwmmax = w;
		}
		raz += pwmmax;
		min += pwmmin;
	}
	raz -= min;
}

void matrices::init_wei(int i, double a, double b, double c, double d)
{
	wei[i][0] = a;
	wei[i][1] = b;
	wei[i][2] = c;
	wei[i][3] = d;
}
void matrices::init_fre(int i, double a, double b, double c, double d)
{
	fre[i][0] = a;
	fre[i][1] = b;
	fre[i][2] = c;
	fre[i][3] = d;
}
void matrices::get_copy(matrices *a)
{
	a->mem_in(len);
	a->len = len;
	a->min = min;
	a->raz = raz;
	int i, j;
	for (i = 0; i<len; i++)
	{
		for (j = 0; j<OLIGNUM; j++)
		{
			a->fre[i][j] = fre[i][j];
			a->wei[i][j] = wei[i][j];
		}
	}
}
#include "pfm_to_pwm.h" //conversion PFM -> PWM

#include "hocomoco_pwm_hs_core.h" //PWMs, human hocomoco core collection
#include "hocomoco_pwm_hs_full.h" //PWMs, human hocomoco full collection
#include "hocomoco_pwm_mm_core.h" //PWMs, mouse hocomoco core collection
#include "hocomoco_pwm_mm_full.h" //PWMs, mouse hocomoco full collection
#include "dapseq_pwm.h" //PWMs, arabidopsis dapseq collection

//profil' raspoznannuh saytov po pyati porogam i odin slitiy profil' (the most permissive threshold)
struct profile {
	int mot;//motif num
	int nam;//nomer poroga
	int nseq;
	int *nsit;//4islo saytov
	int nsit_all;
	int nseq_rec;
	int **sta;//nthr nsit
	char **cep;
	int **cel;// popadanie v interval porogov
	double **sco;
	double **pv;
	int mem_in_sta(void);
	int mem_in_cep(void);
	int mem_in_cel(void);
	int mem_in_sco(void);
	int mem_in_pv(void);
	int mem_in_nsit(void);
	void mem_out_sta(void);
	void mem_out_cep(void);
	void mem_out_cel(void);
	void mem_out_sco(void);
	void mem_out_nsit(void);
	void mem_out_pv(void);
	int get_copy_rand(profile *a, int height);
	int clear_real(void);
	int fprintf_pro(char *mot_db, double thr, char *mode);//mot_db real/rand
	void count_sites(void);
	int test(void);
};
int profile::mem_in_sta(void)
{
	int i;
	sta = new int*[nseq];
	if (sta == NULL) return -1;
	for (i = 0; i<nseq; i++)
	{
		if (nsit[i]>0)
		{
			sta[i] = new int[nsit[i]];
			if (sta[i] == NULL) return -1;
		}
	}
	return 1;
}
int profile::mem_in_cep(void)
{
	int i;
	cep = new char*[nseq];
	if (cep == NULL) return -1;
	for (i = 0; i<nseq; i++)
	{
		if (nsit[i]>0)
		{
			cep[i] = new char[nsit[i]];
			if (cep[i] == NULL) return -1;
		}
	}
	return 1;
}
int profile::mem_in_cel(void)
{
	int i;
	cel = new int*[nseq];
	if (cep == NULL) return -1;
	for (i = 0; i<nseq; i++)
	{
		if (nsit[i]>0)
		{
			cel[i] = new int[nsit[i]];
			if (cel[i] == NULL) return -1;
		}
	}
	return 1;
}
int profile::mem_in_sco(void)
{
	int i;
	sco = new double*[nseq];
	if (sco == NULL) return -1;
	for (i = 0; i<nseq; i++)
	{
		if (nsit[i]>0)
		{
			sco[i] = new double[nsit[i]];
			if (sco[i] == NULL) return -1;
		}
	}
	return 1;
}
int profile::mem_in_pv(void)
{
	int i;
	pv = new double*[nseq];
	if (pv == NULL) return -1;
	for (i = 0; i<nseq; i++)
	{
		if (nsit[i]>0)
		{
			pv[i] = new double[nsit[i]];
			if (pv[i] == NULL) return -1;
		}
	}
	return 1;
}
void profile::mem_out_sta(void)
{
	int i;
	for (i = 0; i<nseq; i++)
	{
		if (nsit[i]>0)delete[] sta[i];
	}
	delete[] sta;
	sta = NULL;
}
void profile::mem_out_cep(void)
{
	int i;
	for (i = 0; i<nseq; i++)if (nsit[i]>0)delete[] cep[i];
	delete[] cep;
	cep = NULL;
}
void profile::mem_out_cel(void)
{
	int i;
	for (i = 0; i<nseq; i++)if (nsit[i]>0)delete[] cel[i];
	delete[] cel;
	cel = NULL;
}
void profile::mem_out_sco(void)
{
	int i;
	for (i = 0; i<nseq; i++)if (nsit[i]>0)delete[] sco[i];
	delete[] sco;
	sco = NULL;
}
void profile::mem_out_pv(void)
{
	int i;
	for (i = 0; i<nseq; i++)if (nsit[i]>0)delete[] pv[i];
	delete[] pv;
	pv = NULL;
}
int profile::mem_in_nsit(void)
{
	nsit = new int[nseq];
	if (nsit == NULL) return -1;
	int i;
	for (i = 0; i<nseq; i++)nsit[i] = 0;
	return 1;
}
void profile::mem_out_nsit(void)
{
	delete[] nsit;
}
int profile::get_copy_rand(profile *a, int height)
{
	int i, j, h, z, ini;
	a->mot = mot;
	a->nam = nam;
	a->nseq = height*nseq;
	a->nsit_all = height*nsit_all;
	a->nseq_rec = height*nseq_rec;
	z = 0;
	for (i = 0; i<nseq; i++)
	{
		int nsi = nsit[i];
		for (h = 0; h<height; h++)
		{
			a->nsit[z++] = nsi;
		}
	}
	ini = a->mem_in_sta();
	if (ini == -1){ fprintf(stderr,"Error: Not enough memory...\n"); return -1; }
	ini = a->mem_in_cep();
	if (ini == -1){ fprintf(stderr,"Error: Not enough memory...\n"); return -1; }
	ini = a->mem_in_cel();
	if (ini == -1){ fprintf(stderr,"Error: Not enough memory...\n"); return -1; }
	ini = a->mem_in_pv();
	if (ini == -1){ fprintf(stderr,"Error: Not enough memory...\n"); return -1; }
	z = 0;
	for (i = 0; i<nseq; i++)
	{
		for (h = 0; h<height; h++)
		{
			for (j = 0; j<a->nsit[z]; j++)
			{
				a->sta[z][j] = sta[i][j];
				a->cep[z][j] = cep[i][j];
				a->cel[z][j] = cel[i][j];
				a->pv[z][j] = pv[i][j];
			}
			z++;
		}
	}
	return 1;
}
int profile::clear_real(void)
{
	int ini = mem_in_sta();
	if (ini == -1){ fprintf(stderr,"Error: Not enough memory...\n"); return -1; }
	ini = mem_in_cep();
	if (ini == -1){ fprintf(stderr,"Error: Not enough memory...\n"); return -1; }
	ini = mem_in_cel();
	if (ini == -1){ fprintf(stderr,"Error: Not enough memory...\n"); return -1; }
	ini = mem_in_sco();
	if (ini == -1){ fprintf(stderr,"Error: Not enough memory...\n"); return -1; }
	ini = mem_in_pv();
	if (ini == -1){ fprintf(stderr,"Error: Not enough memory...\n"); return -1; }
	return 1;
}
int profile::fprintf_pro(char *mot_db, double thr, char *mode)
{
	//int print_sco=1;
	//if(strncmp(mode,"real",4)==0)print_sco=1;
	//	else print_sco=0;
	int i, j;
	char fileo[ARGLEN];
	FILE *out;
	memset(fileo, '\0', sizeof(fileo));
	strcpy(fileo, mode);
	strcat(fileo, "_");
	strcat(fileo, mot_db);//real or random	
	char buf[10];
	sprintf(buf, "%d", mot);
	strcat(fileo, buf);//hocomoco or dapseq		
	strcat(fileo, "_thr");
	memset(buf, '\0', sizeof(buf));
	sprintf(buf, "%d", nam);
	strcat(fileo, buf);//nomer poroga	
	if ((out = fopen(fileo, "wt")) == NULL)
	{
		fprintf(stderr,"Error: Input file %s can't be opened!\n", fileo);
		return -1;
	}
	for (i = 0; i<nseq; i++)
	{
		fprintf(out, ">Seq %d\tThr %f\tNsites %d\n", i + 1, thr, nsit[i]);
		for (j = 0; j<nsit[i]; j++)
		{
			fprintf(out, "%d\t", sta[i][j]);
			fprintf(out, "%f", pv[i][j]);
			fprintf(out, "\t");
			fprintf(out, "%c\n", cep[i][j]);;
		}
	}
	fclose(out);
	return 1;
}
void profile::count_sites(void)
{
	nsit_all = nseq_rec = 0;
	int i;
	for (i = 0; i<nseq; i++)
	{
		int ns = nsit[i];
		if (ns>0)
		{
			nsit_all += ns;
			nseq_rec++;
		}
	}
}
int profile::test(void)
{
	int i, j;
	for (i = 0; i<nseq; i++)
	{
		for (j = 0; j<nsit[i]; j++)
		{
			if (sta[i][j]>5 * SEQLEN || sta[i][j]<0)return -1;
		}
	}
	return 1;
}
#include "pwm_rec.h" //raspoznavanie matricey

//dlya pods4eta zna4imostey CE, full_overlap, partial_overlap, spacer (overlap = partial_overlap OR full_overlap)
struct count {
	int two_sites;
	int any;
	int partial;
	int full;
	int overlap;
	int spacer;
	void ini(void);
};
void count::ini(void)
{
	two_sites = any = partial = full = overlap = spacer = 0;
}
struct result {
	count cell[NUM_THR][NUM_THR];
	count anc;//toward anchor
	count par;//toward partner
	count asy;//toward anyone
	count equ;//neither
	count sit;
	count anc_sit;
	count par_sit;
	count asy_sit;
	count equ_sit;
	void ini(void);
} observed, expected;

void result::ini(void)
{
	int j, k;
	for (j = 0; j < NUM_THR; j++)for (k = 0; k < NUM_THR; k++)cell[j][k].ini();
	anc.ini();
	par.ini();
	asy.ini();
	equ.ini();
	sit.ini();
	anc_sit.ini();
	par_sit.ini();
	asy_sit.ini();
	equ_sit.ini();
}

double pvalue_a[NUM_THR][NUM_THR];
double pvalue_p[NUM_THR][NUM_THR];
double pvalue_f[NUM_THR][NUM_THR];
double pvalue_s[NUM_THR][NUM_THR];
double pvalue_o[NUM_THR][NUM_THR];

double fold_a[NUM_THR][NUM_THR];
double fold_p[NUM_THR][NUM_THR];
double fold_f[NUM_THR][NUM_THR];
double fold_s[NUM_THR][NUM_THR];
double fold_o[NUM_THR][NUM_THR];

struct pvfo
{
	double p;
	double f;
};
struct pval {
	pvfo anchor;
	pvfo partner;
	pvfo asy1;
	pvfo asy2;
	pvfo equ;
	pvfo anc_par;
	void ini(void);
} pv_any, pv_full, pv_partial, pv_overlap, pv_spacer;
void pval::ini(void)
{
	anchor.f = partner.f = anc_par.f = asy1.f = asy2.f = equ.f = 1;
	anchor.p = partner.p = anc_par.p = asy1.p = asy2.p = equ.p = 1;
}

// dlya vyvoda histogram of CE distribution as fanction of mutual orientation and location of anchor/partner motifs
struct combi {
	//	int n_partial;
	//	int n_full;
	//int n_tot;
	double freq[4][MATLEN + SPACLEN];
	double freqa[MATLEN + SPACLEN];// sum of all orientation
	double freqc[MATLEN + SPACLEN];// sum of all orientation
	//	void ini(int len_a, int len_p, int len_sp);
	//	void mem_out(void);
	int fprintf_all(char *file, int mot, char *motif, int len_a, int len_p, int len_sp, char *mode);
};
int combi::fprintf_all(char *file, int mot, char *motif, int len_a, int len_p, int len_sp, char *mode)
{
	char head[6][12];
	strcpy(head[5], "Cumulative");
	strcpy(head[4], "Any");
	strcpy(head[0], "DirectAP");
	strcpy(head[1], "DirectPA");
	strcpy(head[2], "Inverted");
	strcpy(head[3], "Everted");
	char head0[10];
	strcpy(head0, "DirectAA");
	int n_partial;
	n_partial = Min(len_a, len_p);//partial overlap
	n_partial--;
	int n_full = 1 + abs(len_a - len_p) / 2;// full overlap	
	int n_tot = n_full + n_partial + len_sp + 1;
	char sp = 'S';//spacer
	char bo = 'P';//partial
	char in = 'F';//full
	FILE *out;
	if ((out = fopen(file, mode)) == NULL)
	{
		fprintf(stderr,"Error: Input file %s can't be opened!\n", file);
		return -1;
	}
	fprintf(out, "%d\t%s\t", mot, motif);
	int i, j;
	for (j = n_full; j >= 1; j--)fprintf(out, "\t%d%c", j - 1, in);
	for (j = n_partial; j >= 1; j--)fprintf(out, "\t%d%c", j, bo);
	for (j = 0; j <= len_sp; j++)fprintf(out, "\t%d%c", j, sp);
	fprintf(out, "\n");
	if (mot > 0)
	{
		for (i = 3; i >= 0; i--)
		{
			fprintf(out, "\t\t%s", head[i]);
			for (j = 0; j < n_tot; j++)fprintf(out, "\t%f", 100 * freq[i][j]);
			fprintf(out, "\n");
		}
	}
	else
	{
		for (i = 3; i >= 1; i--)
		{
			if (i != 1)fprintf(out, "\t\t%s", head[i]);
			else fprintf(out, "\t\t%s", head0);
			for (j = 0; j < n_tot; j++)fprintf(out, "\t%f", 100 * freq[i][j]);
			fprintf(out, "\n");
		}
	}
	fprintf(out, "\t\t%s", head[4]);
	for (j = 0; j < n_tot; j++)fprintf(out, "\t%f", 100 * freqa[j]);
	fprintf(out, "\n");
	fprintf(out, "\t\t%s", head[5]);
	for (j = 0; j < n_tot; j++)fprintf(out, "\t%f", 100 * freqc[j]);
	fprintf(out, "\n");
	fclose(out);
	return 1;
}
struct asy_plot {
	int **full;
	int **partial;
	int **overlap;
	int **spacer;
	int **any;
	int total_full;
	int total_partial;
	int total_overlap;
	int total_spacer;
	int total_any;
	double min;
	double max;
	double inter;
	int n_karman;
	int mem_in(int karman, double xmin, double xmax, double xinter);
	void mem_out(void);
	void zero(void);
	void sum(void);
};
int asy_plot::mem_in(int karman, double xmin, double xmax, double xinter)
{
	n_karman = karman - 1;//max karman
	min = xmin, max = xmax, inter = xinter;
	int i;
	full = new int*[karman];
	if (full == NULL) return -1;
	for (i = 0; i<karman; i++)
	{
		full[i] = new int[karman];
		if (full[i] == NULL) return -1;
	}
	partial = new int*[karman];
	if (partial == NULL) return -1;
	for (i = 0; i<karman; i++)
	{
		partial[i] = new int[karman];
		if (partial[i] == NULL) return -1;
	}
	overlap = new int*[karman];
	if (overlap == NULL) return -1;
	for (i = 0; i<karman; i++)
	{
		overlap[i] = new int[karman];
		if (overlap[i] == NULL) return -1;
	}
	spacer = new int*[karman];
	if (spacer == NULL) return -1;
	for (i = 0; i<karman; i++)
	{
		spacer[i] = new int[karman];
		if (spacer[i] == NULL) return -1;
	}
	any = new int*[karman];
	if (any == NULL) return -1;
	for (i = 0; i<karman; i++)
	{
		any[i] = new int[karman];
		if (any[i] == NULL) return -1;
	}
	return 1;
}
void asy_plot::zero(void)
{
	int i, j;
	for (i = 0; i <= n_karman; i++)for (j = 0; j <= n_karman; j++)full[i][j] = partial[i][j] = overlap[i][j] = spacer[i][j] = any[i][j] = 0;
}
void asy_plot::sum(void)
{
	int i, j;
	total_full = total_partial = total_overlap = total_spacer = total_any = 0;
	for (i = 0; i <= n_karman; i++)
	{
		for (j = 0; j <= n_karman; j++)
		{
			total_full += full[i][j];
			total_partial += partial[i][j];
			total_overlap += overlap[i][j];
			total_spacer += spacer[i][j];
			total_any += any[i][j];
		}
	}
}
void asy_plot::mem_out(void)
{
	int i;
	for (i = 0; i <= n_karman; i++) delete[] full[i];
	delete[] full;
	for (i = 0; i <= n_karman; i++) delete[] partial[i];
	delete[] partial;
	for (i = 0; i <= n_karman; i++) delete[] overlap[i];
	delete[] overlap;
	for (i = 0; i <= n_karman; i++) delete[] spacer[i];
	delete[] spacer;
	for (i = 0; i <= n_karman; i++) delete[] any[i];
	delete[] any;
}

#include "projoin.h" // soedinenie dvuh profiley saytov dlya odnogo pika
#include "throw_predictions.h" //permutation of sites in a peak
#include "fisher_exact_test.h" // exact fisher test

#include "pfm_list.h" //spiski imen motivov-partnerov
#include "pfm_similarity.h" //permutation test for anchor/partner motif comparison (separate task for all algorithm)

int main(int argc, char *argv[])
{
	int i, j, k, m, n_motifs, mot, *bad_matrix;
	char file_fasta[ARGLEN], mot_db[30], mypath_data[ARGLEN], prom[ARGLEN], partner_db[30], file_pfm_anchor[ARGLEN];
	char ***seq;// peaks

	char file_hist[ARGLEN], file_hist_rand[ARGLEN], file_pval[5][ARGLEN], file_pval_table[ARGLEN];
	char name_anchor[ARGLEN], name_partner[ARGLEN];
	char xreal[] = "real", xrand[] = "rand", xreal_one[] = "real_one";
	char file_fpr[ARGLEN];
	strcpy(file_fpr, "fpr_anchor.txt");

	if (argc != 8)
	{
		fprintf(stderr,"Error: %s 1file_fasta", argv[0]);//1int thresh_num_min 2int thresh_num_max
		//	printf ("4int height_permut 5int size_min_permut 6int size_max_permut 7double pvalue 8double pvalue_mult");
		fprintf(stderr, " 2char anchor_motif 3char partner_db 4int spacer_min 5int spacer_max 6char path_genome 7double pvalue_thr\n");//9char mot_anchor 
		return -1;
	}
	for (i = 1; i < argc; i++)
	{
		int alen = strlen(argv[i]);
		if (alen > ARGLEN)
		{
			fprintf(stderr, "Error: Argument number %d %s\nof command line is too long!\nMaximim %d symbols allowed\n", i,argv[i],ARGLEN);
			return -1;
		}
	}
	int thresh_num_min = 1, thresh_num_max = 5;	// 1 5    or 5 5
	strcpy(file_fasta, argv[1]);
	int height_permut = 100, size_min_permut = 200000, size_max_permut = 300000; //50000 150000 25  parametry permutacii
	double pvalue_mult = 1.5, dpvalue = 0.0000005; // 0.0005 1.5 parametry dlya porogov matric
	int mot_anchor = 0;// 0 = pwm from file >0 pwm from pre-computed database	
	int s_overlap_min = 6, s_ncycle_small = 1000, s_ncycle_large = 10000;//for permutation(motif_comparison) min_length_of_alignment, no. of permutation (test & detailed)
	double s_granul = 0.001;//for permutation(motif_comparison) okruglenie 4astotnyh matric	
	strcpy(file_pfm_anchor, argv[2]);
	strcpy(partner_db, argv[3]); //hs_core, hs_full, mm_core, mm_full, dapseq
	int shift_min = atoi(argv[4]); // minimal spacer length
	int shift_max = atoi(argv[5]); // maximal spacer length
	strcpy(mypath_data, argv[6]); // ./.../hs, mm, at
	double pvalue = atof(argv[7]);
	double pvalue_max_allowed = 0.002;
	double pvalue_min_allowed = 0.0001;

	if (pvalue > pvalue_max_allowed || pvalue < pvalue_min_allowed)
	{
		printf("Allowed pvalue range [%.3f; %.3f]\n", pvalue_min_allowed, pvalue_max_allowed);
		exit(1);
	}

	strcpy(prom, mypath_data);
	int nseq_genome, len_genome;

	strcpy(file_hist, "out_hist");
	strcpy(file_hist_rand, "out_hist_rand");
	strcpy(file_pval[0], "fisher_any_mot");
	strcpy(file_pval[1], "fisher_full_mot");
	strcpy(file_pval[2], "fisher_part_mot");
	strcpy(file_pval[3], "fisher_over_mot");
	strcpy(file_pval[4], "fisher_spac_mot");
	strcpy(file_pval_table, "out_pval");

	{
		memset(name_anchor, '\0', sizeof(name_anchor));
		int len = strlen(file_pfm_anchor);
		k = 0;
		for (j = 0; j<len; j++)
		{
			char cc = file_pfm_anchor[j];
			if (cc == '.' || cc == '\0')
			{
				name_anchor[k] = '\0';
				break;
			}
			name_anchor[k++] = cc;
		}
		TransStrBack(name_anchor);
	}
	if (strstr(partner_db, "hs_core") != NULL) {
		n_motifs = 403;
		//n_motifs=2;
		bad_matrix = bad_matrix_hs_core;
		nseq_genome = 19795;
		len_genome = 2000;
		strcat(prom, "ups2kb.plain");
		strcpy(mot_db, "hocomoco");
	}
	else
	{
		if (strstr(partner_db, "mm_core") != NULL) {
			n_motifs = 359;
			//n_motifs=2
			bad_matrix = bad_matrix_mm_core;
			nseq_genome = 19991;
			len_genome = 2000;
			strcat(prom, "ups2kb.plain");
			strcpy(mot_db, "hocomoco");
		}
		else
		{
			if (strstr(partner_db, "mm_full") != NULL) {
				n_motifs = 532;
				//n_motifs=2;
				bad_matrix = bad_matrix_mm_full;
				nseq_genome = 19991;
				len_genome = 2000;
				strcat(prom, "ups2kb.plain");
				strcpy(mot_db, "hocomoco");
			}
			else
			{
				if (strstr(partner_db, "hs_full") != NULL) {
					n_motifs = 772;
					//n_motifs=2;
					bad_matrix = bad_matrix_hs_full;
					nseq_genome = 19795;
					len_genome = 2000;
					strcat(prom, "ups2kb.plain");
					strcpy(mot_db, "hocomoco");
				}
				else
				{
					if (strstr(partner_db, "dapseq") != NULL) {
						n_motifs = 529;
						//n_motifs=2;
						bad_matrix = bad_matrix_dapseq;
						nseq_genome = 27193;
						len_genome = 1500;
						strcat(prom, "ups1500.plain");
						strcpy(mot_db, "dapseq");
					}
					else
					{
						fprintf(stderr, "Error: Partner database %s\n", partner_db);
						return -1;
					}
				}
			}
		}
	}
	int bad_motifs = 0;
	for (i = 0; bad_matrix[i] != -1; i++)bad_motifs++;

	int nlen = strlen(name_anchor);//name_anchor[nlen]='\0';
	double pvalue_equal = 0.01;
	double pvalue_similarity_tot;

	int length_fasta_max = 0, nseq_real = 0;
	seq = NULL;
	int ftp = fasta_to_plain0(file_fasta, length_fasta_max, nseq_real);
	if (ftp == -1)
	{
		fprintf(stderr, "Error: Fasta file %s error\n", file_fasta);
		return -1;
	}
	int *peak_len_real;
	peak_len_real = new int[nseq_real];
	if (peak_len_real == NULL){ fprintf(stderr,"Error: Out of memory..."); return -1; }

	seq = new char**[2];
	if (seq == NULL){ fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (k = 0; k<2; k++)
	{
		seq[k] = new char*[nseq_real];
		if (seq[k] == NULL){ fprintf(stderr, "Error: Out of memory..."); return -1; }
		for (i = 0; i<nseq_real; i++)
		{
			int length_fasta_max1 = length_fasta_max + 1;
			seq[k][i] = new char[length_fasta_max1];
			if (seq[k][i] == NULL){ fprintf(stderr, "Error: Out of memory..."); return -1; }
			memset(seq[k][i], '\0', length_fasta_max1);
		}
	}
	ftp = fasta_to_plain1(file_fasta, length_fasta_max, nseq_real, seq, peak_len_real);
	if (ftp == -1)
	{
		fprintf(stderr, "Error: File %s error 2nd stage\n", file_fasta);
		return -1;
	}
	double thr[NUM_THR], thr_anchor[NUM_THR];

	profile real_one[2], rand_one[2], rand_hom_one;
	matrices matrix, matrix0;
	combi hist_obs_one, hist_exp_one;

	//for real
	for (j = 0; j<2; j++)
	{
		real_one[j].nseq = nseq_real;
		real_one[j].nam = 1;
		int ini = real_one[j].mem_in_nsit();
		if (ini == -1){ fprintf(stderr, "Error: Not enough memory...\n"); return -1; }
	}

	//for permutation
	srand((unsigned)time(NULL));
	int nseq_rand = nseq_real*height_permut;
	if (nseq_rand<size_min_permut)height_permut = size_min_permut / nseq_real;
	if (nseq_rand>size_max_permut)height_permut = size_max_permut / nseq_real;
	nseq_rand = nseq_real*height_permut;
	double bonferroni_corr, bonferroni_corr_ap, bonferroni_corr_asy;
	rand_hom_one.nseq = nseq_rand;
	rand_hom_one.nam = 1;
	rand_hom_one.mot = 0;
	int ini = rand_hom_one.mem_in_nsit();
	if (ini == -1){ fprintf(stderr, "Error: Not enough memory...\n"); return -1; }
	for (j = 0; j<2; j++)
	{
		rand_one[j].nseq = nseq_rand;
		rand_one[j].nam = 1;
		int ini = rand_one[j].mem_in_nsit();
		if (ini == -1){ fprintf(stderr, "Error: Not enough memory...\n"); return -1; }
	}

	int *peak_len_rand;
	peak_len_rand = new int[nseq_rand];
	if (peak_len_rand == NULL){ fprintf(stderr, "Error: Out of memory..."); return -1; }
	k = 0;
	for (i = 0; i<nseq_real; i++)
	{
		for (m = 0; m<height_permut; m++)peak_len_rand[k++] = peak_len_real[i];
	}
	int *thr_err_real, *thr_err_rand;
	thr_err_real = new int[nseq_real];
	if (thr_err_real == NULL){ fprintf(stderr, "Error: Out of memory..."); return -1; }
	thr_err_rand = new int[nseq_rand];
	if (thr_err_rand == NULL){ fprintf(stderr, "Error: Out of memory..."); return -1; }

	FILE *out_hist, *out_hist_rand;
	if ((out_hist = fopen(file_hist, "wt")) == NULL)
	{
		fprintf(stderr, "Error: Input file %s can't be opened!\n", file_hist);
		return -1;
	}
	fclose(out_hist);
	if ((out_hist_rand = fopen(file_hist_rand, "wt")) == NULL)
	{
		fprintf(stderr, "Error: Input file %s can't be opened!\n", file_hist_rand);
		return -1;
	}
	fclose(out_hist_rand);
	FILE *out_pval_table;
	if ((out_pval_table = fopen(file_pval_table, "wt")) == NULL)
	{
		fprintf(stderr, "Error: Input file %s can't be opened!\n", file_pval_table);
		return -1;
	}
	fprintf(out_pval_table, "# Motif");
	fprintf(out_pval_table, "\tMotif Name");
	fprintf(out_pval_table, "\tFull overlap, -Log10[P-value]");
	fprintf(out_pval_table, "\tPartial overlap,-Log10[P-value]");
	fprintf(out_pval_table, "\tOverlap, -Log10[P-value]");
	fprintf(out_pval_table, "\tSpacer, -Log10[P-value]");
	fprintf(out_pval_table, "\tAny, -Log10[P-value]");
	fprintf(out_pval_table, "\tSimilarity to Anchor, -Log10[P-value]");
	fprintf(out_pval_table, "\tSimilarity to Anchor, SSD");
	fprintf(out_pval_table, "\tSimilarity to Anchor, PCC\t");
	//	fprintf(out_pval_table, "Full overlap, Conservative Anchor, -Log10[P-value]\t");
	//	fprintf(out_pval_table, "Partial overlap, Conservative Anchor, -Log10[P-value]\t");
	fprintf(out_pval_table, "Overlap, Conservative Anchor, -Log10[P-value]\t");
	fprintf(out_pval_table, "Spacer, Conservative Anchor, -Log10[P-value]\t");
	//fprintf(out_pval_table, "Any, Conservative Anchor, -Log10[P-value]\t");
//	fprintf(out_pval_table, "Full overlap, Conservative Partner, -Log10[P-value]\t");
//	fprintf(out_pval_table, "Partial overlap, Conservative Partner, -Log10[P-value]\t");
	fprintf(out_pval_table, "Overlap, Conservative Partner, -Log10[P-value]\t");
	fprintf(out_pval_table, "Spacer, Conservative Partner, -Log10[P-value]\t");
	//	fprintf(out_pval_table, "Any, Conservative Partner, -Log10[P-value]\t");
	fprintf(out_pval_table, "Overlap, Conservative OneMotif, -Log10[P-value]\t");
	fprintf(out_pval_table, "Spacer, Conservative OneMotif, -Log10[P-value]\t");
	fprintf(out_pval_table, "Overlap, Equal Conservation, -Log10[P-value]\t");
	fprintf(out_pval_table, "Spacer, Equal Conservation, -Log10[P-value]\t");
	//	fprintf(out_pval_table, "Any, Conservative OneMotif, -Log10[P-value]\t");
	//	fprintf(out_pval_table, "Full overlap, Asymmetry to Anchor+/Partner-, -Log10[P-value]\t");
		//fprintf(out_pval_table, "Partial overlap, Asymmetry to Anchor+/Partner-, -Log10[P-value]\t");
	fprintf(out_pval_table, "Overlap, Sites Anchor+/Partner-, -Log10[P-value]\t");
	fprintf(out_pval_table, "Spacer, Sites Anchor+/Partner-, -Log10[P-value]\t");
	fprintf(out_pval_table, "Overlap, Sites Asymmetry/Symmetry, -Log10[P-value]\t");
	fprintf(out_pval_table, "Spacer, Sites Asymmetry/Symmetry, -Log10[P-value]\t");
	//	fprintf(out_pval_table, "Any, Asymmetry to Anchor+/Partner-, -Log10[P-value]\t");
	fprintf(out_pval_table, "Bonferroni_CE\tBonferroni_CE(AncPar)\tBonferroni_Asym\n");
	fclose(out_pval_table);

	FILE *out_pval[5];
	double pval_sim[4];
	double dthr_asy = 0.2, thr_asy_min = 0.1 * ((int)(10 * (dthr_asy - log10(pvalue)))), thr_asy_max = thr_asy_min + 2;
	int nthr_asy = (int)((thr_asy_max - thr_asy_min) / dthr_asy) + 2;
	asy_plot real_plot, rand_plot;
	real_plot.mem_in(nthr_asy, thr_asy_min, thr_asy_max, dthr_asy);
	rand_plot.mem_in(nthr_asy, thr_asy_min, thr_asy_max, dthr_asy);

	FILE *out_stat;
	if ((out_stat = fopen("rec_pos.txt", "wt")) == NULL)
	{
		fprintf(stderr, "Error: Input file can't be opened!\n");
		return -1;
	}
	fprintf(out_stat, "# Motif\tMotif Name\t# Threshold\tThreshold\t%% of peaks\tRec. peaks\tTotal peaks\tRate of hits\tRec. hits\tTotal positions\n");
	char file_err[] = "throw_prediction.txt";
	{
		FILE *out_err;
		if ((out_err = fopen(file_err, "wt")) == NULL)
		{
			fprintf(stderr, "Input file %s can't be opened!\n", file_err);
			return -1;
		}
	}
	char file_log[] = "mcot.log";
	{
		FILE *out_log;
		if ((out_log = fopen(file_log, "wt")) == NULL)
		{
			fprintf(out_log, "Error: Input file %s can't be opened!\n", file_log);
			return -1;
		}
		fclose(out_log);
	}
	int len_anchor = 0, len_partner = 0;
	for (mot = 0; mot<n_motifs; mot++)
		//for(mot=0;mot<n_motifs;mot+=144)
	{
		{
			int bad = 0;
			for (i = 0; bad_matrix[i] != -1; i++)
			{
				if (mot == bad_matrix[i])
				{
					bad = 1;
					break;
				}
			}
			if (bad == 1)continue;
		}
		printf("Mot %d\n", mot);		
		int index[NUM_THR];
		double *thr_all;
		int n_thr_all = 0;
		double *fp_rate;
		if (mot == 0)
		{
			double pwm_anchor[MATLEN][OLIGNUM];
			int length = pfm_to_pwm(file_pfm_anchor, &matrix);
			if (length <= 0 || length >= MATLEN)
			{
				fprintf(stderr, "Error: PFM to PWM conversion error, file %s\n", file_pfm_anchor);
				return -1;
			}
			for (i = 0; i<length; i++)for (j = 0; j<OLIGNUM; j++)pwm_anchor[i][j] = matrix.wei[i][j];
			matrix.get_copy(&matrix0);
			int all_pos_rec = 1+(int)(2* pvalue*nseq_genome*(len_genome - matrix.len + 1));
			thr_all = new double[all_pos_rec];
			if (thr_all == NULL) { fprintf(stderr,"Error: Out of memory..."); return -1; }
			fp_rate = new double[all_pos_rec];
			if (fp_rate == NULL){ fprintf(stderr,"Error: Out of memory..."); return -1; }
			for (i = 0; i < all_pos_rec; i++)thr_all[i] = 0;
			int nthr_dist = 0;
			int piptd = pwm_iz_pwm_thr_dist0(pwm_anchor, length, prom, all_pos_rec, nthr_dist, thr_all, fp_rate, mypath_data, nseq_genome, len_genome, pvalue,dpvalue);
			if (piptd == -1)
			{
				fprintf(stderr, "Error: FP rate table error\n");
				return -1;
			}			
			FILE *out_fpr;
			if ((out_fpr = fopen(file_fpr, "wt")) == NULL)
			{
				fprintf(stderr, "Error: Output file %s can't be opened!\n", file_fpr);
				return -1;
			}
			for (i = 0; i < nthr_dist; i++)fprintf(out_fpr, "%.8f\t%g\n", thr_all[i], fp_rate[i]);
			fclose(out_fpr);
			double fpr_select[NUM_THR];
			int stfp = select_thresholds_from_pvalues(nthr_dist, thr_all, fp_rate, pvalue, pvalue_mult, fpr_select, thr, index);
			if (stfp == -1)
			{
				fprintf(stderr, "Error: Bad input matrix of %d motif\n", mot);
				return -1;
			}
			for (i = 0; i < nthr_dist; i++)
			{
				if (fp_rate[i]>0)fp_rate[i] = -log10(fp_rate[i]);
				else fp_rate[i] = 10;
			}
			//			delete [] fp_rate;							
			memset(name_partner, '\0', sizeof(name_partner));
			strcpy(name_partner, name_anchor);
			int nlen = strlen(name_anchor);
			name_partner[nlen] = '\0';
			pvalue_similarity_tot = 1E-300;
			for (i = 0; i<4; i++)pval_sim[i] = pvalue_similarity_tot;
		}
		else
		{
			if (strstr(partner_db, "hs_core") != NULL)
			{
				int mem = motifs.hs_core_ini(hs_core_thr_count, mot);
				if (mem == -1){ fprintf(stderr,"Error: Not enough memory...");	exit(1); }
			}
			else
			{
				if (strstr(partner_db, "mm_core") != NULL)
				{
					int mem = motifs.mm_core_ini(mm_core_thr_count, mot);
					if (mem == -1){ fprintf(stderr,"Error: Not enough memory...");	exit(1); }
				}
				else
				{
					if (strstr(partner_db, "hs_full") != NULL)
					{
						int mem = motifs.hs_full_ini(hs_full_thr_count, mot);
						if (mem == -1){ fprintf(stderr,"Error: Not enough memory...");	exit(1); }
					}
					else
					{
						if (strstr(partner_db, "dapseq") != NULL)
						{
							int mem = motifs.dapseq_ini(dapseq_thr_count, mot);
							if (mem == -1){ fprintf(stderr,"Error: Not enough memory...");	exit(1); }
						}
					}
				}
			}
			for (j = 0; j<NUM_THR; j++)
			{
				thr[j] = motifs.thr_selected[j];
				index[j] = motifs.inx_selected[j];
			}
			n_thr_all = motifs.nthr;
			fp_rate = new double[n_thr_all];
			if (fp_rate == NULL){ fprintf(stderr,"Error: Out of memory..."); return -1; }
			thr_all = new double[n_thr_all];
			if (thr_all == NULL){ fprintf(stderr,"Error: Out of memory..."); return -1; }
			for (j = 0; j<n_thr_all; j++)
			{
				fp_rate[j] = motifs.fpr[j];
				thr_all[j] = motifs.thr[j];
			}
			motifs.mem_out();
			memset(name_partner, '\0', sizeof(name_partner));
			if (strcmp(partner_db, "hs_core") == 0)
			{
				matrix.init_hs_core(mot);
				strcpy(name_partner, hs_core_names[mot]);
			}
			else
			{
				if (strcmp(partner_db, "mm_core") == 0)
				{
					matrix.init_mm_core(mot);
					strcpy(name_partner, mm_core_names[mot]);
				}
				else
				{
					if (strcmp(partner_db, "hs_full") == 0)
					{
						matrix.init_hs_full(mot);
						strcpy(name_partner, hs_full_names[mot]);
					}
					else
					{
						if (strcmp(partner_db, "mm_full") == 0)
						{
							matrix.init_mm_full(mot);
							strcpy(name_partner, mm_full_names[mot]);
						}
						else
						{
							if (strcmp(partner_db, "dapseq") == 0)
							{
								matrix.init_dapseq(mot);
								strcpy(name_partner, dapseq_names[mot]);
							}
							else
							{
								fprintf(stderr, "Error: Partner database %s\n", partner_db);
								return -1;
							}
						}
					}
				}
			}
			int nlen = strlen(name_partner);
			name_partner[nlen] = '\0';
			for (i = 0; i<4; i++)pval_sim[i] = 1;
			pvalue_similarity_tot = pfm_similarity(&matrix, &matrix0, s_granul, s_overlap_min, s_ncycle_small, s_ncycle_large, pval_sim);
		}
		matrix.norm();
		int ap;
		if (mot == 0)
		{
			ap = 0;
			for (j = 0; j<NUM_THR; j++)thr_anchor[j] = thr[j];
		}
		else
		{
			ap = 1;
		}
		//int pwm_rec0(matrices *mat, double thr, int len_pro, int nseq_pro, char ***seq, profile *real)  count all sites
		////recognition 1st		
		int all_pos = 0;//total number of available positions
		int wm_rec = pwm_rec0(&matrix, thr[NUM_THR - 1], length_fasta_max, nseq_real, seq, &real_one[ap], all_pos);
		if (wm_rec == -1)
		{
			fprintf(stderr, "Error: Motif %d recognition 1st stage error\n", mot);
			return -1;
		}
		//memory allocation for all sites
		real_one[ap].clear_real();
		real_one[ap].mot = mot;
		real_one[ap].nam = NUM_THR;
		real_one[ap].count_sites();
		//recognition 2nd
		wm_rec = pwm_rec1(&matrix, thr[NUM_THR - 1], length_fasta_max, nseq_real, seq, &real_one[ap]);
		if (wm_rec == -1)
		{
			fprintf(stderr, "Error: Motif %d recognition 2nd stage error\n", mot);
			return -1;
		}
		//count nsites for various thresholds
		for (i = 0; i<nseq_real; i++)
		{
			for (k = 0; k<real_one[ap].nsit[i]; k++)
			{
				double sco = real_one[ap].sco[i][k];
				for (j = 0; j<NUM_THR; j++)
				{
					if (sco >= thr[j])
					{
						real_one[ap].cel[i][k] = j;
						break;
					}
				}
			}
		}
		//initiation of profiles for various thresholds
		int rec_seq[NUM_THR], rec_pos[NUM_THR];
		for (i = 0; i<NUM_THR; i++)rec_seq[i] = rec_pos[i] = 0;
		for (i = 0; i<nseq_real; i++)
		{
			int inx[NUM_THR];
			for (j = 0; j<NUM_THR; j++)inx[j] = 0;
			for (k = 0; k<real_one[ap].nsit[i]; k++)
			{
				double sco = real_one[ap].sco[i][k];
				for (j = 0; j<NUM_THR; j++)
				{
					if (sco >= thr[j])
					{
						rec_pos[j]++;
						inx[j]++;
						if (inx[j] == 1)rec_seq[j]++;
						break;
					}
				}
			}
		}
		//transform score to fprate
		for (i = 0; i<nseq_real; i++)
		{
			for (k = 0; k<real_one[ap].nsit[i]; k++)
			{
				double sco = real_one[ap].sco[i][k];
				int inter = real_one[ap].cel[i][k];
				int inx2, inx1 = index[inter];
				if (inter == 0)inx2 = 0;
				else inx2 = index[inter - 1];
				double pv_sc = 0;
				for (j = inx1; j >= inx2; j--)
				{
					if (thr_all[j + 1] <= sco && thr_all[j] >= sco)
					{
						pv_sc = fp_rate[j];
						break;
					}
				}
				real_one[ap].pv[i][k] = pv_sc;
			}
		}
		delete[] fp_rate;
		delete[] thr_all;
		int fprint_pro = real_one[ap].fprintf_pro(mot_db, thr[NUM_THR - 1], xreal);
		if (fprint_pro == -1)
		{
			fprintf(stderr, "Error: Real print profile error, motif %d\n", mot);
			return -1;
		}
		for (j = 0; j<NUM_THR; j++)
		{
			if (mot == 0)fprintf(out_stat, "Anchor");
			else fprintf(out_stat, "Partner %d", mot);
			fprintf(out_stat, "\t%s\t%d\t%f\t", name_partner, j + 1, thr[j]);
			fprintf(out_stat, "%f\t%d\t%d\t", 100 * (double)rec_seq[j] / nseq_real, rec_seq[j], nseq_real);
			fprintf(out_stat, "%g\t%d\t%d\n", (double)rec_pos[j] / all_pos, rec_pos[j], all_pos);
		}
		if (ap == 0)
		{ 
			len_anchor = len_partner = matrix.len; 
		}
		else len_partner = matrix.len;
		matrix.mem_out(matrix.len);
		//int len_ap=matrix.len;					
		/*int test;
		test = real_one[ap].test();
		if(test==-1)
		{
		printf("Real One %d error Mot %d Thr One\n",ap, mot);
		return -1;
		}*/

		//one thresh  rand - &rand_hom_one,&rand_one[ap]   real - real_one[0],real_one[ap]
		//anchor
		if (mot == 0)
		{
			int cop = real_one[0].get_copy_rand(&rand_hom_one, height_permut);
			if (cop == -1)
			{
				fprintf(stderr, "Error: Rand Copy error Mot %d Thr One\n", mot);
				return -1;
			}
			/*		int fprint_pro=rand_hom_one.fprintf_pro(mot_db,thr[NUM_THR-1],xrand);
			if(fprint_pro==-1)
			{
			printf("Rand print profile error, motif %d\n",mot);
			return -1;
			}*/

			/*
			test = rand_hom_one.test();
			if(test==-1)
			{
			printf("Rand Hom One error Mot %d Thr One\n", mot);
			return -1;
			}
			for(j=0;j<NUM_THR;j++)
			{
			test = rand_hom[j].test();
			if(test==-1)
			{
			printf("Rand Hom error Mot %d Thr %d\n",mot,j+1);
			return -1;
			}
			}*/

		}
		//partner
		{
			int cop = real_one[ap].get_copy_rand(&rand_one[ap], height_permut);
			if (cop == -1)
			{
				fprintf(stderr, "Error: Rand Copy error Mot %d Thr One\n", mot);
				return -1;
			}
			/*
			test = rand_one[ap].test();
			if(test==-1)
			{
			printf("Rand One %d error Mot %d Thr One\n",ap, mot);
			return -1;
			}
			for(j=0;j<NUM_THR;j++)
			{
			test = rand[ap][j].test();
			if(test==-1)
			{
			printf("Rand %d error Mot %d Thr %d\n",ap, mot,j+1);
			return -1;
			}
			}*/
		}
		{// one threshold
			for (m = 0; m<nseq_real; m++)thr_err_real[m] = 0;
			rand_one[ap].mot = mot;
			int throwp = throw_predictions(peak_len_rand, &rand_hom_one, &rand_one[ap], len_anchor, len_partner, 0, thr_err_real, nseq_real, nseq_rand, seq[0], height_permut, file_err);
			if (throwp == -1)
			{
				fprintf(stderr, "Error: Throw Prediction error One - Anc 0 Par %d\n", mot);
				return -1;
			}
			real_plot.zero();
			rand_plot.zero();
			observed.ini();
			expected.ini();
			k = 0;
			for (i = 0; i<nseq_real; i++)
			{
				for (j = 0; j<height_permut; j++)
				{
					thr_err_rand[k++] = thr_err_real[i];
				}
			}
			//printf("Mot %d Enter projoin\n",mot);
			int nseq_two_sites_real = 0, nseq_two_sites_rand = 0;
			int proj = projoin(xrand, mot_db, rand_hom_one, rand_one[ap], shift_min, shift_max, len_anchor, len_partner, thr_err_rand, nseq_rand, seq, &expected, &hist_exp_one, peak_len_rand, &rand_plot, nseq_two_sites_rand);
			if (proj == -1)
			{
				fprintf(stderr, "Error: Projoin Rand error Anc 0 Par %d\n", mot);
				return -1;
			}
			proj = projoin(xreal, mot_db, real_one[0], real_one[ap], shift_min, shift_max, len_anchor, len_partner, thr_err_real, nseq_real, seq, &observed, &hist_obs_one, peak_len_real, &real_plot, nseq_two_sites_real);
			if (proj == -1)
			{
				fprintf(stderr, "Error: Projoin Real error Anc 0 Par %d\n", mot);
				return -1;
			}
			{
				bonferroni_corr = (double)nseq_two_sites_real*nseq_two_sites_rand;
				bonferroni_corr *= 5;//potoki
				bonferroni_corr *= (n_motifs - bad_motifs - 1);
				bonferroni_corr_ap = bonferroni_corr;
				bonferroni_corr_asy = bonferroni_corr;
				bonferroni_corr *= NUM_THR;
				bonferroni_corr *= NUM_THR;
				bonferroni_corr_ap *= 2;
				double pv_standard = -log10(0.05);
				bonferroni_corr = pv_standard + log10(bonferroni_corr);
				bonferroni_corr_ap = pv_standard + log10(bonferroni_corr_ap);
				bonferroni_corr_asy = pv_standard + log10(bonferroni_corr_asy);
			}

			//printf("Mot %d Enter hist\n",mot);
			char modew[] = "wt", modea[] = "at";
			char file_hist_one[ARGLEN];
			strcpy(file_hist_one, file_hist);
			char buf[4];
			memset(buf, '\0', sizeof(buf));
			sprintf(buf, "%d", mot);
			strcat(file_hist_one, buf);
			hist_obs_one.fprintf_all(file_hist, mot, name_partner, len_anchor, len_partner, shift_max, modea);
			hist_obs_one.fprintf_all(file_hist_one, mot, name_partner, len_anchor, len_partner, shift_max, modew);
			hist_exp_one.fprintf_all(file_hist_rand, mot, name_partner, len_anchor, len_partner, shift_max, modea);
			real_plot.sum();
			rand_plot.sum();
			//printf("Mot %d Enter plot\n",mot);
			//if (mot != 0)
			{
				char flow[5][8] = { "Full", "Partial", "Overlap", "Spacer", "Any" };
				FILE *out_plot[5];
				char file_plot[5][ARGLEN];				
				char rzd = ',';				
				if (real_plot.total_full != 0 && rand_plot.total_full != 0)
				{
					strcpy(file_plot[0], "plot_");
					strcat(file_plot[0], flow[0]);
					strcat(file_plot[0], "_Anchor");
					if (mot != 0)strcat(file_plot[0], "_Partner");
					else strcat(file_plot[0], "_Anchor");
					strcat(file_plot[0], buf);
					if ((out_plot[0] = fopen(file_plot[0], "wt")) == NULL)
					{
						fprintf(stderr, "Input file %s can't be opened!\n", file_plot[0]);
						return -1;
					}
					double val = real_plot.min;
					fprintf(out_plot[0], "%c<%.1f", rzd, val);
					for (i = 1; i < real_plot.n_karman; i++)
					{
						fprintf(out_plot[0], "%c%.1f", rzd, val);
						val += real_plot.inter;
						fprintf(out_plot[0], "..%.1f", val);
					}
					fprintf(out_plot[0], "%c>%.1f", rzd, val);
					fprintf(out_plot[0], "\n");
					val = real_plot.min;
					for (i = 0; i <= real_plot.n_karman; i++)
					{
						if (i == 0)fprintf(out_plot[0], "<%.1f", val);
						else
						{
							if (i == real_plot.n_karman)fprintf(out_plot[0], ">%.1f", val);
							else
							{
								fprintf(out_plot[0], "%.1f", val);
								val += real_plot.inter;
								fprintf(out_plot[0], "..%.1f", val);
							}
						}
						for (j = 0; j <= real_plot.n_karman; j++)
						{
							fprintf(out_plot[0], "%c", rzd);
							{
								double mill = 1000 * ((double)real_plot.full[i][j] / real_plot.total_full - (double)rand_plot.full[i][j] / rand_plot.total_full);
								if (fabs(mill) > 0.5)fprintf(out_plot[0], "%.f", mill);
							}
						}
						fprintf(out_plot[0], "\n");
					}
					fclose(out_plot[0]);
				}
				if (real_plot.total_partial != 0 && rand_plot.total_partial != 0)
				{
					strcpy(file_plot[1], "plot_");
					strcat(file_plot[1], flow[1]);
					strcat(file_plot[1], "_Anchor");
					if (mot != 0)strcat(file_plot[1], "_Partner");
					else strcat(file_plot[1], "_Anchor");
					strcat(file_plot[1], buf);
					if ((out_plot[1] = fopen(file_plot[1], "wt")) == NULL)
					{
						fprintf(stderr, "Input file %s can't be opened!\n", file_plot[1]);
						return -1;
					}
					double val = real_plot.min;
					fprintf(out_plot[1], "%c<%.1f", rzd, val);
					for (i = 1; i < real_plot.n_karman; i++)
					{
						fprintf(out_plot[1], "%c%.1f", rzd, val);
						val += real_plot.inter;
						fprintf(out_plot[1], "..%.1f", val);
					}
					fprintf(out_plot[1], "%c>%.1f", rzd, val);
					fprintf(out_plot[1], "\n");
					val = real_plot.min;
					for (i = 0; i <= real_plot.n_karman; i++)
					{
						if (i == 0)fprintf(out_plot[1], "<%.1f", val);
						else
						{
							if (i == real_plot.n_karman)fprintf(out_plot[1], ">%.1f", val);
							else
							{
								fprintf(out_plot[1], "%.1f", val);
								val += real_plot.inter;
								fprintf(out_plot[1], "..%.1f", val);
							}
						}
						for (j = 0; j <= real_plot.n_karman; j++)
						{
							fprintf(out_plot[1], "%c", rzd);
							{
								double mill = 1000 * ((double)real_plot.partial[i][j] / real_plot.total_partial - (double)rand_plot.partial[i][j] / rand_plot.total_partial);
								if (fabs(mill) > 0.5)fprintf(out_plot[1], "%.f", mill);
							}
						}
						fprintf(out_plot[1], "\n");
					}
					fclose(out_plot[1]);
				}
				if (real_plot.total_overlap != 0 && rand_plot.total_overlap != 0)
				{
					strcpy(file_plot[2], "plot_");
					strcat(file_plot[2], flow[2]);
					strcat(file_plot[2], "_Anchor");
					if (mot != 0)strcat(file_plot[2], "_Partner");
					else strcat(file_plot[2], "_Anchor");
					strcat(file_plot[2], buf);
					if ((out_plot[2] = fopen(file_plot[2], "wt")) == NULL)
					{
						fprintf(stderr, "Error: Input file %s can't be opened!\n", file_plot[2]);
						return -1;
					}
					double val = real_plot.min;
					fprintf(out_plot[2], "%c<%.1f", rzd, val);
					for (i = 1; i < real_plot.n_karman; i++)
					{
						fprintf(out_plot[2], "%c%.1f", rzd, val);
						val += real_plot.inter;
						fprintf(out_plot[2], "..%.1f", val);
					}
					fprintf(out_plot[2], "%c>%.1f", rzd, val);
					fprintf(out_plot[2], "\n");
					val = real_plot.min;
					for (i = 0; i <= real_plot.n_karman; i++)
					{
						if (i == 0)fprintf(out_plot[2], "<%.1f", val);
						else
						{
							if (i == real_plot.n_karman)fprintf(out_plot[2], ">%.1f", val);
							else
							{
								fprintf(out_plot[2], "%.1f", val);
								val += real_plot.inter;
								fprintf(out_plot[2], "..%.1f", val);
							}
						}
						for (j = 0; j <= real_plot.n_karman; j++)
						{
							fprintf(out_plot[2], "%c", rzd);
							{
								double mill = 1000 * ((double)real_plot.overlap[i][j] / real_plot.total_overlap - (double)rand_plot.overlap[i][j] / rand_plot.total_overlap);
								if (fabs(mill) > 0.5)fprintf(out_plot[2], "%.f", mill);
							}
						}
						fprintf(out_plot[2], "\n");
					}
					fclose(out_plot[2]);
				}
				if (real_plot.total_spacer != 0 && rand_plot.total_spacer != 0)
				{
					strcpy(file_plot[3], "plot_");
					strcat(file_plot[3], flow[3]);
					strcat(file_plot[3], "_Anchor");
					if (mot != 0)strcat(file_plot[3], "_Partner");
					else strcat(file_plot[3], "_Anchor");
					strcat(file_plot[3], buf);
					if ((out_plot[3] = fopen(file_plot[3], "wt")) == NULL)
					{
						fprintf(stderr, "Error: Input file %s can't be opened!\n", file_plot[3]);
						return -1;
					}
					double val = real_plot.min;
					fprintf(out_plot[3], "%c<%.1f", rzd, val);
					for (i = 1; i < real_plot.n_karman; i++)
					{
						fprintf(out_plot[3], "%c%.1f", rzd, val);
						val += real_plot.inter;
						fprintf(out_plot[3], "..%.1f", val);
					}
					fprintf(out_plot[3], "%c>%.1f", rzd, val);
					fprintf(out_plot[3], "\n");
					val = real_plot.min;
					for (i = 0; i <= real_plot.n_karman; i++)
					{
						if (i == 0)fprintf(out_plot[3], "<%.1f", val);
						else
						{
							if (i == real_plot.n_karman)fprintf(out_plot[3], ">%.1f", val);
							else
							{
								fprintf(out_plot[3], "%.1f", val);
								val += real_plot.inter;
								fprintf(out_plot[3], "..%.1f", val);
							}
						}
						for (j = 0; j <= real_plot.n_karman; j++)
						{
							fprintf(out_plot[3], "%c", rzd);
							{
								double mill = 1000 * ((double)real_plot.spacer[i][j] / real_plot.total_spacer - (double)rand_plot.spacer[i][j] / rand_plot.total_spacer);
								if (fabs(mill) > 0.5)fprintf(out_plot[3], "%.f", mill);
							}
						}
						fprintf(out_plot[3], "\n");
					}
					fclose(out_plot[3]);
				}
				if (real_plot.total_any != 0 && rand_plot.total_any != 0)
				{
					strcpy(file_plot[4], "plot_");
					strcat(file_plot[4], flow[4]);
					strcat(file_plot[4], "_Anchor");
					if (mot != 0)strcat(file_plot[4], "_Partner");
					else strcat(file_plot[4], "_Anchor");
					strcat(file_plot[4], buf);
					if ((out_plot[4] = fopen(file_plot[4], "wt")) == NULL)
					{
						fprintf(stderr, "Error: Input file %s can't be opened!\n", file_plot[4]);
						return -1;
					}
					double val = real_plot.min;
					fprintf(out_plot[4], "%c<%.1f", rzd, val);
					for (i = 1; i < real_plot.n_karman; i++)
					{
						fprintf(out_plot[4], "%c%.1f", rzd, val);
						val += real_plot.inter;
						fprintf(out_plot[4], "..%.1f", val);
					}
					fprintf(out_plot[4], "%c>%.1f", rzd, val);
					fprintf(out_plot[4], "\n");
					val = real_plot.min;
					for (i = 0; i <= real_plot.n_karman; i++)
					{
						if (i == 0)fprintf(out_plot[4], "<%.1f", val);
						else
						{
							if (i == real_plot.n_karman)fprintf(out_plot[4], ">%.1f", val);
							else
							{
								fprintf(out_plot[4], "%.1f", val);
								val += real_plot.inter;
								fprintf(out_plot[4], "..%.1f", val);
							}
						}
						for (j = 0; j <= real_plot.n_karman; j++)
						{
							fprintf(out_plot[4], "%c", rzd);
							{
								double mill = 1000 * ((double)real_plot.any[i][j] / real_plot.total_any - (double)rand_plot.any[i][j] / rand_plot.total_any);
								if (fabs(mill) > 0.5)fprintf(out_plot[4], "%.f", mill);
							}
						}
						fprintf(out_plot[4], "\n");
					}
					fclose(out_plot[4]);
				}
			}
		}// one threshold
		//many thresholds
		for (i = 0; i<NUM_THR; i++)for (j = 0; j<NUM_THR; j++)pvalue_a[i][j] = pvalue_f[i][j] = pvalue_p[i][j] = pvalue_o[i][j] = pvalue_s[i][j] = 1;
		for (j = 0; j<NUM_THR; j++)
		{
			for (k = 0; k<NUM_THR; k++)
			{
				//printf("Mot %d J %d K %d\n",mot,j+1,k+1);								
				int fisher = fisher_exact_test(observed.cell[j][k].any, observed.cell[j][k].two_sites, expected.cell[j][k].any, expected.cell[j][k].two_sites, pvalue_a[j][k], fold_a[j][k], 0);
				if (fisher == -1)
				{
					fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher = fisher_exact_test(observed.cell[j][k].full, observed.cell[j][k].two_sites, expected.cell[j][k].full, expected.cell[j][k].two_sites, pvalue_f[j][k], fold_f[j][k], 0);
				if (fisher == -1)
				{
					fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher = fisher_exact_test(observed.cell[j][k].partial, observed.cell[j][k].two_sites, expected.cell[j][k].partial, expected.cell[j][k].two_sites, pvalue_p[j][k], fold_p[j][k], 0);
				if (fisher == -1)
				{
					fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher = fisher_exact_test(observed.cell[j][k].overlap, observed.cell[j][k].two_sites, expected.cell[j][k].overlap, expected.cell[j][k].two_sites, pvalue_o[j][k], fold_o[j][k], 0);
				if (fisher == -1)
				{
					fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
				fisher = fisher_exact_test(observed.cell[j][k].spacer, observed.cell[j][k].two_sites, expected.cell[j][k].spacer, expected.cell[j][k].two_sites, pvalue_s[j][k], fold_s[j][k], 0);
				if (fisher == -1)
				{
					fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
					return -1;
				}
			}
		}
		{
			pv_any.ini();
			pv_full.ini();
			pv_partial.ini();
			pv_overlap.ini();
			pv_spacer.ini();
			//any
			int fisher;
			fisher = fisher_exact_test(observed.anc.any, observed.anc.two_sites, expected.anc.any, expected.anc.two_sites, pv_any.anchor.p, pv_any.anchor.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.par.any, observed.par.two_sites, expected.par.any, expected.par.two_sites, pv_any.partner.p, pv_any.partner.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.asy.any, observed.asy.two_sites, expected.asy.any, expected.asy.two_sites, pv_any.asy1.p, pv_any.asy1.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.equ.any, observed.equ.two_sites, expected.equ.any, expected.equ.two_sites, pv_any.equ.p, pv_any.equ.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.anc_sit.any, observed.sit.any, expected.anc_sit.any, expected.sit.any, pv_any.anc_par.p, pv_any.anc_par.f, 1);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.asy_sit.any, observed.sit.any, expected.asy_sit.any, expected.sit.any, pv_any.asy2.p, pv_any.asy2.f, 1);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			//full
			fisher = fisher_exact_test(observed.anc.full, observed.anc.two_sites, expected.anc.full, expected.anc.two_sites, pv_full.anchor.p, pv_full.anchor.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.par.full, observed.par.two_sites, expected.par.full, expected.par.two_sites, pv_full.partner.p, pv_full.partner.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.asy.full, observed.asy.two_sites, expected.asy.full, expected.asy.two_sites, pv_full.asy1.p, pv_full.asy1.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.equ.full, observed.equ.two_sites, expected.equ.full, expected.equ.two_sites, pv_full.equ.p, pv_full.equ.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.anc_sit.full, observed.sit.full, expected.anc_sit.full, expected.sit.full, pv_full.anc_par.p, pv_full.anc_par.f, 1);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.asy_sit.full, observed.sit.full, expected.asy_sit.full, expected.sit.full, pv_full.asy2.p, pv_full.asy2.f, 1);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			//partial
			fisher = fisher_exact_test(observed.anc.partial, observed.anc.two_sites, expected.anc.partial, expected.anc.two_sites, pv_partial.anchor.p, pv_partial.anchor.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.par.partial, observed.par.two_sites, expected.par.partial, expected.par.two_sites, pv_partial.partner.p, pv_partial.partner.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.asy.partial, observed.asy.two_sites, expected.asy.partial, expected.asy.two_sites, pv_partial.asy1.p, pv_partial.asy1.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.equ.partial, observed.equ.two_sites, expected.equ.partial, expected.equ.two_sites, pv_partial.equ.p, pv_partial.equ.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.anc_sit.partial, observed.sit.partial, expected.anc_sit.partial, expected.sit.partial, pv_partial.anc_par.p, pv_partial.anc_par.f, 1);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.asy_sit.partial, observed.sit.partial, expected.asy_sit.partial, expected.sit.partial, pv_partial.asy2.p, pv_partial.asy2.f, 1);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			//overlap
			fisher = fisher_exact_test(observed.anc.overlap, observed.anc.two_sites, expected.anc.overlap, expected.anc.two_sites, pv_overlap.anchor.p, pv_overlap.anchor.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.par.overlap, observed.par.two_sites, expected.par.overlap, expected.par.two_sites, pv_overlap.partner.p, pv_overlap.partner.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.asy.overlap, observed.asy.two_sites, expected.asy.overlap, expected.asy.two_sites, pv_overlap.asy1.p, pv_overlap.asy1.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.equ.overlap, observed.equ.two_sites, expected.equ.overlap, expected.equ.two_sites, pv_overlap.equ.p, pv_overlap.equ.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.anc_sit.overlap, observed.sit.overlap, expected.anc_sit.overlap, expected.sit.overlap, pv_overlap.anc_par.p, pv_overlap.anc_par.f, 1);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.asy_sit.overlap, observed.sit.overlap, expected.asy_sit.overlap, expected.sit.overlap, pv_overlap.asy2.p, pv_overlap.asy2.f, 1);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			//spacer
			fisher = fisher_exact_test(observed.anc.spacer, observed.anc.two_sites, expected.anc.spacer, expected.anc.two_sites, pv_spacer.anchor.p, pv_spacer.anchor.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.par.spacer, observed.par.two_sites, expected.par.spacer, expected.par.two_sites, pv_spacer.partner.p, pv_spacer.partner.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.asy.spacer, observed.asy.two_sites, expected.asy.spacer, expected.asy.two_sites, pv_spacer.asy1.p, pv_spacer.asy1.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.equ.spacer, observed.equ.two_sites, expected.equ.spacer, expected.equ.two_sites, pv_spacer.equ.p, pv_spacer.equ.f, 0);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.anc_sit.spacer, observed.sit.spacer, expected.anc_sit.spacer, expected.sit.spacer, pv_spacer.anc_par.p, pv_spacer.anc_par.f, 1);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
			fisher = fisher_exact_test(observed.asy_sit.spacer, observed.sit.spacer, expected.asy_sit.spacer, expected.sit.spacer, pv_spacer.asy2.p, pv_spacer.asy2.f, 1);
			if (fisher == -1)
			{
				fprintf(stderr, "Error: Fisher test error Anc 0 Par %d\n", mot);
				return -1;
			}
		}
		for (i = 0; i<5; i++)
		{
			char file_pval0[ARGLEN];
			strcpy(file_pval0, file_pval[i]);
			char buf[10];
			sprintf(buf, "%d", mot);
			strcat(file_pval0, buf);
			if ((out_pval[i] = fopen(file_pval0, "wt")) == NULL)
			{
				fprintf(stderr, "Error: Input file %s can't be opened!\n", file_pval0);
				return -1;
			}
			fprintf(out_pval[i], "Anchor Thr\tPartner Thr\t\tReal CE+\tReal Total\tRand CE+\tRand Total\tFold\tP-value\n");
		}
		for (j = 0; j<NUM_THR; j++)
		{
			for (k = 0; k<NUM_THR; k++)
			{
				for (i = 0; i<5; i++)fprintf(out_pval[i], "A %d\tP %d\t\t", j + 1, k + 1);
				fprintf(out_pval[0], "%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.cell[j][k].any, observed.cell[j][k].two_sites, expected.cell[j][k].any, expected.cell[j][k].two_sites, fold_a[j][k], pvalue_a[j][k]);
				fprintf(out_pval[1], "%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.cell[j][k].full, observed.cell[j][k].two_sites, expected.cell[j][k].full, expected.cell[j][k].two_sites, fold_f[j][k], pvalue_f[j][k]);
				fprintf(out_pval[2], "%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.cell[j][k].partial, observed.cell[j][k].two_sites, expected.cell[j][k].partial, expected.cell[j][k].two_sites, fold_p[j][k], pvalue_p[j][k]);
				fprintf(out_pval[3], "%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.cell[j][k].overlap, observed.cell[j][k].two_sites, expected.cell[j][k].overlap, expected.cell[j][k].two_sites, fold_o[j][k], pvalue_o[j][k]);
				fprintf(out_pval[4], "%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.cell[j][k].spacer, observed.cell[j][k].two_sites, expected.cell[j][k].spacer, expected.cell[j][k].two_sites, fold_s[j][k], pvalue_s[j][k]);
			}
		}		
		for (i = 0; i<5; i++)fprintf(out_pval[i], "\n");
		//any vs rand
		fprintf(out_pval[0], "Anchor\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc.any, observed.anc.two_sites, expected.anc.any, expected.anc.two_sites, pv_any.anchor.f, pv_any.anchor.p);
		fprintf(out_pval[0], "Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.par.any, observed.par.two_sites, expected.par.any, expected.par.two_sites, pv_any.partner.f, pv_any.partner.p);
		fprintf(out_pval[0], "Asymmetry\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.asy.any, observed.asy.two_sites, expected.asy.any, expected.asy.two_sites, pv_any.asy1.f, pv_any.asy1.p);
		fprintf(out_pval[0], "Symmetry\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.equ.any, observed.equ.two_sites, expected.equ.any, expected.equ.two_sites, pv_any.equ.f, pv_any.equ.p);
		//any vs real
		fprintf(out_pval[0], "Anchor_Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc_sit.any, observed.sit.any, expected.anc_sit.any, expected.sit.any, pv_any.anc_par.f, pv_any.anc_par.p);
		fprintf(out_pval[0], "Asym_Sym\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.asy_sit.any, observed.sit.any, expected.asy_sit.any, expected.sit.any, pv_any.asy2.f, pv_any.asy2.p);
		//full vs rand		
		fprintf(out_pval[1], "Anchor\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc.full, observed.anc.two_sites, expected.anc.full, expected.anc.two_sites, pv_full.anchor.f, pv_full.anchor.p);
		fprintf(out_pval[1], "Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.par.full, observed.par.two_sites, expected.par.full, expected.par.two_sites, pv_full.partner.f, pv_full.partner.p);
		fprintf(out_pval[1], "Asymmetry\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.asy.full, observed.asy.two_sites, expected.asy.full, expected.asy.two_sites, pv_full.asy1.f, pv_full.asy1.p);
		fprintf(out_pval[1], "Symmetry\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.equ.full, observed.equ.two_sites, expected.equ.full, expected.equ.two_sites, pv_full.equ.f, pv_full.equ.p);
		//full vs real
		fprintf(out_pval[1], "Anchor_Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc_sit.full, observed.sit.full, expected.anc_sit.full, expected.sit.full, pv_full.anc_par.f, pv_full.anc_par.p);
		fprintf(out_pval[1], "Asym_Sym\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.asy_sit.full, observed.sit.full, expected.asy_sit.full, expected.sit.full, pv_full.asy2.f, pv_full.asy2.p);
		//partial vs rand
		fprintf(out_pval[2], "Anchor\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc.partial, observed.anc.two_sites, expected.anc.partial, expected.anc.two_sites, pv_partial.anchor.p, pv_partial.anchor.p);
		fprintf(out_pval[2], "Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.par.partial, observed.par.two_sites, expected.par.partial, expected.par.two_sites, pv_partial.partner.f, pv_partial.partner.p);
		fprintf(out_pval[2], "Asymmetry\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.asy.partial, observed.asy.two_sites, expected.asy.partial, expected.asy.two_sites, pv_partial.asy1.f, pv_partial.asy1.p);
		fprintf(out_pval[2], "Symmetry\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.equ.partial, observed.equ.two_sites, expected.equ.partial, expected.equ.two_sites, pv_partial.equ.f, pv_partial.equ.p);
		//partial vs real
		fprintf(out_pval[2], "Anchor_Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc_sit.partial, observed.sit.partial, expected.anc_sit.partial, expected.sit.partial, pv_partial.anc_par.f, pv_partial.anc_par.p);
		fprintf(out_pval[2], "Asym_Sym\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.asy_sit.partial, observed.sit.partial, expected.asy_sit.partial, expected.sit.partial, pv_partial.asy2.f, pv_partial.asy2.p);
		//overlap vs rand
		fprintf(out_pval[3], "Anchor\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc.overlap, observed.anc.two_sites, expected.anc.overlap, expected.anc.two_sites, pv_overlap.anchor.f, pv_overlap.anchor.p);
		fprintf(out_pval[3], "Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.par.overlap, observed.par.two_sites, expected.par.overlap, expected.par.two_sites, pv_overlap.partner.f, pv_overlap.partner.p);
		fprintf(out_pval[3], "Asymmetry\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.asy.overlap, observed.asy.two_sites, expected.asy.overlap, expected.asy.two_sites, pv_overlap.asy1.f, pv_overlap.asy1.p);
		fprintf(out_pval[3], "Symmetry\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.equ.overlap, observed.equ.two_sites, expected.equ.overlap, expected.equ.two_sites, pv_overlap.equ.f, pv_overlap.equ.p);
		//overlap vs real
		fprintf(out_pval[3], "Anchor_Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc_sit.overlap, observed.sit.overlap, expected.anc_sit.overlap, expected.sit.overlap, pv_overlap.anc_par.f, pv_overlap.anc_par.p);
		fprintf(out_pval[3], "Asym_Sym\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.asy_sit.overlap, observed.sit.overlap, expected.asy_sit.overlap, expected.sit.overlap, pv_overlap.asy2.f, pv_overlap.asy2.p);
		//spacer vs rand
		fprintf(out_pval[4], "Anchor\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc.spacer, observed.anc.two_sites, expected.anc.spacer, expected.anc.two_sites, pv_spacer.anchor.f, pv_spacer.anchor.p);
		fprintf(out_pval[4], "Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.par.spacer, observed.par.two_sites, expected.par.spacer, expected.par.two_sites, pv_spacer.partner.f, pv_spacer.partner.p);
		fprintf(out_pval[4], "Asymmetry\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.asy.spacer, observed.asy.two_sites, expected.asy.spacer, expected.asy.two_sites, pv_spacer.asy1.f, pv_spacer.asy1.p);
		fprintf(out_pval[4], "Symmetry\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.equ.spacer, observed.equ.two_sites, expected.equ.spacer, expected.equ.two_sites, pv_spacer.equ.f, pv_spacer.equ.p);
		//spacer vs real
		fprintf(out_pval[4], "Anchor_Partner\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.anc_sit.spacer, observed.sit.spacer, expected.anc_sit.spacer, expected.sit.spacer, pv_spacer.anc_par.f, pv_spacer.anc_par.p);
		fprintf(out_pval[4], "Asym_Sym\t\t\t%d\t%d\t%d\t%d\t%.3f\t%g\n", observed.asy_sit.spacer, observed.sit.spacer, expected.asy_sit.spacer, expected.sit.spacer, pv_spacer.asy2.f, pv_spacer.asy2.p); 
		for (i = 0; i<5; i++)fclose(out_pval[i]);
		double pval_tot_min[5] = { 0, 0, 0, 0, 0 };
		double limit = 300;
		double pv_limit = 1E-300;
		double lgpv[5][NUM_THR][NUM_THR];
		for (j = 0; j<NUM_THR; j++)
		{
			for (k = 0; k<NUM_THR; k++)
			{
				for (i = 0; i<5; i++)lgpv[i][k][j] = 0;
				if (pvalue_a[k][j] <= pv_limit)lgpv[0][k][j] = limit;
				else lgpv[0][k][j] = -log10(pvalue_a[k][j]);
				if (pvalue_f[k][j] <= pv_limit)lgpv[1][k][j] = limit;
				else lgpv[1][k][j] = -log10(pvalue_f[k][j]);
				if (pvalue_p[k][j] <= pv_limit)lgpv[2][k][j] = limit;
				else lgpv[2][k][j] = -log10(pvalue_p[k][j]);
				if (pvalue_o[k][j] <= pv_limit)lgpv[3][k][j] = limit;
				else lgpv[3][k][j] = -log10(pvalue_o[k][j]);
				if (pvalue_s[k][j] <= pv_limit)lgpv[4][k][j] = limit;
				else lgpv[4][k][j] = -log10(pvalue_s[k][j]);
			}
		}
		for (i = 0; i<5; i++)
		{
			int limit_found = 0;
			int mnoj = NUM_THR*NUM_THR;
			for (j = 0; j<NUM_THR; j++)
			{
				for (k = 0; k<NUM_THR; k++)
				{
					double val = lgpv[i][k][j];
					if (val == limit)
					{
						limit_found = 1;
						break;
					}
				}
				if (limit_found == 1)break;
			}
			if (limit_found == 1)
			{
				pval_tot_min[i] = limit;
			}
			else
			{
				for (j = 0; j<NUM_THR; j++)
				{
					for (k = 0; k<NUM_THR; k++)
					{
						double val = lgpv[i][k][j];
						if (val>pval_tot_min[i])pval_tot_min[i] = val;
					}
				}
				//pval_tot_min[i]=pow(10,-pval_tot_min[i]);
			}
		}
		{
			pv_any.anc_par.p = -log10(pv_any.anc_par.p);
			pv_full.anc_par.p = -log10(pv_full.anc_par.p);
			pv_partial.anc_par.p = -log10(pv_partial.anc_par.p);
			pv_overlap.anc_par.p = -log10(pv_overlap.anc_par.p);
			pv_spacer.anc_par.p = -log10(pv_spacer.anc_par.p);
			if (pv_any.anc_par.p > bonferroni_corr_asy)
			{
				if (pv_any.anc_par.f < 1)pv_any.anc_par.p *= -1;
			}
			else pv_any.anc_par.p = 0;
			if (pv_full.anc_par.p > bonferroni_corr_asy)
			{
				if (pv_full.anc_par.f < 1)pv_full.anc_par.p *= -1;
			}
			else pv_full.anc_par.p = 0;
			if (pv_partial.anc_par.p > bonferroni_corr_asy)
			{
				if (pv_partial.anc_par.f < 1)pv_partial.anc_par.p *= -1;
			}
			else pv_partial.anc_par.p = 0;
			if (pv_overlap.anc_par.p > bonferroni_corr_asy)
			{
				if (pv_overlap.anc_par.f < 1)pv_overlap.anc_par.p *= -1;
			}
			else pv_overlap.anc_par.p = 0;
			if (pv_spacer.anc_par.p > bonferroni_corr_asy)
			{
				if (pv_spacer.anc_par.f < 1)pv_spacer.anc_par.p *= -1;
			}
			else pv_spacer.anc_par.p = 0;
		}
		{
			pv_any.asy2.p = -log10(pv_any.asy2.p);
			pv_full.asy2.p = -log10(pv_full.asy2.p);
			pv_partial.asy2.p = -log10(pv_partial.asy2.p);
			pv_overlap.asy2.p = -log10(pv_overlap.asy2.p);
			pv_spacer.asy2.p = -log10(pv_spacer.asy2.p);
			if (pv_any.asy2.p > bonferroni_corr_asy)
			{
				if (pv_any.asy2.f < 1)pv_any.asy2.p *= -1;
			}
			else pv_any.asy2.p = 0;
			if (pv_full.asy2.p > bonferroni_corr_asy)
			{
				if (pv_full.asy2.f < 1)pv_full.asy2.p *= -1;
			}
			else pv_full.asy2.p = 0;
			if (pv_partial.asy2.p > bonferroni_corr_asy)
			{
				if (pv_partial.asy2.f < 1)pv_partial.asy2.p *= -1;
			}
			else pv_partial.asy2.p = 0;
			if (pv_overlap.asy2.p > bonferroni_corr_asy)
			{
				if (pv_overlap.asy2.f < 1)pv_overlap.asy2.p *= -1;
			}
			else pv_overlap.asy2.p = 0;
			if (pv_spacer.asy2.p > bonferroni_corr_asy)
			{
				if (pv_spacer.asy2.f < 1)pv_spacer.asy2.p *= -1;
			}
			else pv_spacer.asy2.p = 0;
		}
		if ((out_pval_table = fopen(file_pval_table, "at")) == NULL)
		{
			fprintf(stderr, "Error: Input file %s can't be opened!\n", file_pval_table);
			return -1;
		}
		if (mot == 0)
		{
			fprintf(out_pval_table, "Anchor %d", mot);
			fprintf(out_pval_table, "\tAnchor");
		}
		else
		{
			fprintf(out_pval_table, "Partner %d", mot);
			fprintf(out_pval_table, "\t%s", name_partner);
		}
		for (i = 1; i<5; i++)
		{
			if (pval_tot_min[i]>bonferroni_corr)fprintf(out_pval_table, "\t%.2f", pval_tot_min[i]);
			else fprintf(out_pval_table, "\t0");
		}
		if (pval_tot_min[0]>bonferroni_corr)fprintf(out_pval_table, "\t%.2f", pval_tot_min[0]);
		else fprintf(out_pval_table, "\t0");
		{
			pv_full.anchor.p = -log10(pv_full.anchor.p);
			pv_full.partner.p = -log10(pv_full.partner.p);
			pv_full.asy1.p = -log10(pv_full.asy1.p);
			pv_full.equ.p = -log10(pv_full.equ.p);
			pv_partial.anchor.p = -log10(pv_partial.anchor.p);
			pv_partial.partner.p = -log10(pv_partial.partner.p);
			pv_partial.asy1.p = -log10(pv_partial.asy1.p);
			pv_partial.equ.p = -log10(pv_partial.equ.p);
			pv_overlap.anchor.p = -log10(pv_overlap.anchor.p);
			pv_overlap.partner.p = -log10(pv_overlap.partner.p);
			pv_overlap.asy1.p = -log10(pv_overlap.asy1.p);
			pv_overlap.equ.p = -log10(pv_overlap.equ.p);
			pv_spacer.anchor.p = -log10(pv_spacer.anchor.p);
			pv_spacer.partner.p = -log10(pv_spacer.partner.p);
			pv_spacer.asy1.p = -log10(pv_spacer.asy1.p);
			pv_spacer.equ.p = -log10(pv_spacer.equ.p);
			pv_any.anchor.p = -log10(pv_any.anchor.p);
			pv_any.partner.p = -log10(pv_any.partner.p);
			pv_any.asy1.p = -log10(pv_any.asy1.p);
			pv_any.equ.p = -log10(pv_any.equ.p);
			if (mot != 0)
			{				
				fprintf(out_pval_table, "\t%.2f\t%.2f\t%.2f", -log10(pvalue_similarity_tot), -log10(pval_sim[0]), -log10(pval_sim[1]));
				if (pv_overlap.anchor.p > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_overlap.anchor.p);
				else fprintf(out_pval_table, "\t0");
				if (pv_spacer.anchor.p > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_spacer.anchor.p);
				else fprintf(out_pval_table, "\t0");
				if (pv_overlap.partner.p > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_overlap.partner.p);
				else fprintf(out_pval_table, "\t0");
				if (pv_spacer.partner.p > bonferroni_corr_ap)fprintf(out_pval_table, "\t%.2f", pv_spacer.partner.p);
				else fprintf(out_pval_table, "\t0");
			}
			else fprintf(out_pval_table, "\t\t\t\t\t\t\t");
			if (pv_overlap.asy1.p != 0)fprintf(out_pval_table, "\t%.2f", pv_overlap.asy1.p);
			else fprintf(out_pval_table, "\t0");
			if (pv_spacer.asy1.p != 0)fprintf(out_pval_table, "\t%.2f", pv_spacer.asy1.p);
			else fprintf(out_pval_table, "\t0");
			if (pv_overlap.equ.p != 0)fprintf(out_pval_table, "\t%.2f", pv_overlap.equ.p);
			else fprintf(out_pval_table, "\t0");
			if (pv_spacer.equ.p != 0)fprintf(out_pval_table, "\t%.2f", pv_spacer.equ.p);
			else fprintf(out_pval_table, "\t0");
			if (mot != 0)
			{
				if (pv_overlap.anc_par.p != 0)fprintf(out_pval_table, "\t%+.2f", pv_overlap.anc_par.p);
				else fprintf(out_pval_table, "\t0");
				if (pv_spacer.anc_par.p != 0)fprintf(out_pval_table, "\t%+.2f", pv_spacer.anc_par.p);
				else fprintf(out_pval_table, "\t0");
			}
			else fprintf(out_pval_table, "\t\t");
			if (pv_overlap.asy2.p != 0)fprintf(out_pval_table, "\t%+.2f", pv_overlap.asy2.p);
			else fprintf(out_pval_table, "\t0");
			if (pv_spacer.asy2.p != 0)fprintf(out_pval_table, "\t%+.2f", pv_spacer.asy2.p);
			else fprintf(out_pval_table, "\t0");
			fprintf(out_pval_table, "\t%.2f", bonferroni_corr);
			fprintf(out_pval_table, "\t%.2f", bonferroni_corr_ap);
			fprintf(out_pval_table, "\t%.2f", bonferroni_corr_asy);
		}
		fprintf(out_pval_table, "\n");
		fclose(out_pval_table);
		{
			rand_one[ap].mem_out_sta();
			rand_one[ap].mem_out_cep();
			rand_one[ap].mem_out_cel();
			rand_one[ap].mem_out_pv();
			for (i = 0; i<nseq_rand; i++)rand_one[ap].nsit[i] = 0;
		}
		if (ap == 1)
		{
			real_one[ap].mem_out_sta();
			real_one[ap].mem_out_cep();
			real_one[ap].mem_out_cel();
			real_one[ap].mem_out_sco();
			real_one[ap].mem_out_pv();
			for (i = 0; i<nseq_real; i++)real_one[ap].nsit[i] = 0;
		}
		{
			FILE *out_log;
			if ((out_log = fopen(file_log, "wt")) == NULL)
			{
				fprintf(out_log, "Input file %s can't be opened!\n", file_log);
				return -1;
			}
			fprintf(out_log, "Calculations are completed for %d motifs out of total %d\n", mot + 1, n_motifs);
			fclose(out_log);
		}
	}
	real_plot.mem_out();
	rand_plot.mem_out();
	delete[] thr_err_real;
	delete[] thr_err_rand;
	delete[] peak_len_real;
	delete[] peak_len_rand;
	for (k = 0; k<2; k++)
	{
		for (i = 0; i<nseq_real; i++)
		{
			delete[] seq[k][i];
		}
		delete[] seq[k];
	}
	delete[] seq;
	for (i = 0; i<2; i++)real_one[i].mem_out_nsit();
	for (i = 0; i<2; i++)rand_one[i].mem_out_nsit();
	real_one[0].mem_out_sta();
	real_one[0].mem_out_cep();
	real_one[0].mem_out_cel();
	real_one[0].mem_out_sco();
	real_one[0].mem_out_pv();
	rand_hom_one.mem_out_sta();
	rand_hom_one.mem_out_cep();
	rand_hom_one.mem_out_cel();
	rand_hom_one.mem_out_pv();
	rand_hom_one.mem_out_nsit();
	matrix0.mem_out(matrix0.len);
	return 0;
}
