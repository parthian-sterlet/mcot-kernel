# Readme for MCOT software package

## General description

MCOT (Motifs Co-Occurrence Tool) is a software package for recognition of 
composite elements (CEs) in a single ChIP-seq dataset. CEs  detected by MCOT 
include two potential binding sites of transcription factors (TFs) in all 
possible mutual orientations. MCOT considers CEs with a full/partial overlap of 
motifs or with a spacer in a certain range. Each potential CE recognized by MCOT 
contains the motif of immunoprecipitated TF in respective ChIP-seq experiment 
(anchor motif) and another motif (partner). Identical/distinct anchor and 
partner motifs imply the search for CEs of homotypic or heterotypic type 
(respectively).

## Implementation

MCOT implemented in C++ and it can be conventionally compiled in Linux or 
Windows operating system. To run MCOT user should compile the corresponding 
source code file, files mcot_anchor.cpp and mcot.cpp respect to one-partner and 
many partners options.

## Installation
(Linux) Run in terminal (Packages “build-essential” and “cmake” 
should be installed on Ubuntu system):
```
cd <project>
mkdir tmp
cd tmp
cmake ..
make
```

All executable files for one-partner and many partners options will be in src/anchor\_vs\_one and src/anchor\_vs\_many

(Windows) Run in terminal (Win -> Visual Studio 2017 -> Visual Studio Tools -> 
VC -> Native Tools x64. Else, you should install “CMake” module while VS 2017 installing)

```
cd <project>
mkdir tmp
cd tmp
cmake ..
MSBuild mcot-kernel.sln /p:Configuration=Release /p:Platform=Win32
```
Programs `anchor_vs_many` and `anchor_vs_one` for one\_partner and many\_partners 
options should be in `src/anchor_vs_one_partner/Release/`
and `src/anchor_vs_many_partners/Release/`


## Command line arguments

The command line for one-partner option:


`./anchor_vs_one <1 fasta> <2 anchor.motif> <3 partner.motif> <4 minimal spacer length> <5 maximal spacer length> <6 path to whole-genome promoters>`


The command line for many-partner option:


`./anchor_vs_many <1 fasta> <2 anchor.motif> <3 partners.library> <4 minimal spacer length> <5 maximal spacer length> <6 path to whole-genome promoters>`



`<fasta>` = DNA sequences of ChIP-seq peaks in fasta format


`<anchor.motif>`, `<partner.motif>` = frequency matrices of motifs in standard format, e.g.  

\***

```
> Motif name
0.001   0.001    0.942    0.056
0.001    0.001    0.997    0.001
0.001    0.001    0.997    0.001
0.548    0.001    0.400    0.051
0.589    0.087    0.245    0.079
0.231    0.001    0.001    0.767
0.001    0.001    0.001    0.997
0.001    0.381    0.001    0.617
0.001    0.964    0.001    0.034
0.001    0.891    0.001    0.107
***

```

`<partners.library>` = for this parameter five options are available: “hs_core”, “mm_core”, “hs_full”, “mm_full” and “dapseq”. These values respect to libraries of human/mouse core (396/353) and full (747/509) collections of motifs (Kulakovskiy et al., 2018, http://hocomoco11.autosome.ru/); and the library of 514 motifs from Plant Cistrome of A. thaliana motifs (http://neomorph.salk.edu/dap_web/pages/index.php, O’Malley et al., 2016).


`<minimal spacer length>` = integer value from 0 to \<maximal spacer length>  (the default value 0 is recommended, any positive value restricts short spacers)



`<maximal spacer length>` = integer value from 0 to 100 (the default value 29)



`<path to whole-genome promoters>` =  a path to the whole-genome dataset of promoters, three folders “hs”, “mm” and “at” that imply application of H.sapiens, M.musculus and A.thaliana promoter datasets for setting of thresholds for input motifs.




## Input data

MCOT requires (a) DNA sequences of ChIP-seq peaks and (b) anchor and partner motifs. We recommend application of a conventional de novo motif search tool, e.g. HOMER (Heinz et al., 2010) to define the anchor motif. 

MCOT have two options for definition of the partner motif:


* User defines a matrix of partner motif (one-partner option);

* User points to a library of known motifs (many partners option).

MCOT allows the variation of the upper limit of spacer length from zero to 100 base pairs.

## Motifs recognition

MCOT applies the recognition model of Position Weight Matrix (PWM) for mapping motifs  in peaks. For each matrix, MCOT uses five thresholds {T[1],...T[5]} according to the unified set of five expected false positive rates (FPR) for a whole-genome dataset of promoters, {5.24E-5, 1.02E-04, 1.9E-4, 3.33E-4, 5E-4}. The profile of the most stringent hits contains matrix scores T ≥ T[1], the next profile comprises PWM scores {T} in the range T[2] ≥ T > T[1], etc. Hence, MCOT computes five profiles of hits with certain level of conservation for each input motif.




## Composite elements search and annotation

MCOT classifies CEs structure according to the following criteria:
* _Orientation_. Four types of distinct mutual orientations are considered: in the same DNA strand (Direct Anchor/Partner and Direct Partner/Anchor), in opposite strands (Everted and Inverted);
* _Overlap or Spacer_. There are three distinct cases of mutual locations: Full overlap (one motif located entirely within another one); Partial overlap; and Spacer. To describe each case MCOT uses the following characteristics: the distance between nearest borders of two motifs (Full); the length of overlapped region (Partial); and the length of spacer;
* _Asymmetry of сonservation_.  All predicted CEs are subdivided into two classes: those with more conservative anchor and partner motifs. The conservation of motif hit is estimated as  -Log10(FPR), where FPR is computed by the respective score of recognition model

For each of 25 combinations of motifs conservation and each computation flow MCOT compiles the 2x2 contingency table (Table 1) and compute the significance of Fisher’s exact test p-value(CE) that compares the content of CEs in ChIP-seq peaks with that for background model.

**Table 1.** 2x2 contingency table for calculation of CE significance.

|                             | Sequences without CEs               | Sequences with CEs  | Total, sequences with both motifs|
|-----------------------------|-------------------------------------|---------------------|----------------------------------|
|Observed (peaks)             |$`Obs_{CE-} = Obs_{Tot} - Obs_{CE+}`$|$`Obs_{CE+}`$        | $`Obs_{Tot}`$                    |
|Expected (permuted sequences)|$`Exp_{CE-} = Exp_{Tot} - Exp_{CE+}`$|$`Exp_{CE+}`$        | $`Exp_{Tot}`$                    |

The background model implies the preservation of content for each motif in each peak. MCOT generates background profiles of hits of anchor and partner motifs iteratively for each peak with a special permutation procedure.

MCOT computes CE significance separately for five computation flows:  Any (Spacer or Overlap), Full (overlap), Partial (overlap), Overlap (full and partial), Spacer.

Fisher’s exact test computes the CE enrichment for each of 5x5 combinations of conservation of motifs. Since MCOT checks five conservation ranges for each motifs, the Bonferroni-corrected

p-value < 0.05/(5*5) = 0.002 is used to mark the adjusted threshold of CE significance in each computation flows. I.e. the fall of the minimal p-value among 25 estimates below 0.002 implied the significance. If simultaneously 396/353 matrices of human/mouse TFs were applied (option ‘many partners’) then p-value threshold fell below to 0.002/396 ≈ 5E-6.


Finally, MCOT estimates the similarity of anchor and partner motifs with the motifs similarity p-value. MCOT used matrix column similarities measures PCC (Pearson Correlation Coefficient, Pietrokovski, 1996) and SSD (Sum of Squared Distances, Wasserman and Sandelin, 2004) to compute two p-values. If maximum among them was less than 0.05, then pair of motifs had the significant similarity. Hence, a respective CE may be a possible false positive prediction.

Additionally, MCOT computed significance of CEs with more conserved anchor motif and more conserved partner motif. These calculations are performed according to scheme represented above in Table 1, but MCOT counts either CEs with more conserved Anchor or Partner motifs, i.e. either -Log10[FPR(Anchor)] > -Log10[FPR(Partner)] or -Log10[FPR(Anchor)] ≤ -Log10[FPR(Partner)].

Finally, to estimate the asymmetry of CE, MCOT partitions all CEs on those with more conserved Anchor or Partner motifs and compute the significance of asymmetry with Fisher’s exact test (Table 2). The asymmetry significance -Log10[P-value] is printed as positive (with sign ‘+’) in the case of enrichment toward the Anchor motif, otherwise, MCOT sign ‘-’ (negative value) denotes the enrichment toward the Partner motif.  These calculation yielded for each computation flow one p-value(CE, Asymmetry).

**Table 2.** 2x2 contingency table for calculation of CE asymmetry.

|                                    | CEs with more conserved anchor motif| CEs with more conserved partner motif  | CEs, total                              |
|------------------------------------|-------------------------------------|----------------------------------------|-----------------------------------------|
|Observed (CEs in peaks)             |$`Obs_{CE, Anchor}`$                 |$`Obs_{CE, Partner}`$                   | $`Obs_{CE, Anchor} + Obs_{CE, Partner}`$|
|Expected (CEs in permuted sequences)|$`Obs_{CE, Anchor}`$                 |$`Obs_{CE, Partner}`$                   | $`Obs_{CE, Anchor} + Obs_{CE, Partner}`$|


Detailed enrichment or depletion of CEs with specific combinations of motifs conservation are represented in the scatterplot text file plot_*, see below.

## Output data

MCOT gives the following output data:



* __Files <*_thr5>, recognition profiles of motifs__ Each file respects to one motif. A file has fasta-like format, i.e. for each peak the header line starts with ‘>’ symbol. Next, each subsequent line represents one hit in a peak, particularly it position, respective conservation value -Log10(FPR) and DNA strand.

```
Example
>Seq 1  Thr 0.864497    Nsites 1
205 3.304404    -
>Seq 2    Thr 0.864497    Nsites 0
>Seq 3    Thr 0.864497    Nsites 2
88    3.607778    -
160    3.338550    -
```

Here and below ChIP-seq data for mouse FoxA2 and CE FoxA2-HNF1β (Wederell et al., 2008) illustrate MCOT application. The anchor FoxA2 motif we deduced from de novo search (Homer, Heinz et al., 2009) and the partner HNF1β motif we took from the mouse Hocomoco core collection (Kulakovskiy et al., 2018).



* __File <rec_pos.txt>, the detailed recognition statistics__. For each motif and each recognition threshold MCOT provides (1) the number and the name of the motif (anchor motif is designated as ‘Anchor’; numbers 1,2, ... belong to partner motifs), (2) the number and the value of the threshold; (3) the percentage of peaks containing at least one hit of the motif, the number of peaks with recognized motif and the total number of peaks, (4) the number of recognized hits per base pair, the number of recognized hits and the total number of available locations for the motif.



Example


| Motif Num | Motif Name | Threshold | % of peaks | Rec. peaks | Total peaks | Rate of hits | Rec. hits | Total positions |
|-----------|------------|-----------|------------|------------|-------------|--------------|-----------|-----------------|
|0          | Anchor     | 0.955618  | 32.93      | 2484       | 7543        | 7.23E-04     | 2887      | 3994992         |
|0          | Anchor     | 0.945548  | 19.86      | 1498       | 7543        | 4.20E-04     | 1677      | 3994992         |
|0          | Anchor     | 0.934918  | 25.30      | 1908       | 7543        | 5.50E-04     | 2199      | 3994992         |
|0          | Anchor     | 0.923776  | 29.76      | 2245       | 7543        | 6.69E-04     | 2672      | 3994992         |
|0          | Anchor     | 0.913846  | 28.85      | 2176       | 7543        | 6.63E-04     | 2647      | 3994992         |
|1          |Partner     | 0.959810  | 9.35       | 705        | 7543        | 1.90E-04     | 760       | 4002535         |
|1          |Partner     | 0.947124  | 7.42       | 560        | 7543        | 1.49E-04     | 596       | 4002535         |
|1          |Partner     | 0.934158  | 11.93      | 900        | 7543        | 2.44E-04     | 978       | 4002535         |
|1          |Partner     | 0.923150  | 16.56      | 1249       | 7543        | 3.51E-04     | 1403      | 4002535         |
|1          |Partner     | 0.959810  | 18.60      | 1403       | 7543        | 3.94E-04     | 1575      | 4002535         |


* __File <out_pval>, the summary for statistical significances for all pairs of anchor-partner motifs__ represents the calculation results for different potential CE variants: a homotypic CE (# Motif = Anchor) and  one/several heterotypic CE(s) (Partner 1, Partner 2, etc.) for five computation flows (Full/Partial overlap, Overlap, Spacer and Any).

The first block of output data for each pair of motifs represents (a) five P-values of CE enrichment in five computation flows; (b) for each heterotypic pair three P-value for similarity between Anchor and Partner motifs.

__TODO__ 
Example:




|Partner Num|   Partner Name | Any Pv | Full Pv | Partial Pv | Overlap Pv | Spacer Pv |    |Sim   |Sim SSD| Sim PCC |
|-----------|----------------|--------|---------|------------|------------|-----------|----|------|-------|---------|
|0          | Anchor         |7.8E-16 | 1       | 5.9E-08    | 3.2E-07    | 4.0E-13   |    |      |       |         |
|1          | Partner 1      |4.6E-10 | 1       | 4.0E-09    | 3.0E-06    | 6.6E-05   |    | 0.09 |0.09   | 1       |


The second block of data consequently for five computation flows (Any, Full, etc.) represents significances for CEs with more conserved anchor motif, for those with more conserved partner motif and for those with similar conservation of motifs.

Example:



|Partner Number|..|Anchor Any|Partner Any|Equal Any| | Anchor Full|Partner Full|Equal Full|etc.|
|--------------|--|----------|-----------|---------|-|------------|------------|----------|----|
|0             |  |5.0E-11   |7.2E-11    |6.6E-07  | |2           |2           | 1        |    |
|1             |  |1.9E-22   |1.2E-53    |1.3E-22  | |3.4E-45     |3.4E-171    |4.1E-23   |    |


The third block of data consequently for five computation flows (Any, Full, etc.) represents three significances (-Log10[p-value]) of asymmetry in CEs: ‘Anchor vs. Partner’, ‘Anchor vs Equal’ and ‘Partner vs. Equal’. The positive/negative values reflect an enrichment toward the first/second participant in a comparison

Example:

|Partner Number|..|Anchor Partner Any|Anchor Equal Any|Partner Equal Any| |Anchor Partner Full|Anchor Equal Full|Partner Equal Full|etc.|
|--------------|--|------------------|----------------|-----------------|-|-------------------|-----------------|------------------|----|
|0             |  |0.00              |0.41            |0.32             | |-0.30              |-0.30            |-0.30             |    |
|1             |  |-2.57             |1.16            |6.44             | |-10.14             |2.70             |22.44             |    |


* __File <out_hist>, the abundance of various CE types as a function of mutual orientation and location of the motifs__ The percentage of peaks containing CE variants specific in mutual orientation (four types) and mutual locations from a few possible full overlaps (‘F’), through a variety of partial overlaps (‘P’) and finally from the minimal to the maximal spacer length (‘S’).

Example:




|              | 1F | 10P | 9P | 8P | 7P | 6P | 5P | 4P | 3P | 2P | 1P | 0S | 1S | 2S |
|--------------|----|-----|----|----|----|----|----|----|----|----|----|----|----|----|
| Everted      |    |     |    |    |    |    |0.93|0.03|0.13|0.20|0.19|0.09|0.19|0.08|
| Inverted     |    |     |    |0.72|    |0.04|0.12|0.04|0.16|0.20|0.19|0.19|0.20|0.17|
| DirectPA     |    |     |    |    |    |    |    |0.05|0.01|0.08|0.04|0.07|0.09|0.17|
| DirectAP     |    |     |    |3.31|    |0.03|0.34|0.15|0.11|0.30|0.23|0.13|0.08|0.16|


* __Files <fisher\_*>, the 2x2 tables of CE significance for five computation flows for all motifs pairs and for  all 5x5 combination of motifs conservation__. Each line of output file contains data concerning one 2x2 contingency table, in particular (1) the designation of conservation (PWM threshold (T), indices from 1 to 5 mean the change from the most stringent to the most permissive, see above); (2) four counts for 2x2 contingency table, ‘the number of peaks containing at least one CE (CE+) & ‘the number of peaks containing at least one hit of each motif (Total)’ for peaks (Real) and permuted (Rand) datasets. Finally, the table contains significance of CEs (p-values) computed by Fisher’s exact test for 25 cells of 5x5 tables of combinations of thresholds. Next, the respective data are shown for  (a) significances of CEs with more conserved anchor motif, with more conserved partner motif and with balanced conservation of two motifs; and (b) significances of asymmetry in CEs ‘Anchor vs. Partner’, Anchor vs. Equal’ and ‘Partner vs. Equal’, in the latter three cases the Fold for each Fisher’s test are also provided. For example, for ‘Anchor vs. Partner’ test criteria fols>1 & fold<1 imply the enrichment of CEs with more conserved anchor & partner motifs.

Example



|Anchor Thr|Partner Thr|  |Real CE+|Real Total|Rand CE+|Rand Total|Fold|P-value |
|----------|-----------|--|--------|----------|--------|----------|----|--------|
|1         |1          |  |5       |135       |298     |5400      |    |1       |
|1         |2          |  |4       |157       |310     |6480      |    |1       |
|1         |3          |  |5       |202       |382     |8080      |    |1       |
|1         |4          |  |23      |290       |622     |11600     |    |0.08    |
|1         |5          |  |61      |304       |637     |12160     |    |1.62E-18|
|2         |1          |  |2       |123       |252     |4920      |    |1       |
|2         |2          |  |2       |121       |211     |4840      |    |1       |
|2         |3          |  |7       |168       |318     |6720      |    |1       | 
|2         |4          |  |35      |234       |472     |9360      |    |3.27E-08|
|2         |5          |  |45      |230       |501     |9200      |    |4.59E-13|
|3         |1          |  |7       |196       |406     |7840      |    |1       |
|3         |2          |  |8       |184       |413     |7360      |    |1       |
|3         |3          |  |19      |278       |560     |11120     |    |0.21    |
|3         |4          |  |42      |386       |865     |15440     |    |1.07E-04|
|3         |5          |  |61      |383       |839     |15320     |    |4.67E-13|
|4         |1          |  |6       |209       |492     |8360      |    |1       |
|4         |2          |  |3       |184       |359     |7360      |    |1       |
|4         |3          |  |4       |248       |562     |9920      |    |1       |
|4         |4          |  |9       |346       |724     |13840     |    |1       |
|4         |5          |  |25      |353       |785     |14120     |    |0.28    |
|5         |1          |  |4       |192       |464     |7680      |    |1       |
|5         |2          |  |2       |203       |468     |8120      |    |1       |
|5         |3          |  |6       |296       |663     |11840     |    |1       |
|5         |4          |  |23      |424       |979     |16960     |    |1       |
|5         |5          |  |31      |384       |824     |15360     |    |0.04    |
|          |           |  |        |          |        |          |    |        |
|Anchor    |           |  |172     |1286      |3637    |51440     |    |6.82E-15|
|Partner   |           |  |54      |1172      |3511    |46880     |    |1       |
|Equal     |           |  |66      |1062      |2523    |42480     |    |0.74    |
|Anchor    | Partner   |  |172     |1286      |54      |1172      |2.90|2.42E-14|
|Anchor    | Equal     |  |172     |1286      |66      |1062      |2.15|8.46E-09|
|Partner   | Equal     |  |54      |1172      |66      |1062      |0.74|0.13    |



* __Files <*.best>, the list of predicted CEs__. For each recognized CE MCOT provides (1) the header of a peak, (2) the start and the end positions of each motif in a peak, (3) mutual location (Full / Partial / Spacer types and the respective base pair measures, see above), (4) the strands of the motifs in a peak and the mutual orientation of the motifs (one of four types), (5) PWM scores and DNA sequences of the motifs.



|\#Seq    |Anchor start|Anchor end|Partner start|Partner end|Mutual Loc|Loc Type|Strands|Mutual Ori|Anchor score|Partner score|Anchor seq  |Partner seq  | 
|---------|------------|----------|-------------|-----------|----------|--------|-------|----------|------------|-------------|------------|-------------|
|Seq 491  |135         |146       |134          |146        |0F        |Full    |+-     |Evert     |0.918       |0.937        |ccaggatgtcaa|ttgacatcctggg|
|Seq 493  |354         |365       |378          |390        |12S       |Spacer  |+-     |Invert    |0.936       |0.955        |tgaggaagtgaa|ttgacattcttcc|
|Seq 513  |156         |167       |148          |160        |5P        |Partial |--     |DirectAP  |0.949       |0.959        |agaggaaatgac|atgacagattggg|


* __Files <plot_*>, scatterplots that show the CE asymmetry, i.e. the abundance of CEs with various ratios of conservation of Anchor and Partner motifs__ For each of five computation flows (Full, Partial, Overlap, Spacer and Any) one scatterplot is computed. For foreground and background data (peaks and sequences with permuted hits) the respective list of predicted CEs are subdivided on two fractions: those with more conserved Anchor and Partner motifs. The conservation of a hit is estimated with the respective -Log10(FPR) value. The minimal conservation value is equal to -Log10(5E-4) ~ 3.3. Next, MCOT computes two matrices {OBSi,j} and {EXPi,j} of absolute frequencies of conservation of Anchor and Partner motifs for observed and expected data. Here indices i and j imply the conservation -Log10(FPR) of Anchor and Partner motifs. This conservation is falling within the ranges [<3.5], [3.5..3.7], [3.7..3.9] etc. up to [5.3..5.5] and [>5.5]. Finally, the per mille measure transforms the absolute frequencies to relative ones as follow:
{1000\*OBSi,j/N(OBS)} and {1000\*EXPi,j/N(EXP)}, where N(OBS) and N(EXP) are total counts of predicted CEs in observed and expected lists. The output scatterplot shows the difference between relative frequencies of observed and expected CEs. Example below shows the asymmetry toward the Partner (HNF1B) motif in Overlap computation flow.

Example


```
,<3.5,3.5..3.7,3.7..3.9,3.9..4.1,4.1..4.3,4.3..4.5,4.5..4.7,4.7..4.9,4.9..5.1,5.1..5.3,5.3..5.5,>5.5
<3.5,14,14,-8,-14,-26,-12,1,,,,,
3.5..3.7,7,,-15,,-7,-7,,,,,,-2
3.7..3.9,37,-4,-2,-6,4,-11,-3,,,,,-3
3.9..4.1,17,-1,14,1,-4,-9,,,,,,
4.1..4.3,39,,-2,8,5,-8,-2,,,,,
4.3..4.5,,3,-6,4,5,5,,,,,,
4.5..4.7,,-10,22,,,-5,,,,,,
4.7..4.9,13,-4,,,-3,,,,,,,
4.9..5.1,-23,,-3,,5,-4,1,,,,,
5.1..5.3,1,,,,-2,,,,,,,
5.3..5.5,-4,-3,-2,,-3,,,,,,,
>5.5,4,2,-1,-2,-3,,,,,,,
```

## References

[Levitsky V.G., Zemlyanskaya E.V., Oshchepkov D.Yu., Podkolodnaya O.A., Ignatieva E.V., Grosse I., Mironova V.V., Merkulova T.I. A single ChIP-seq dataset is sufficient for comprehensive analysis of motifs co-occurrence with MCOT package. Nucleic Acids Research, 2019](https://doi.org/10.1093/nar/gkz800)

[Heinz,S., Benner,C., Spann,N., Bertolino,E., Lin,Y.C., Laslo,P., Cheng,J.X., Murre,C., Singh,H. and Glass,C.K. (2010) Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Mol Cell, 38, 576-589.](https://doi.org/10.1016/j.molcel.2010.05.004)

[Kulakovskiy,I.V., Vorontsov,I.E., Yevshin,I.S., Sharipov,R.N., Fedorova,A.D., Rumynskiy,E.I., Medvedeva,Y.A., Magana-Mora,A., Bajic,V.B., Papatsenko,D.A., et al. (2018) HOCOMOCO: expansion and enhancement of the collection of transcription factor binding sites models. Nucleic Acids Res., 46, D252-D259.](https://academic.oup.com/nar/article/46/D1/D252/4616875)

[O'Malley,R.C., Huang,S.C., Song,L., Lewsey,M.G., Bartlett,A., Nery,J.R., Galli,M., Gallavotti,A., Ecker,R. (2016) Cistrome and epicistrome features shape the regulatory DNA landscape. Cell,165, 1280-1292.](https://doi.org/10.1016/j.cell.2016.08.063)

[Pietrokovski,S. (1996) Searching databases of conserved sequence regions by aligning protein multiple-alignments. Nucleic Acids Res., 24, 3836-3845.](https://academic.oup.com/nar/article/24/19/3836/2384639)

[Sandelin,A., Wasserman,W.W. (2004) Constrained binding site diversity within families of transcription factors enhances pattern discovery bioinformatics. J Mol Biol., 338, 207-215.](https://doi.org/10.1016/j.jmb.2004.02.048)

[Wederell,E.D., Bilenky,M., Cullum,R., Thiessen,N., Dagpinar,M., Delaney,A., Varhol R, Zhao Y, Zeng T, Bernier B, et al. (2008). Global analysis of in vivo Foxa2-binding sites in mouse adult liver using massively parallel sequencing. Nucleic Acids Res., 36, 4549-4564.](https://doi.org/10.1093/nar/gkn382) 
