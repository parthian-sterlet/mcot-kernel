# Readme for MCOT software package

## General description

MCOT (Motifs Co-Occurrence Tool) is a software package for recognition of 
composite elements (CEs) in a single ChIP-seq dataset. CEs  detected by MCOT 
include  two transcription factor (TF) binding sites in all possible mutual 
orientations, whose motifs are overrepresented in a given ChIP-seq dataset. 
MCOT considers CEs with an overlap of motifs or with a spacer. Each potential 
CE recognized by MCOT contains the motif of immunoprecipitated TF in respective 
ChIP-seq experiment (anchor motif) and another motif (partner). 
Identical/distinct anchor and partner motifs imply the search for CEs of 
homotypic or heterotypic type (respectively).

## Implementation and installation

MCOT implemented in C++ and it can be conventionally compiled in Linux or 
Windows operating system. To run MCOT user should compile the corresponding 
source code file, files mcot_anchor.cpp and mcot.cpp respect to one-partner and 
many partners options.

(Linux) Run in terminal (Packages “build-essential” and “cmake” 
should be installed on Ubuntu system):
```
cd <project>
mkdir tmp
cd tmp
cmake ..
make
```

All executable files will be in tmp/src/anchor_vs_<one_partner, many_partners>/anchor_vs_<one_partner, many_partners>/

(Windows) Run in terminal (Win -> Visual Studio 2017 -> Visual Studio Tools -> 
VC -> Native Tools x64. Else, you should install “CMake” module while VS 2017 installing)

```
cd <project>
mkdir tmp
cd tmp
cmake ..
MSBuild mcot-kernel.sln /p:Configuration=Release /p:Platform=Win32
```
Programs should be in tmp/src/anchor_vs_<one_partner, many_partners>/Release/

## Command line arguments

The command line for one-partner option:


`./anchor_vs_one_partner <1 fasta> <2 anchor.motif> <3 partner.motif> <4 minimal spacer length> <5 maximal spacer length> <6 path to whole-genome promoters>`


The command line for many-partner option:


`./anchor_vs_many_partners <1 fasta> <2 anchor.motif> <3 partners.library> <4 minimal spacer length> <5 maximal spacer length> <6 path to whole-genome promoters>`



`<fasta>` = DNA sequences of ChIP-seq peaks in fasta format


`<anchor.motif>`, `<partner.motif>` = frequency matrices of motifs in standard format, e.g.  ***

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


`<minimal spacer length> = integer value from 0 to <maximal spacer length>  (the default value 0 is recommended, any positive value restricts short spacers)`


`<upper limit of spacer length>` = integer value from 0 to 100 (the default value 30 is recommended)


`<path to whole-genome promoters>` = the distribution contains three folders “hs”, “mm” and “at” that imply application of H.sapiens, M.musculus and A.thaliana whole genome promoter datasets for setting of thresholds for input motifs.


## Input data

MCOT requires (a) DNA sequences of ChIP-seq peaks and (b) anchor and partner motifs. We recommend application of a conventional de novo motif search tool, e.g. HOMER (Heinz et al., 2010) or ChIPMunk (Kulakovskiy et al., 2010) to define the anchor motif. 

MCOT have two options for definition of the partner motif:

* User defines a matrix of partner motif (one-partner option);

* User points to a library of known motifs (many partners option).

MCOT allows the variation of the upper limit of spacer length from zero to 100 base pairs.

## Motifs recognition

MCOT applies the recognition model of Position Weight Matrix (PWM) for mapping motifs  in peaks. For each matrix, MCOT uses five thresholds {t1,...t5} according to the unified set of five expected false positive rates for a whole-genome dataset of promoters, {5E-4, 3.33E-4, 1.9E-4, 1.02E-04, 5.24E-5} . The profile of the most stringent hits contains matrix scores t ≥ t1, the next profile comprises PWM scores in the range t2 ≥ t > t1, etc. Each profile respects to certain level of a motif conservation.



## Composite elements search and annotation

MCOT classifies CEs structure according to the following criteria:
* _The same or opposite  DNA strands_. Four types of distinct mutual orientations are considered: in the same DNA strand (Direct Anchor/Partner (1) and Direct Partner/Anchor (2)), in opposite DNA strands (Everted (3)  and Inverted (4));
* _Overlap or Spacer_. There are three distinct cases of mutual locations: (1) Full overlap (one motif located entirely within another one); (2) Partial overlap; and (3) Spacer. To describe each case MCOT uses the following characteristics: (1) the distance between nearest borders of two motifs; (2) the number of nucleotides located at the intersection of motifs; and (3) the length of spacer;
* _Asymmetry of similarity to recognition model_.  25 variants of CEs with a certain ratio of similarities of the anchor and partner motifs to their recognition models follow from 5x5 combinations of conservations of two motifs.


For each combination of motifs conservation and each computation flow  MCOT assigns  the p-value of Fisher’s exact test that compares the content of CEs in ChIP-seq peaks with that for background model. The background model implies the preservation of content for each motif in each peak. MCOT generates background profiles of hits of anchor and partner motifs iteratively for each peak with a special permutation procedure.

MCOT computes CE enrichment separately for five computation flows: Full overlap, Partial overlap, any Overlap (full and partial), Spacer, Any co-occurrence of the motifs (with a spacer or with an overlap).

Fisher’s exact test computes the CE enrichment for each of 5x5 combinations of stringencies of recognition models. Since MCOT checks five thresholds for each motifs, the Bonferroni-corrected 
 p-value<0.05/(5*5) = 0.002 is used to mark the adjusted threshold of CE significance in each computation flows. 

Finally, MCOT estimates  the similarity of anchor and partner motifs with the motifs similarity p-value. The criterion p-value < 0.05 supports the hypothesis that certain pair of motifs has the significant similarity and marks the respective pair as a possible false positive prediction.

## Output data

MCOT gives the following output data:

* __File <rec_pos.txt>, detailed recognition statistics__. For each motif and each recognition threshold (five in total) MCOT provides (1) the number* and the name of the motif, (2) the percentage of peaks containing at least one hit of the motif, the number of peaks with recognized motif and the total number of peaks, (4) the number of recognized hits per base pair, the number of recognized hits and the total number of available locations for the motif.

*If the option ‘one partner’ is applied, than 0,1 are the numbers of the motifs. If the option ‘many partners’ is applied, than 0 is an anchor motif and 1,2, ... is the list of potential partners to be tested. 

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


* __File <out_pval>, the summary for statistical significances for all pairs of anchor-partner motifs__ represents the calculation results for different potential CE variants: a homotypic CE,  one or several heterotypic CE(s) depending on the option one/many partners. Any - any co-occurrence of the motifs, Full - full overlap, Part - partial overlap, Over - Partial or Full overlap, Spac - spacer. The list contains (a) for each pair of motifs five p-values (Pv) of CE enrichment in five computation flows; (b) for each heterotypic pair the p-value for similarity of the anchor and the partner motifs. Next, for each computation flow the asymmetry (Asy) coefficient reflects to the tendency to have more conserved either the anchor or the partner motif within the CE.

Example for one-partner option:




|Motif Num|   Motif Name | Any Pv | Full Pv | Part Pv | Over Pv | Spac Pv |    | Any Asy | Full Asy |  Part Asy | Over Asy | Spac Asy |   | Sim Pv |
|---------|--------------|--------|---------|---------|---------|---------|----|---------|----------|-----------|----------|----------|---|--------|
|0,0      | First First  |7.8E-16 | 1       |5.9E-08  | 3.2E-07 | 4.0E-13 |    | 0.47    |0.00      | -0.07     | -0.06    | 0.95     |   |        |
|0,1      | First Second |4.6E-10 |1        | 4.0E-09 | 3.0E-06 | 6.6E-05 |    | 0.43    | 0.00     | 0.00      | -0.03    | 0.19     |   | 1      |
|1,1      | Second Second|5.6E-02 |1        |7.8E-03  | 1.0E-02 | 4.3E-02 |    | 0.00    | 0.03     | 0.09      | 0.08     | 0.04     |   |        |

 

Example for many partners option:



|Motif Num|   Motif Name | Any Pv | Full Pv | Part Pv | Over Pv | Spac Pv |    | Any Asy | Full Asy |  Part Asy | Over Asy | Spac Asy |   | Sim Pv |
|---------|--------------|--------|---------|---------|---------|---------|----|---------|----------|-----------|----------|----------|---|--------|
|0,0      | Anchor       |8.1E-16 | 1       | 2.7E-07 | 5.6E-07 | 2.8E-14 |    | 0.52    | 0.00     | -0.09     | -0.09    | 0.82     |   |        |
|0,1      | Partner 1    |4.7E-02 | 1       | 0.21    | 0.92    | 3.0E-02 |    | 0.01    | 0.00     | -0.01     | 0.00     | 0.31     |   | 1      |
|0,2      | Partner 2    |0.59    | 0.29    | 0.36    | 0.16    | 0.11    |    | 0.03    | -0.01    | -0.04     | -0.05    | 0.09     |   | 1      |



* __File <out_hist>, the abundance of various CE types as a function of mutual orientation and location of the motifs__ The percentage of peaks containing CE variants specific in mutual orientation (four types) and mutual locations from a few possible full overlaps (‘F’), through a variety of partial overlaps (‘P’) and finally from zero to the maximal spacer length (‘S’).

Example:




|0,1  |              | 1F | 10P | 9P | 8P | 7P | 6P | 5P | 4P | 3P | 2P | 1P | 0S | 1S | 2S |
|-----|--------------|----|-----|----|----|----|----|----|----|----|----|----|----|----|----|
|     | Everted      |    |     |    |    |    |    |0.93|0.03|0.13|0.20|0.19|0.09|0.19|0.08|
|     | Inverted     |    |     |    |0.72|    |0.04|0.12|0.04|0.16|0.20|0.19|0.19|0.20|0.17|
|     | DirectPA     |    |     |    |    |    |    |    |0.05|0.01|0.08|0.04|0.07|0.09|0.17|
|     | DirectAP     |    |     |    |3.31|    |0.03|0.34|0.15|0.11|0.30|0.23|0.13|0.08|0.16|


* __Files <fisher\_*>, the 2x2 tables of CE significance for five computation flows for all motifs pairs and for  all 5x5 combination of motifs conservation__. Each line of output file contains data concerning one 2x2 contingency table, in particular (1) the designation of conservation (PWM threshold (t), 1/5 are the most stringent/permissive) for each motif (see above); (2) four counts for 2x2 contingency table, ‘the number of peaks containing at least one CE (CE+) & ‘the number of peaks containing at least one hit of each motif (Total)’ for peaks (Real) and permuted (Rand) datasets. Finally, the table contains the fold and p-value computed by Fisher’s exact test.

Example for a one computation flow



|Anchor Thr|Partner Thr|  |Real CE+|Real Total|Rand CE+|Rand Total|Fold|P-value|
|----------|-----------|--|--------|----------|--------|----------|----|-------|
|1         |1          |  |29      |248       |467     |9672      |2.42|3.0E-05|
|1         |2          |  |25      |202       |370     |7878      |2.64|2.6E-05|
|1         |3          |  |64      |326       |599     |12714     |4.17|4.8E-21|
|1         |4          |  |62      |451       |940     |17589     |2.57|6.1E-11|
|1         |5          |  |51      |476       |945     |18564     |2.10|1.7E-06|
|2         |1          |  |21      |139       |223     |5421      |3.67|9.0E-07|
|2         |2          |  |12      |119       |207     |4641      |2.26|1.7E-02|
|2         |3          |  |24      |173       |290     |6747      |3.23|1.4E-06|
|2         |4          |  |21      |268       |508     |10452     |1.61|4.3E-02|
|2         |5          |  |17      |309       |584     |12051     |1.14|0.68   |
|3         |1          |  |15      |171       |321     |6669      |1.82|4.0E-02|
|3         |2          |  |7       |139       |270     |5421      |1.01|1.00   |
|3         |3          |  |15      |234       |483     |9126      |1.21|0.54   |
|3         |4          |  |19      |330       |664     |12870     |1.12|0.70   |
|3         |5          |  |20      |353       |670     |13767     |1.16|0.53   |
|4         |1          |  |20      |208       |420     |8112      |1.86|1.4E-02|
|4         |2          |  |9       |186       |376     |7254      |0.93|1      |
|4         |3          |  |12      |283       |569     |11037     |0.82|1      |
|4         |4          |  |22      |379       |776     |14781     |1.11|0.72   |
|4         |5          |  |20      |441       |847     |17199     |0.92|1      |
|5         |1          |  |19      |220       |418     |8580      |1.77|2.4E-02|
|5         |2          |  |10      |172       |310     |6708      |1.26|0.56   |
|5         |3          |  |13      |280       |485     |10920     |1.05|1.00   |
|5         |4          |  |25      |379       |804     |14781     |1.21|0.36   |
|5         |5          |  |20      |430       |829     |16770     |0.94|1      |



* __Files <*.best>, the list of predicted CEs__. For each recognized CE MCOT provides (1) the header of a peak, (2) the start and the end positions of each motif in a peak, (3) mutual location (Full / Partial / Spacer types and the respective base pair measures, see above), (4) the strands of the motifs in a peak and the mutual orientation of the motifs (one of four types), (5) PWM scores and DNA sequences of the motifs.



|\#Seq    |Anchor start|Anchor end|Partner start|Partner end|Mutual Loc|Loc Type|Strands|Mutual Ori|Anchor score|Partner score|Anchor seq  |Partner seq  | 
|---------|------------|----------|-------------|-----------|----------|--------|-------|----------|------------|-------------|------------|-------------|
|Seq 491  |135         |146       |134          |146        |0F        |Full    |+-     |Evert     |0.918       |0.937        |ccaggatgtcaa|ttgacatcctggg|
|Seq 493  |354         |365       |378          |390        |12S       |Spacer  |+-     |Invert    |0.936       |0.955        |tgaggaagtgaa|ttgacattcttcc|
|Seq 513  |156         |167       |148          |160        |5P        |Partial |--     |DirectAP  |0.949       |0.959        |agaggaaatgac|atgacagattggg|


## References

Heinz,S., Benner,C., Spann,N., Bertolino,E., Lin,Y.C., Laslo,P., Cheng,J.X., Murre,C., Singh,H. and Glass,C.K. (2010) Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Mol Cell, 38, 576-589.

Kulakovskiy,I.V., Boeva,V.A., Favorov,A.V., Makeev,V.J. (2010) Deep and wide digging for binding motifs in ChIP-Seq data. Bioinformatics, 26, 2622-2623.

Kulakovskiy,I.V., Vorontsov,I.E., Yevshin,I.S., Sharipov,R.N., Fedorova,A.D., Rumynskiy,E.I., Medvedeva,Y.A., Magana-Mora,A., Bajic,V.B., Papatsenko,D.A., et al. (2018) HOCOMOCO: expansion and enhancement of the collection of transcription factor binding sites models. Nucleic Acids Res., 46, D252-D259.

O'Malley,R.C., Huang,S.C., Song,L., Lewsey,M.G., Bartlett,A., Nery,J.R., Galli,M., Gallavotti,A., Ecker,R. (2016) Cistrome and epicistrome features shape the regulatory DNA landscape. Cell, 165, 1280-1292.




