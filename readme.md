# Readme for MCOT software package

## General description

MCOT (Motifs Co-Occurrence Tool) is a software package for recognition of composite elements (CEs) in a single ChIP-seq dataset (Levitsky et al., [2019](https://doi.org/10.1093/nar/gkz800), [2020](https://doi.org/10.3390/ijms21176023), [2023](https://doi.org/10.20944/preprints202311.1617.v1)). CEs detected by MCOT include two potential binding sites of transcription factors (TFs) in all possible mutual orientations. MCOT considers CEs with a full/partial overlap of motifs or with a spacer in a certain range. Each potential CE recognized by MCOT contains the motif of immunoprecipitated TF in respective ChIP-seq experiment (anchor motif) and another motif (partner). Identical/distinct anchor and partner motifs imply the search for CEs of homotypic or eterotypic type (respectively). The web service [WebMCOT](https://webmcot.sysbio.cytogen.ru/) [Levitsky et al., 2022](https://doi.org/10.3390/ijms23168981) implements MCOT software package.

## Implementation

MCOT implemented in C++ and it can be conventionally compiled in Linux or Windows operating system. To run MCOT user should compile the corresponding source code file. Files mcot_anchor.cpp and mcot.cpp respect to one-partner and many partners options for Position Weight Matrix (PWM) model of a binding site. File anchor_pro.cpp respects to one-partner option, but it runs with arbitrary models of site, including not-PWM ones, e.g. [BaMM](https://github.com/soedinglab/BaMM_webserver) [(Siebert and Söding, 2016)](https://doi.org/10.1093/nar/gkw521), and [SiteGA](https://github.com/parthian-sterlet/sitega) [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545)

## Install
[Use <>code button](https://docs.github.com/en/repositories/working-with-files/using-files/downloading-source-code-archives) to dowload this repository

(Linux) 
```
cd <project>
cd run
chmod a+x unpack.sh
./unpack.sh
cd ..
```
This unpacks required input data.

Run in terminal (Packages “build-essential” and “cmake” 
should be installed on Ubuntu system):
```
mkdir cmake
cd cmake
cmake ..
make
```

All executable files for one-partner, many-partners and anchor-pro options will be in src/anchor\_vs\_one, src/anchor\_vs\_many and src/anchor\_pro

(Windows) 
```
cd <project>
cd run
unpack.bat
cd ..
```
This unpacks required input data.

Run in terminal (Win -> Visual Studio 2017 -> Visual Studio Tools -> 
VC -> Native Tools x64. Else, you should install “CMake” module while VS 2017 installing)

```
mkdir cmake
cd cmake
cmake ..
MSBuild mcot-kernel.sln /p:Configuration=Release /p:Platform=Win32
```
Programs `anchor_vs_many`, `anchor_vs_one` and `anchor_pro` for one\_partner, many\_partners and arbitrary\_models\_one\_partner
options should be in `src/anchor_vs_one_partner/Release/`,  `src/anchor_vs_many_partners/Release/` and  `src/anchor_pro/Release/`


## Command line arguments

The command line for one-partner option:


`./mcot_anchor.exe <1 fasta> <2 anchor.motif> <3 partner.motif> <4 minimal spacer length> <5 maximal spacer length> <6 file of whole-genome promoters with its path> <7 pvalue_thr> <8 -log10[p-value]_thr> <9 asymmetry_fold(-log10(ERR))>`


The command line for many-partner option:


`./mcot.exe <1 fasta> <2 anchor.motif> <3 partners.library> <4 minimal spacer length> <5 maximal spacer length> <6 path to whole-genome promoters with its path> <7 pvalue_thr> <8 -log10[p-value]_thr> <9 asymmetry_ratio(-log10(ERR))>`


`<1 fasta>` = DNA sequences of peaks in FASTA format, a minimum recommended number of peaks is about 300-500, the maximum number is not restricted, however 5000-10000 or higher number of peaks requires a higher computation time than several thousands of peaks, hence about 1000-2000 peaks are enough. Sequences should have lengths substatially higher than lengths of recognition models for anchor and partner motifs to contain possible composite elememnts with an overlap or spacer.


`<2 anchor.motif>`, `<3 partner.motif>` = frequency matrices of the first and second motifs motifs in [the standard format](https://github.com/parthian-sterlet/mcot-kernel/blob/master/examples/one/jun.motif), e.g.  

\***

```
> Motif name
0.001    0.001    0.942    0.056
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

`<3 partners.library>` = binary file od library of partner motifs. The file name must contain one of the following substrings: "h12core_hg38" , "h12core_mm10" , "h11core_hg38" , "h11core_mm10" , "dapseq", "jaspar24_at10", "jaspar24_dm6". These values respect to the accepted for analysis libraries of motifs for human/murine TFs derived from _in vivo_ and _in vitro_ data from [Hocomoco version 12](http://hocomoco12.autosome.ru/) (1420/1142) and [Hocomoco version 11](https://hocomoco11.autosome.org/) (401/346) ([Vorontsov et al., 2024](https://doi.org/10.1093/nar/gkad1077), [Kulakovskiy et al., 2018](https://doi.org/10.1093/nar/gkx1106)); the library of 510 motifs derived from DAP-seq data [Plant Cistrome](http://neomorph.salk.edu/dap_web/pages/index.php) for *A.thaliana* TFs ([O’Malley et al., 2016](https://doi.org/10.1016/j.cell.2016.08.063)), and two libraries of 555/151 motifs for plants/insects from JASPAR2024 [Rauluseviciute et al., 2024](https://doi.org/10.1093/nar/gkad1059). See lists of partner motifs accepted for analysis for these libraries [here](https://github.com/parthian-sterlet/mcot-kernel/blob/master/include/pfm_list.h) For each partner motif the position frequency/weight matrices (PFM/PWM), and the table "Thresholds vs ERRs for various ERRs are stored in a binary format, so that each library represents one file, see [partners](https://github.com/parthian-sterlet/mcot-kernel/blob/master/partners/) folder, e.g. [h11core_hg38.binary](https://github.com/parthian-sterlet/mcot-kernel/blob/master/partners/h11core_hg38.binary). The lists of matrices for all libraries (PFM/PWM, position frequency/weight matrices), lists of thresholds and ERRs for all matrices, and notations about whether each matrix is taken or rejected by the [MCOT filter](https://doi.org/10.1093/nar/gkz800) are in [summary XLSX file](https://github.com/parthian-sterlet/mcot-kernel/blob/master/partners/MCOT_dec2023_motifs_libraries.xlsx). For each library, this file lists all partner motifs, names of corresponding TFs, TF classes and TF families. For each motif the minimal ERR value is computed, its [threshold values of 2E-5](https://doi.org/10.1093/nar/gkz800) mark partner motifs accepted to the analysis. The lists of PFMs from [Hocomoco12](https://hocomoco12.autosome.org/)/[Hocomoco11](https://hocomoco11.autosome.org/), [Plant Cistrome](http://neomorph.salk.edu/dap_web/pages/index.php), [JASPAR2024](https://doi.org/10.1093/nar/gkad1059) and corresponding PWMs are [here](https://github.com/parthian-sterlet/mcot-kernel/blob/master/partners/MCOT_dec2023_pcm_pwm.tar.gz). See below [detailed description and source code](https://github.com/parthian-sterlet/mcot-kernel/blob/master/readme.md#generation-of-partner-library) for generation of binary files according to input partner PFMs.

`<4 minimal spacer length>` = integer value from 0 to \<maximal spacer length>  (the default value 0 is recommended, any positive value restricts short spacers)

`<5 maximal spacer length>` = integer value from 0 to 100 (the default value 30)

`<6 file of whole-genome promoters with its path>` =  fasta file of whole-genome dataset of promoters. Four files in folders “hs”, “mm”, “at” and “dm” imply application of promoters of protein coding genes for *H. sapiens*, *M. musculus*, *A. thaliana* and *D. melanogater* for setting of thresholds for input motifs.

`<7 pvalue_thr>` = recognition threshold of motifs, Expected Recognition Rate (ERR), ERR is computed as a recognition rate for the whole-genome set of promoters of protein-coding genes, default value 0.0005, the maximal allowable value is 0.0025.

`<8 -log10[p-value]_thr>` = threshold of -Log10[p-value] to display the significances of enrichment of CEs in output data, p-values are computed fy Fisher exact test (the default value 10)

`<9 asymmetry_ratio(-log10(ERR))>` = the threshold for asymmetry ratio (TAR) restricting relationship of ERR values of two motifs in CEs, e.g. for two mofs in CE values ERR1 and ERR2 means the assymetry ratio (AR) is equal to Max\{-log10(ERR1), -log10(ERR2)\} / Min\{-log10(ERR1), -log10(ERR2)\}, and criteria AR  > TAR and AR <= TAR mean asymmetric and symmetric CEs, respectively.

The command line for anchor_pro option:

`./mcor_anchor_pro.exe <1 file_fasta> <2 motif_model1.profile> <3 motif_model2.profile> <4 int motif_model1.length> <5 int motif_model2.length> <6 int motif_model1.table_thr_err> <7 int motif_model2.table_thr_err> <8 int spacer_min> <9 int spacer_max> <10double pvalue_thr> <11double -log10[p-value]_thr> <12double asymmetry_ratio(-log10(ERR))>`

`<1 file_fasta>` = DNA sequences of tested peaks in FASTA format

`<2 motif1.profile>` = Profile for the first model, see [example profile of model 1](https://github.com/parthian-sterlet/mcot-kernel/blob/master/examples/pro/creb1_49_pwm)

`<3 motif2.profile>` = Profile for the second model, see [example profile of model 2](https://github.com/parthian-sterlet/mcot-kernel/blob/master/examples/pro/creb1_49_sga)

`<4 int motif_model1.length>` = integer value, length of the first model

`<5 int motif_model2.length>` = integer value, length of the second model

`<6 int motif_model1.table_thr_err>` = Table **Threshold vs. -Log10(ERR)** for the first motif, see [example of distribution for model 1](https://github.com/parthian-sterlet/mcot-kernel/blob/master/examples/pro/GSM2827249_CREB1_hg38_pwm.dist)

`<7 int motif_model2.table_thr_err>` = Table **Threshold vs. -Log10(ERR)** for the second motif, see [example of distribution for model 2](https://github.com/parthian-sterlet/mcot-kernel/blob/master/examples/pro/GSM2827249_CREB1_hg38_sga.dist)

`<8 int spacer_min>` = integer value from 0 to \<maximal spacer length>  (the default value 0 is recommended, any positive value restricts short spacers)

`<9 int spacer_max>` = integer value from 0 to 100 (the default value 29)

`<10 pvalue_thr>` = recognition threshold of motifs transformed to the logarithmic scale, -log10(ERR), default value 0.0005

`<11 -log10[p-value]_thr>` = threshold to display the significances of enrichment of CEs in output data (the default value 10)

`<12 asymmetry_ratio(-log10(ERR))>` = the ratio ratio_thr restricting ERR values of two motifs in asymmetric CEs, e.g. for two value ERR1 and ERR2 for certain CE means ratio_thr = Max\{-log10(ERR1), -log10(ERR2)\} / Min\{-log10(ERR1), -log10(ERR2)\}, and ratio > ratio_thr and ratio < ratio_thr mean asymmetric and symmetric CEs, respectively.

## Input data

MCOT requires (a) DNA sequences of ChIP-seq peaks and (b) anchor and partner motifs. We recommend application of a conventional *de novo* motif search tool, e.g. [Homer](http://homer.ucsd.edu/homer/) ([Heinz et al., 2010](https://doi.org/10.1016/j.molcel.2010.05.004)) and [STREME](https://meme-suite.org/meme/tools/streme) ([Bailey et al., 2021](https://doi.org/10.1093/bioinformatics/btab203)) to define an anchor motif. 

MCOT have two options for definition of the partner motif:

* a matrix of partner motif (one-partner option);

* a library of known motifs (many partners option).

`<anchor_pro>` requires input files [**Table Threshold vs. ERR**](https://github.com/parthian-sterlet/mcot-kernel/blob/master/examples/pro/GSM2827249_CREB1_hg38_pwm.dist) for both models. This file contains the list of pairs {Threshold, -Log10(ERR)} values.
For a PWM model the respictive file can be taken as the output files of runs with `<anchor_vs_one>` or `<anchor_vs_many>` options, respecting to the anchor motif <err\*\.txt>. For a non-PWM model, the corresponding table should be deduced from the recognition profile of potential hits for the whole genome dataset of promoters of protein-coding genes, e.g. the [SiteGA](https://github.com/parthian-sterlet/sitega) [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545) tool has a special option to compute the required table

Advanced options include: 
* the minimal and maximal spacer lengths restricting the search of spaced locations of motifs The default range of spacer length is [0; 29] nt. MCOT allows the variation of the upper limit of spacer length up to 100 base pairs;
* Expected Recognition Rate (ERR), the frequency of motifs occurrence is required to set the recognition thresholds for motifs, default value 5E-4, the maximal allowable value is 2.5E-3;
* CE enrichment, threshold –Log10(p-value(CE)) for CE enrichment, default value 5;
* Asymmetry ratio, for two value ERR1 and ERR2 for certain CE means Max{-log10(ERR1), -log10(ERR2)} / Min{-log10(ERR1), -log10(ERR2)}

## Motifs recognition

MCOT with options one\_partner, many\_partners applies the recognition model of PWM for mapping motifs in peaks, otherwise for option anchor\_pro MCOT takes ready mapping of predicted hits from a file. 
For each model, MCOT uses five thresholds {T[1],...T[5]} according to the unified set of ERRs for a whole-genome dataset of promoters, e.g. for the option of command line `<pvalue_thr>` = 5E-4, five ERRs are equal to {5.24E-5, 1.02E-04, 1.9E-4, 3.33E-4, 5E-4}. The profile of the most stringent hits contains matrix scores T ≥ T[1], the next profile comprises PWM scores {T} in the range T[2] ≥ T > T[1], etc. Hence, MCOT computes five profiles of hits with certain level of conservation for each input motif. Note that the change of the recognition threshold `<pvalue_thr>` defines the most mild (the fifth) threshold and proportionally shifts the rest four thresholds. The next thresholds (ERR[3] ... ERR[0]) are defined as follows: k[4]= 1.5, {ERR[i]=ERR[i+1]/k[i+1], k[i] = 1 + k[i+1]/2}. If ERR[4] = 5E-4 than {ERR[3], ERR[2], ERR[1], ERR[0]} = {3.33E-4, 1.9E-4, 1.02E-04, 5.24E-5}. If ERR[4] = 1E-3 than {ERR[3], ERR[2], ERR[1], ERR[0]} = {6.67E-4, 3.81E-4, 2.03E-04, 1.05E-4}.

## Composite elements search and annotation

MCOT classifies CEs structure according to the following attributes:
* _Homotypic_ and _Heterotypic_ CEs respect the same or distinct motif models. This means very CEs of binding sites of structurally similar or differnt transcription factors.
* _Orientation_. Four types of distinct mutual orientations of two sites in CEs are considered: in the same DNA strand (Direct Anchor/Partner and Direct Partner/Anchor), in opposite strands (Everted and Inverted);
* _Overlap or Spacer_. There are three distinct cases of mutual locations of two sites within a CE: Full overlap (one motif located entirely within another one); Partial overlap; and Spacer. To describe each case MCOT uses the following characteristics: the distance between nearest borders of two motifs (Full); the length of overlapped region (Partial); and the length of spacer;
* _Asymmetry of сonservation_. All predicted CEs are subdivided into two types: those with more notable or less notable difference in conservation of two motifs, these types are called asymmetric and symmetric. Heterotypic asymmetric CEs are further subdivided into those with more conserved anchor and partner motifs. The conservation of a site is estimated as -Log10(ERR), where ERR is computed by the respective score of recognition model. See principal explanation of asymmetry within CE in [Levitsky et al. (2023)](https://doi.org/10.20944/preprints202311.1617.v1)

For each of 25 combinations of motifs conservation and each computation flow MCOT compiles the 2x2 contingency table (Table 1) and compute the significance of Fisher’s exact test p-value(CE) that compares the content of CEs in ChIP-seq peaks with that for background model.

**Table 1.** 2x2 contingency table for calculation of CE significance.

|                             | Sequences without CEs               | Sequences with CEs  | Total, sequences with both motifs|
|-----------------------------|-------------------------------------|---------------------|----------------------------------|
|Observed (peaks)             |$`Obs_{CE-} = Obs_{Tot} - Obs_{CE+}`$|$`Obs_{CE+}`$        | $`Obs_{Tot}`$                    |
|Expected (permuted sequences)|$`Exp_{CE-} = Exp_{Tot} - Exp_{CE+}`$|$`Exp_{CE+}`$        | $`Exp_{Tot}`$                    |

The background model implies the preservation of content for each motif in each peak. MCOT generates background profiles of hits of anchor and partner motifs iteratively for each peak with a special permutation procedure.

MCOT computes CE significance separately for five computation flows:  Any (Spacer or Overlap), Full, Partial, Overlap (Full and Partial), Spacer. 
Fisher’s exact test computes the CE enrichment for each of 5x5 combinations of conservation of motifs. 

Finally, MCOT estimates the similarity of anchor and partner motifs with the motifs similarity p-value. MCOT used matrix column similarities measures PCC ([Pearson Correlation Coefficient, Pietrokovski, 1996](https://academic.oup.com/nar/article/24/19/3836/2384639)) and SSD ([Sum of Squared Distances, Sandelin and Wasserman, 2004](https://doi.org/10.1016/j.jmb.2004.02.048)) to compute two p-values. If maximum among them was less than 0.05, then pair of motifs had the significant similarity. Hence, a respective CE may be a possible false positive prediction.

Additionally, MCOT computed significance of CEs with more conserved anchor motif and more conserved partner motif. These calculations are performed according to scheme represented above in Table 1, but MCOT counts either CEs with more conserved Anchor or Partner motifs, i.e. either -Log10[ERR(Anchor)] > -Log10[ERR(Partner)] or -Log10[ERR(Anchor)] ≤ -Log10[ERR(Partner)].

Finally, to estimate the asymmetry of CE, MCOT partitions all CEs on those with more conserved Anchor or Partner motifs and compute the significance of asymmetry with Fisher’s exact test (Table 2). The asymmetry significance -Log10[P-value] is printed as positive (with sign ‘+’) in the case of enrichment toward the Anchor motif, otherwise, MCOT sign ‘-’ (negative value) denotes the enrichment toward the Partner motif.  These calculation yielded for each computation flow one p-value(CE, Asymmetry).

**Table 2.** 2x2 contingency table for calculation of CE asymmetry.

|                                    | CEs with more conserved anchor motif| CEs with more conserved partner motif  | CEs, total                              |
|------------------------------------|-------------------------------------|----------------------------------------|-----------------------------------------|
|Observed (CEs in peaks)             |$`Obs_{CE, Anchor}`$                 |$`Obs_{CE, Partner}`$                   | $`Obs_{CE, Anchor} + Obs_{CE, Partner}`$|
|Expected (CEs in permuted sequences)|$`Obs_{CE, Anchor}`$                 |$`Obs_{CE, Partner}`$                   | $`Obs_{CE, Anchor} + Obs_{CE, Partner}`$|


Detailed enrichment or depletion of CEs with specific combinations of motifs conservation are represented in the scatterplot text file plot_*, see below.

To take into account multiple comparisons we applied the Bonferroni’s correction and used the following critical values to filter out not significant results:
* significance of CEs regardless motifs conservation, Bonferroni_CE = 0.05/(Nfor\*Nback\*Nflow\*Nthr\*Nthr);
* significance of asymmetric CEs toward one of motifs, Bonferroni_CE(AncPar) = 0.05/(Nfor\*Nback\*Nflow\*2);
* CE asymmetry, Bonferroni_Asym = 0.05/(Nfor\*Nback\*Nflow). 
Here Nfor and Nback means the size of foreground and background datasets (i.e., the number of peaks and random sequences, which generated in MCOT, 
Nflow = 5 designates the number of MCOT computation flows and Nthr = 5 means the number of thresholds for each motif. 

## Output data

MCOT gives the following output data:


* __Files <\*\_thr5>, recognition profiles of motifs__ . Each file respects to one motif. A file has fasta-like format, i.e. for each peak the header line starts with ‘>’ symbol. Next, each subsequent line represents one hit in a peak, particularly it position, respective conservation value -Log10(ERR) and DNA strand.

```
Example
>Seq 1  Thr 0.864497    Nsites 1
205 3.304404    -
>Seq 2    Thr 0.864497    Nsites 0
>Seq 3    Thr 0.864497    Nsites 2
88    3.607778    -
160    3.338550    -
```

Here and below ChIP-seq data for mouse FoxA2 and CE FoxA2-HNF1β ([Wederell et al., 2008](https://doi.org/10.1093/nar/gkn382)) illustrate MCOT application. The anchor FoxA2 motif we deduced from de novo search [Homer](http://homer.ucsd.edu/homer/) ([Heinz et al., 2010](https://doi.org/10.1016/j.molcel.2010.05.004)) and the partner HNF1β motif we took from the mouse Hocomoco core collection ([Kulakovskiy et al., 2018](https://doi.org/10.1093/nar/gkx1106)).


* __Files <err\*\.txt>, Table Threshold vs. ERR (Expected Recognition Rate)__ . File contains two columns: threshold and ERR estimated as the site density for the whole genome dataset of protein-coding genes

```
Example
0.99745670	6.786437089
0.99660894	6.724287943
0.99612450	6.395229799
...	
0.81195753	1.746523876
0.81191716	1.746252613
0.81187679	1.745894418
```

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

Example:


|# Motif  | Motif Name | Full  overlap, -Log10[P-value] | Partial overlap,-Log10[P-value]| Overlap, -Log10[P-value] | Spacer, -Log10[P-value] | Any, -Log10[P-value] | Similarity to Anchor, -Log10[P-value] | Similarity to Anchor, SSD | Similarity to Anchor, PCC |...|
|---------|------------|--------------------------------|--------------------------------|--------------------------|-------------------------|----------------------|---------------------------------------|---------------------------|---------------------------|---|
|Anchor| FOXA2      | 0.0                            | 97.4                           | 92.4                     | 9.8                     | 35.4                 | n/a                                   | n/a                       | n/a                       |   |
|Partner| HNF1B      | 76.6                           | 53.2                           | 64.3                     | 4.3                     | 35.3                 | 0                                     | 0                         | 0                         |   |




The second block of data consequently for five computation flows represents significances for CEs with more conserved Anchor and Partner motifs.

Example:

| # Motif | Motif Name | ... | Full overlap, Conservative Anchor, -Log10[P-value] | Full overlap, Conservative Partner, -Log10[P-value] | Partial overlap, Conservative Anchor, -Log10[P-value] | Partial overlap, Conservative Partner, -Log10[P-value] | ... |
|---------|------------|-----|----------------------------------------------------|-----------------------------------------------------|-------------------------------------------------------|--------------------------------------------------------|-----|
| Anchor| FOXA2      |     | n/a                                                | n/a                                                 | n/a                                                   | n/a                                                    |     |
|Partner| HNF1B      |     | 50.8                                               | 0.6                                                 | 40.1                                                  | 19.2                                                   |     |




The third block of data consequently for five computation flows represents five significances of asymmetry (p-value) in CEs: ‘Anchor vs. Partner’. The positive/negative values reflect the enrichment toward the Anchor/Partner motifs.

Example:


| # Motif  | Motif Name | ... | Full overlap, Asymmetry to Anchor+/Partner-, -Log10[P-value] | Partial overlap, Asymmetry to Anchor+/Partner-, -Log10[P-value] | Overlap, Asymmetry to Anchor+/Partner-, -Log10[P-value] | ... |
|----------|------------|-----|--------------------------------------------------------------|-----------------------------------------------------------------|---------------------------------------------------------|-----|
| Anchor| FOXA2      |     | n/a                                                          | n/a                                                             | n/a                                                     |     |
| Partner| HNF1B      |     | -46.9                                                        | -73.1                                                           | -113.1                                                  |     |

The final fourth block shows recommended Bonferroni’s correction thresholds (-Log10[P-value]) for significances:

| # Motif   | ... | Bonferroni_CE | Bonferroni_CE(AncPar) | Bonferroni_Asym |
|-----------|---|---------|-------|--------------------------------------------------------------|
| Anchor  || 11.8       |       |       |
| Partner || 11.71      | 10.61 | 10.31 |


* __File <out_hist\*>, the abundance of various CE types as a function of mutual orientation and location of the motifs__ The percentage of peaks containing CE variants specific in mutual orientation (four types) and mutual locations from a few possible full overlaps (‘F’), through a variety of partial overlaps (‘P’) and finally from the minimal to the  spacer length (‘S’).

Example:

| |1F|0F|11P|10P|9P|8P|7P|6P|5P|4P|3P|2P|1P|0S|1S|
|-|--|--|---|---|--|--|--|--|--|--|--|--|--|--|--|
|Everted||2.54||||3.88|||0.16|0.94||0.07|0.18|0.76|0.09|
|Inverted||0.07|2.51|||0.07|1.08||0.04|0.13|0.72|0.16|0.18|0.31|0.67|
|DirectPA||0.04|||1.86||||0.09|0.45||0.02|0.11|0.13|0.09|
|DirectAP||0.04|0.18|||0.07|0.16|0.11||0.02|0.18|0.11|0.11|0.2|0.22|
|Any||6.84|6.84|0.00|4.77|10.23|3.10|0.29|0.75|3.79|2.24|0.92|1.49|3.51|2.70|
|Cumulative||6.84|9.60|9.60|14.02|19.43|19.83|20.06|20.34|20.92|21.21|21.95|22.87|23.79|24.66|

This example shows distribution for heterotypic CE, since notations DirectAP / DirectPA imply Anchor-Partner / Partner-Anchor cases. A homotypic CE has only one direct orientation of CE, Anchor-Anchor (notation DirectAA). 'Any' row implies frequency of CEs with a certain spacer or an overlap for all four orientations. 'Cumulative' row means the sum of frequencies for all mutual locations which respect to the current and closer positioning of motifs, i.e. the sum of values of 'All' row from the most left column to the current one.


* __Files <fisher\_*>, the 2x2 tables of CE significance for five computation flows for all motifs pairs and for  all 5x5 combination of motifs conservation__. 
Each line of output file contains data concerning one 2x2 contingency table, in particular (1) the designation of conservation (PWM threshold (T), indices from 
1 to 5 mean the change from the most stringent to the most permissive, see above); (2) four counts for 2x2 contingency table (see Table 1 above), 
‘the number of peaks containing at least one CE (CE+) & ‘the number of peaks containing at least one hit of each motif (Total)’ for peaks (Real) 
and permuted (Rand) datasets. Finally, the table contains significance of CEs (p-values) computed by Fisher’s exact test for 25 cells of 5x5 tables 
of combinations of thresholds. Next, the respective data are shown for (a) significances of CEs with more conserved Anchor and Partner motifs (lines 'Anchor', Partner'), significances of any asymmetric CEs with more conserved either Anchor or Partner motifs (line 'Asymmetry'), significances of symmetric CEs with the same  conservation of Anchor and Partner motifs (line 'Symmetry'), and (b) significances of asymmetry in CEs ‘Anchor vs. Partner’ with the positive/negative Fold respecting to the enrichment toward the Anchor/Partner motifs, significances of asymmetry in CEs ‘Asymmetry vs. Symmetry with the positive/negative Fold respecting to the enrichment toward the Asymmetry/Symmetry in the motifs conservation. In these calculation points (a) and (b) imply counting of peaks and CEs, respectively. 

Example below shows FOXA2 (Anchor) and HNF1B (Partner) motifs for Overlap computation flow.

Anchor Thr | Partner Thr |  | Real CE+ | Real Total | Rand CE+ | Rand Total | Fold | P-value
|----------|-------------|--|----------|------------|----------|------------|------|-------|
A 1 | P 1 |  | 13 | 211 | 20 | 14137 | 43.550 | 4.59e-16
A 1 | P 2 |  | 10 | 167 | 12 | 11189 | 55.833 | 2.00e-13
A 1 | P 3 |  | 9 | 142 | 18 | 9514 | 33.500 | 9.33e-11
A 1 | P 4 |  | 22 | 312 | 27 | 20904 | 54.593 | 7.98e-28
A 1 | P 5 |  | 32 | 359 | 38 | 24053 | 56.421 | 2.97e-40
A 2 | P 1 |  | 10 | 125 | 18 | 8375 | 37.222 | 3.45e-12
A 2 | P 2 |  | 13 | 78 | 8 | 5226 | 108.875 | 9.84e-20
A 2 | P 3 |  | 10 | 91 | 5 | 6097 | 134.000 | 8.07e-16
A 2 | P 4 |  | 21 | 172 | 23 | 11524 | 61.174 | 1.42e-27
A 2 | P 5 |  | 23 | 191 | 25 | 12797 | 61.640 | 4.14e-30
A 3 | P 1 |  | 7 | 160 | 18 | 10720 | 26.056 | 5.02e-08
A 3 | P 2 |  | 11 | 106 | 17 | 7102 | 43.353 | 7.15e-14
A 3 | P 3 |  | 15 | 107 | 18 | 7169 | 55.833 | 9.87e-20
A 3 | P 4 |  | 18 | 229 | 37 | 15343 | 32.595 | 4.69e-20
A 3 | P 5 |  | 23 | 263 | 32 | 17621 | 48.156 | 3.31e-28
A 4 | P 1 |  | 11 | 192 | 35 | 12864 | 21.057 | 4.45e-11
A 4 | P 2 |  | 11 | 133 | 13 | 8911 | 56.692 | 9.71e-15
A 4 | P 3 |  | 24 | 124 | 29 | 8308 | 55.448 | 5.58e-31
A 4 | P 4 |  | 36 | 303 | 58 | 20301 | 41.586 | 7.52e-42
A 4 | P 5 |  | 37 | 286 | 39 | 19162 | 63.564 | 5.87e-48
A 5 | P 1 |  | 17 | 171 | 13 | 11457 | 87.615 | 3.18e-24
A 5 | P 2 |  | 20 | 124 | 15 | 8308 | 89.333 | 1.23e-28
A 5 | P 3 |  | 80 | 151 | 59 | 10117 | 90.847 | 1.07e-118
A 5 | P 4 |  | 55 | 229 | 37 | 15343 | 99.595 | 7.19e-79
A 5 | P 5 |  | 53 | 250 | 33 | 16750 | 107.606 | 9.56e-77
Anchor |  |  | 170 | 1392 | 237 | 93264 | 48.059 | 1.76e-199
Partner |  |  | 95 | 1057 | 106 | 70819 | 60.047 | 3.61e-118
Asymmetry |  |  | 399 | 1741 | 598 | 116647 | 44.704 | 1e-300
Symmetry |  |  | 0 | 0 | 0 | 0 | 1.000 | 1
Anchor_Partner |  |  | 280 | 1418 | 335 | 904 | 0.533 | 7.08e-20
Asym_Sym |  |  | 1418 | 1418 | 904 | 904 | 1.000 | 2


* __Files <\*.best>, the list of predicted CEs__. 
For each recognized CE MCOT provides (1) the header of a peak, (2) the start and the end positions of each motif in a peak, (3) mutual location (Full / Partial / Spacer types and the respective base pair measures, see above), (4) the strands of the Anchor/Partner motifs in a peak and the mutual orientation of the motifs (one of four types), (5) conservation scores and DNA sequences of the motifs. To provide detailed information on asymmetry ‘Anchor vs. Partner’ in comparison of peaks and permuted sequences MCOT provides files <real*.best> and rand*.best> that show the lists of predicted CE for peaks and permuted sequences.


|#Seq|Anchor start|Anchor end|Partner start|Partner end|Mutual Loc|Loc Type|Strands|Mutual Ori|Anchor hit conservation, -Log10(ERR)|Partner hit conservation, -Log10(ERR)|Anchor seq|Partner seq|
|----|------------|----------|-------------|-----------|----------|--------|-------|----------|------------------------------------|-------------------------------------|----------|-----------|
|Seq 18|271|282|268|282|0F|Full|-+|Inverted|3.365|3.447|tgtttatctttc|agtgaaagataaaca|
|Seq 23|565|576|569|583|8P|Partial|-+|Everted|3.457|4.506|tattgacttacc|agtcaataagttaca|
|Seq 33|145|156|179|193|22S|Spacer|++|DirectAP|3.448|4.289|tgttgacagact|ggttaatgctttcct|



* __Files <plot\_\*>, heatmaps that show the CE asymmetry, i.e. the abundance of CEs with various ratios of conservation of Anchor and Partner motifs__ 
For each of five computation flows (Full, Partial, Overlap, Spacer and Any) one heatmap is computed. For foreground and background data (peaks and sequences with permuted hits) the respective list of predicted CEs are subdivided on two fractions: those with more conserved Anchor and Partner motifs. The conservation of a hit is estimated with the respective -Log10(ERR) value. The minimal conservation value is equal to -Log10(5E-4) ~ 3.3. Next, MCOT computes two matrices {OBSi,j} and {EXPi,j} of absolute frequencies of conservation of Anchor and Partner motifs for observed and expected data. Here indices i and j imply the conservation -Log10(ERR) of Anchor and Partner motifs. E.g., for ERR = 5E-4 (-Log10(5E-4) ~ 3.30) this conservation is falling within the ranges [<3.5], [3.5..3.7], [3.7..3.9] etc. up to [5.3..5.5] and [>5.5]. Finally, the per mille measure transforms the absolute frequencies to relative ones as follow:

$`{1000*OBS_{i,j}/N(OBS)}`$ and $`{1000*EXP_{i,j}/N(EXP)}`$, 

where $`N(OBS)`$ and $`N(EXP)`$ are total counts of predicted CEs in observed and expected lists. The output heatmap shows the difference between relative frequencies of observed and expected CEs. Example below shows the asymmetry toward the Partner (HNF1B) motif in Overlap computation flow.


| |<3.5|3.5..3.7|3.7..3.9|3.9..4.1|4.1..4.3|4.3..4.5|4.5..4.7|4.7..4.9|4.9..5.1|5.1..5.3|5.3..5.5|>5.5|
|-|----|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|----|
|<3.5|14|14|-8|-14|-26|-12|1||||||
|3.5..3.7|7||-15||-7|-7||||||-2|
|3.7..3.9|37|-4|-2|-6|4|-11|-3|||||-3|
|3.9..4.1|17|-1|14|1|-4|-9|||||||
|4.1..4.3|39||-2|8|5|-8|-2||||||
|4.3..4.5||3|-6|4|5|5|||||||
|4.5..4.7||-10|22|||-5|||||||
|4.7..4.9|13|-4|||-3||||||||
|4.9..5.1|-23||-3||5|-4|1||||||
|5.1..5.3|1||||-2||||||||
|5.3..5.5|-4|-3|-2||-3||||||||
|>5.5|4|2|-1|-2|-3||||||||

The same calculations are performed for Anchor-Anchor CEs, in this case the enrichment of symmetric vs. asymmetric CEs are tested.

## Generation of partner library
MCOT use several partner libraries including hundreds of motifs for several specific species. Each library includes the list of motifs, for each motifs its PFM and PWM are provided, and the table of all recognition threshold values of this PWM and respective -Log10(ERR) values are provides too. A file of library is written in binary format. The script prepares this library using the list of PFMs/PWMs and the promoters of all protein coding genes from whole-genome for a specific species.
* The program [pfm_to_pwm_mat.cpp](https://github.com/parthian-sterlet/mcot-kernel/blob/master/src/pfm_to_pwm/pfm_to_pwm_mat.cpp) converts PFM file to PWM file, two formats of PFM file are acceptable, the [default](https://github.com/parthian-sterlet/mcot-kernel/blob/master/examples/library/partner1.pfm), and [transposed](https://github.com/parthian-sterlet/mcot-kernel/blob/master/examples/library/partner2.pfm) matrices.

**Command line arguments:**

`<1 input PFM>` = input PFM format file

`<2 input PWM>` = input PWM format file

* The program [pwm_iz_pwm_thr_dist0.cpp](https://github.com/parthian-sterlet/mcot-kernel/blob/master/src/pwm_thr_err/pwm_iz_pwm_thr_dist0.cpp) takes PFM file, PWM file and FASTA file of promoter to compute a table of PWM thresholds and corresponding -Log10(ERR) values estimating Expected Recognition Rates for PWM. Sequential application of this program to the list of motifs creates the binary file required for the [many-partner option program](https://github.com/parthian-sterlet/mcot-kernel/blob/master/src/anchor_vs_many_partners/mcot.cpp). [command_line_library](https://github.com/parthian-sterlet/mcot-kernel/blob/master/run/command_line_library) shows an example run to construct a library

**Command line arguments:**

`<1 input PFM>` = input PFM format file

`<2 input PWM>` = input PWM format file

`<3 input file of whole-genome promoters with its path>` =  fasta file of whole-genome dataset of promoters, three files in folders “hs”, “mm”, “at” and “dm” imply application of *H. sapiens*, *M. musculus*, *A. thaliana* and *D. melanogater* promoters for setting of thresholds for input motifs.

`<4 output text file pwm.table_thr_err>` = Table **Threshold vs. -Log10(ERR)** for pwm in text format, see [example of distribution](https://github.com/parthian-sterlet/mcot-kernel/blob/master/examples/pro/GSM2827249_CREB1_hg38_pwm.dist)

`<5 output binary file>` = concatenated files of PFM, PWM, and Table **Threshold vs. -Log10(ERR)** in binary format, this file contains only motifs passing criterion on the minimal ERR values of the best PWM hit (see labels Good/Bad in the log output file)

`<6 pvalue_thr>` = maximal Expected Recognition Rate (ERR), ERR is computed as a recognition rate for the whole-genome set of promoters of protein-coding genes, default value is 0.0025.

`<7 dpvalue_thr>` = step for (ERR), this paramater is required to avoid too long output tables, this parameter means the rounding accuracy of the ERR values, default value is 5E-7.

`<8 output file log>` = log file reporting the number of distinct thresholds of the PWM, and the minimal ERR value, which respects its best hit.

`<9 double best_hit_thr>` = double value of the threshold for the minimal ERR values of the best PWM hit, according to [Levitsky et al., 2019](https://doi.org/10.1093/nar/gkz800), default value is 2E-5

## References
[Bailey, T.L. (2021) STREME: accurate and versatile sequence motif discovery. Bioinformatics. 37, 2834–40](https://doi.org/10.1093/bioinformatics/btab203)

[Levitsky, V., Zemlyanskaya, E., Oshchepkov, D., Podkolodnaya, O., Ignatieva, E., Grosse, I., et al. (2019) A single ChIP-seq dataset is sufficient for comprehensive analysis of motifs co-occurrence with MCOT package. Nucleic Acids Res. 47, e139.](https://doi.org/10.1093/nar/gkz800)

[Levitsky, V., Oshchepkov, D., Zemlyanskaya, E., Merkulova, T. (2020) Asymmetric conservation within pairs of co-occurred motifs mediates weak direct transcription factor binding in ChIP-seq data. Int J Mol Sci. 21, 6023.](https://doi.org/10.3390/ijms21176023)

[Levitsky, V.G., Mukhin, A.M., Oshchepkov, D.Y., Zemlyanskaya, E.V., Lashin, S.A. (2022) Web-MCOT Server for Motif Co-Occurrence Search in ChIP-Seq Data. Int. J. Mol. Sci., 23, 8981.](https://doi.org/10.3390/ijms23168981)

[Levitsky, V.G., Raditsa, V.V., Tsukanov, A., Mukhin, A.M., Merkulova, T.I. (2023) Asymmetry of Motifs Conservation within Homotypic Composite Elements Differentiates DNA-Binding Domains of Target Transcription Factors in ChIP-Seq Data. Preprints 2023111617.]( https://doi.org/10.20944/preprints202311.1617.v1)

[Heinz, S., Benner, C., Spann, N., Bertolino, E., Lin, Y.C., Laslo, P., et al. (2010) Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Mol Cell, 38, 576-589.](https://doi.org/10.1016/j.molcel.2010.05.004)

[Kulakovskiy, I.V., Vorontsov, I.E., Yevshin, I.S., Sharipov, R.N., Fedorova, A.D., Rumynskiy, E.I., et al. (2018) HOCOMOCO: expansion and enhancement of the collection of transcription factor binding sites models. Nucleic Acids Res., 46, D252-D259.](https://doi.org/10.1093/nar/gkx1106)

[O'Malley, R.C., Huang, S.C., Song, L., Lewsey, M.G., Bartlett, A., Nery, J.R., et al. (2016) Cistrome and epicistrome features shape the regulatory DNA landscape. Cell, 165, 1280-1292.](https://doi.org/10.1016/j.cell.2016.08.063)

[Pietrokovski, S. (1996) Searching databases of conserved sequence regions by aligning protein multiple-alignments. Nucleic Acids Res., 24, 3836-3845.](https://academic.oup.com/nar/article/24/19/3836/2384639)

[Rauluseviciute, I., Riudavets-Puig, R., Blanc-Mathieu, R., Castro-Mondragon, J. A., Ferenc, K., Kumar, V., Lemma, R. B., Lucas, J., Chèneby, J., Baranasic, D., et al., (2024). JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. Nucleic Acids Res, gkad1059](https://doi.org/10.1093/nar/gkad1059)

[Sandelin, A., Wasserman, W.W. (2004) Constrained binding site diversity within families of transcription factors enhances pattern discovery bioinformatics. J Mol Biol., 338, 207-215.](https://doi.org/10.1016/j.jmb.2004.02.048)

[Siebert, M., Söding, J. (2016). Bayesian Markov models consistently outperform PWMs at predicting motifs in nucleotide sequences. Nucleic Acids Res., 44, 6055–69.](https://doi.org/10.1093/nar/gkw521)

[Tsukanov, A. V., Mironova, V. V., Levitsky, V. G. (2022) Motif models proposing independent and interdependent impacts of nucleotides are related to high and low affinity transcription factor binding sites in Arabidopsis. Front Plant Sci., 13, 938545.](https://doi.org/10.3389/fpls.2022.938545)

[Vorontsov, I. E., Eliseeva, I. A., Zinkevich, A., Nikonov, M., Abramov, S., Boytsov, A., Kamenets, V., Kasianova, A., Kolmykov, S., Yevshin, I. S., et al. (2024). HOCOMOCO in 2024: a rebuild of the curated collection of binding models for human and mouse transcription factors. Nucleic Acids Res, gkad1077](https://doi.org/10.1093/nar/gkad1077)

[Wederell, E.D., Bilenky, M., Cullum, R., Thiessen, N., Dagpinar, M., Delaney, A., et al. (2008) Global analysis of in vivo Foxa2-binding sites in mouse adult liver using massively parallel sequencing. Nucleic Acids Res., 36, 4549-4564.](https://doi.org/10.1093/nar/gkn382) 
