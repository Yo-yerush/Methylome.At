# Methylome.At

Methylome.At is a comprehensive, R-based pipeline for *Arabidopsis thaliana* that processes post-alignment WGBS data, identifies DMRs (including non-CG contexts, using [DMRcaller](https://github.com/nrzabet/DMRcaller) package), integrates multiple genomic resources for functional interpretation, and generates extensive visualizations and annotations to advance understanding of plant epigenetic regulation.

### Publication
[soon...](https://yo)

Methylome.At will produce 10 analysis and each analysis contains CG, CHG and CHH contexts:

![Methylome.At Flow](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/pipeline_scheme.png)

* Chloroplast Conversion Rate
* Total Methylation Levels
* Methylation Distribution
* Gene Body and Transposable Elements Meta-plots
* DMRs Identification (using [DMRcaller](https://github.com/nrzabet/DMRcaller))
* DMRs Distribution Mapping
* Genome Annotation for DMRs and DMPs
* Gene Ontology (GO) Analysis
* KEGG Pathway Enrichment Test
* BedGraph Files Generation for DMRs

# <a name="SystemRequirements"></a>System Requirements
#### Using Conda environment
* Linux Environment
* [Conda](https://docs.conda.io/en/latest/miniconda.html) ([download](https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh))
* CPU: No special restrictions, but it is recommended to work with more than 10 cores for improved efficiency.
### using local R environment
* [R](https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html)>=4.3
* [RCurl](https://cran.r-project.org/web/packages/RCurl/index.html)
* [textshaping](https://github.com/r-lib/textshaping)
* R packages (install 'install_R_packages.R' script)
 ```
dplyr
tidyr
ggplot2
lattice
PeakSegDisk
geomtextpath
parallel
BiocManager
DMRcaller
rtracklayer
topGO
KEGGREST
Rgraphviz
org.At.tair.db
GenomicFeatures
plyranges
 ```


# Installation
1. Download the source code.
  ```
git clone https://github.com/Yo-yerush/Methylome.At.git
cd PATH/TO/Methylome.At/
  ```
2. Setup [conda](https://docs.conda.io/en/latest/miniconda.html) environment and install all R packages
* using build-in setup script:
```
./setup_env.sh
```
* or manually:
```
conda create --name Methylome.At_env
conda activate Methylome.At_env
conda install -c conda-forge r-base r-curl r-rcurl r-textshaping
Rscript scripts/install_R_packages.R
```
# Input files
1. CX_report file is post-alignment data the output of [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/).
> CX_report
```
Chr1     3563    +       0       6       CHG     CCG
Chr1     3564    +       5       2       CG      CGA
Chr1     3565    -       2       3       CG      CGG
Chr1     3571    -       0       5       CHH     CAA
Chr1     3577    +       1       5       CHH     CTA
```

2. Samples table file ([example](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/sample_table_example.txt))

****tab-delimited***, no header
```
treatment	PATH/TO/CX_files/
```


3. Annotation and description files

* By default, Methylome.At provides *gene annotation file* in GFF3 format, *transposable elements annotation file* in text file (provide by [TAIR10](https://www.arabidopsis.org/)), and a *description file* integrated from multiple databases.
* users can alternativley use other *annotation files* (.gtf/.gff/.gff3) and *description file* (.csv/.txt)
> description file columns
```
gene_id Symbol	Short_description	Gene_description	Computational_description	AraCyc.Db	AraCyc.Name	gene_model_type	Protein.families	GO.biological.process	GO.cellular.component	GO.molecular.function	note	Derives_from	old_symbols	EC	KEGG_pathway	refseq_id	PMID

```

# Tutorial 

## Run Methylome.At

**Usage:**
```
$ ./Methylome.At.sh --help

Usage: ./Methylome.At.sh --samples_file FILE [options]

Required argument:
  --samples_file                Path to samples file [required]

Optional arguments:
  --minProportionDiff_CG        Minimum proportion difference for CG [default: 0.4]
  --minProportionDiff_CHG       Minimum proportion difference for CHG [default: 0.2]
  --minProportionDiff_CHH       Minimum proportion difference for CHH [default: 0.1]
  --binSize                     DMRs bin size [default: 100]
  --minCytosinesCount           Minimum cytosines count [default: 4]
  --minReadsPerCytosine         Minimum reads per cytosine [default: 6]
  --pValueThreshold             P-value threshold [default: 0.05]
  --n_cores                     Number of cores [default: 10]
  --GO_analysis                 Perform GO analysis [default: TRUE]
  --KEGG_pathways               Perform KEGG pathways analysis [default: TRUE]
  --annotation_file             Genome Annotation file [default: Methylome.At annotations file (TAIR10 based)]
  --description_file            Description file [default: Methylome.At description file]
  --TEs_file                    Transposable Elements file [default: TAIR10 'Transposable Elements' annotations]
  --Methylome_At_path           Path to Methylome.At [default: PATH/TO/Methylome.At]
```

## Run Methylome.At for **metaPlots**

**Usage:**
```
$ ./Methylome.At_metaPlots.sh --help

Usage: ./Methylome.At_metaPlots.sh --samples_file FILE [options]

Required argument:
  --samples_file                Path to samples file [required]

Optional arguments:
  --Genes_n_TEs                 Analyze Genes and TEs metaPlot [logical. default: TRUE]
  --Gene_features               Analyze Gene Features metaPlot [logical. default: TRUE]
  --minReadsPerCytosine         Minimum reads per cytosine [default: 6]
  --metaPlot_random_genes       Number of random genes for metaPlots. 'all' for all the coding-genes and TEs [default: 10000]
  --n_cores                     Number of cores [default: 20]
  --bin_size_features           Bin-size (set only for 'Gene_features' analysis!) [default: 10]
  --annotation_file             Genome Annotation file [default: Methylome.At annotations file (TAIR10 based)]
  --TEs_file                    Transposable Elements file [default: TAIR10 'Transposable Elements' annotations]
  --Methylome_At_path           Path to 'Methylome.At' directory [default: /home/yoyerush/yo/test_111224/Methylome.At]

```

## .log file and run time
An automated log file will be created during the process. ([example](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/log_file.log))

Run time (40 cores): 
> Methylome.At: 8.09 hours

> Methylome.At_metaPlots: 4.51 hours

## Output Figures
*each analysis also output the data as a table

*all the results taken from [paper]()

**Total Methylation Levels**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Whole_Genome.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Heterochromatin_region.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Euchromatin_region.svg)


**Methylation Distribution**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/ChrPlot_CHG.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/ChrPlot_difference_CHG.svg)


**Gene Body and Transposable Elements Meta-plots**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/TEs_CHG_metaPlot.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Genes_CHG_metaPlot.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Genes_features_CHG_metaPlot.svg)

**DMRs direction (Gain or Loss)**

![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/pie_CHH_gainORloss.svg)
![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/ratio.distribution_CHH_gainORloss.svg)

**DMRs Distribution Mapping**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/ChrPlot_DMRs_CHG.svg)

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/DMRs_Density.svg)


**Genome Annotation for DMRs**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/CHG_genom_annotations.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/legend_genom_annotations.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/CHG_TE.vs.ProteinCodingGenes.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/CHG_TE_Super_Family.svg)

**Gene Ontology (GO) Analysis**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/CHG_Genes_GO.svg)

**KEGG Pathway Enrichment Test**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Promoters.KEGG.svg)

## Output Files tree

####  multiple files key:

> CNTX = CG, CHG, CHH

> FEATURE = gene, promoter, CDS, intron, 5'UTR, 3'UTR

> GO_TYPE = BP (biological process), MF (molecular function), CC (cellular component)

#### short:
```
mto1_vs_wt
 ¦--ChrPlot_CX
 ¦--ChrPlot_DMRs       
 ¦--conversion_rate.csv
 ¦--DMRs_CNTX_mto1_vs_wt.csv
 ¦--DMRs_Density_mto1_vs_wt.svg
 ¦--gain_OR_loss
 ¦--genome_annotation
 ¦   ¦--CNTX
 ¦   ¦--CNTX_genom_annotations.svg
 ¦   ¦--CX_annotation
 ¦   °--TEs_addiotionnal_results
 ¦--GO & KEGG analysis
 ¦   ¦--CNTX
 ¦   ¦   °--FEATURE
 ¦   °--CNTX_FEATURE.svg
 ¦--metaPlots
 ¦   ¦--TEs/Genes/Gene_features
 ¦   ¦   ¦--metaPlot_tables/
 ¦   ¦   °--CNTX_metaPlot.svg
 °--methylation_levels
```
#### long:
```
mto1_vs_wt
 ¦--ChrPlot_CX
 ¦   ¦--ChrPlot_CNTX_mto1_vs_wt.svg
 ¦   °--ChrPlot_difference_CNTX_mto1_vs_wt.svg
 ¦--ChrPlot_DMRs
 ¦   °--ChrPlot_DMRs_CNTX_mto1_vs_wt.svg          
 ¦--conversion_rate.csv
 ¦--DMRs_CNTX_mto1_vs_wt.csv
 ¦--DMRs_Density_mto1_vs_wt.svg
 ¦--gain_OR_loss
 ¦   ¦--pie_CNTX_gainORloss.svg
 ¦   °--ratio.distribution_CNTX_gainORloss.svg
 ¦--genome_annotation
 ¦   ¦--CNTX
 ¦   ¦   °--FEATURE_CNTX_genom_annotations.csv
 ¦   ¦--CNTX_genom_annotations.svg
 ¦   ¦--CX_annotation
 ¦   ¦   °--CNTX_CX_annotation.csv
 ¦   °--TEs_addiotionnal_results
 ¦       ¦--CNTX_TE.vs.ProteinCodingGenes.svg
 ¦       ¦--CNTX_TE_Super_Family.svg
 ¦       °--CNTX_TE_Super_Family_Freq.csv
 ¦--GO_analysis
 ¦   ¦--CNTX
 ¦   ¦   °--FEATURE
 ¦   ¦       ¦--GO_TYPE.gain.topGO.csv
 ¦   ¦       ¦--GO_TYPE.loss.topGO.csv
 ¦   ¦       ¦--GO_TYPE_gain_weight01.pdf
 ¦   ¦       °--GO_TYPE_loss_weight01.pdf
 ¦   °--CNTX_FEATURE_GO.svg
 ¦--KEGG_pathway
 ¦   ¦--CNTX
 ¦   ¦   °--FEATURE
 ¦   ¦       ¦--gain.KEGG.csv
 ¦   ¦       °--loss.KEGG.csv
 ¦   °--FEATURE.KEGG.svg
 ¦--metaPlots
 ¦   ¦--TEs
 ¦   ¦   ¦--metaPlot_tables
 ¦   ¦   ¦   ¦--mto1.CNTX.STREAM.csv
 ¦   ¦   ¦   °--wt.CNTX.STREAM.csv
 ¦   ¦   °--TEs_CNTX_metaPlot.svg
 ¦   ¦--Genes
 ¦   ¦   ¦--metaPlot_tables
 ¦   ¦   ¦   ¦--mto1.CNTX.STREAM.csv
 ¦   ¦   ¦   °--wt.CNTX.STREAM.csv
 ¦   ¦   °--Genes_CNTX_metaPlot.svg
 ¦   ¦--Gene_features
 ¦   ¦   ¦--metaPlot_tables
 ¦   ¦   ¦   ¦--mto1.CNTX.STREAM.csv
 ¦   ¦   ¦   °--wt.CNTX.STREAM.csv
 ¦   ¦   °--Gene_features_CNTX_metaPlot.svg
 °--methylation_levels
     ¦--Euchromatin_region.csv
     ¦--Euchromatin_region.svg
     ¦--Heterochromatin_region.csv
     ¦--Heterochromatin_region.svg
     ¦--Whole_Genome.csv                   
     °--Whole_Genome.svg
```
