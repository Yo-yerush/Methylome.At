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

# System Requirements
* Using Conda environment
* * Linux Environment
* * [Conda](https://docs.conda.io/en/latest/miniconda.html) ([download](https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh))
* * [Whiptail](https://linux.die.net/man/1/whiptail) (for UI tutorial)
* * CPU: No special restrictions, but it is recommended to work with more than 10 cores for improved efficiency.
####
* using local R environment
* * [R](https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html) (≥ 4.4.0)
* * [RCurl](https://cran.r-project.org/web/packages/RCurl/index.html)
* * [textshaping](https://github.com/r-lib/textshaping) (≥ 0.4.1)
* * R packages (install 'install_R_packages.R' script)
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
#### 1. Download the source code.
```
git clone https://github.com/Yo-yerush/Methylome.At.git
```
```
cd ./Methylome.At
```
#### 2. Setup [conda](https://docs.conda.io/en/latest/miniconda.html) environment and install all the packages.
* using build-in setup script:
```
chmod +x ./setup_env.sh
./setup_env.sh
```
#####  *Check the 'setup_env' log file to verify packages installation*

### 
* or manually:
```
conda create --name Methylome.At_env
conda activate Methylome.At_env
conda install -c conda-forge -c bioconda r-base=4.4.2 r-curl r-rcurl r-devtools zlib r-textshaping harfbuzz fribidi freetype libpng pkg-config libxml2 r-xml bioconductor-rsamtools
Rscript scripts/install_R_packages.R
```

# Input files
#### 1. Samples table file ([example](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/sample_table_example.txt))

##### ****tab-delimited***, no header
```
treatment	PATH/TO/CX_report.txt
```

#### 2. '**CX_report**' file is an post-alignment methylation status for every cytosine in the genome, output from [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) program.
##### ****tab-delimited***, no header. See [columns definition](https://support.illumina.com/help/BaseSpace_App_MethylSeq_help/Content/Vault/Informatics/Sequencing_Analysis/Apps/swSEQ_mAPP_MethylSeq_CytosineReport.htm).
```
Chr1     3563    +       0       6       CHG     CCG
Chr1     3564    +       5       2       CG      CGA
Chr1     3565    -       2       3       CG      CGG
Chr1     3571    -       0       5       CHH     CAA
Chr1     3577    +       1       5       CHH     CTA
```

##### Convert 'CGmap' to 'CX_report' file:
```
./scripts/cgmap_to_cx.sh PATH/TO/input_file.CGmap PATH/TO/output_CX_report.txt
```

#### 3. Annotation and description files

* By default, Methylome.At provides *gene annotation file* in GFF3 format, *transposable elements annotation file* in text file (provide by [TAIR10](https://www.arabidopsis.org/)), and a *description file* integrated from multiple databases.
* users can alternativley use other *annotation files* (.gtf/.gff/.gff3) and *description file* (.csv/.txt)
##### description file columns:
```
gene_id Symbol	Short_description	Gene_description	Computational_description	AraCyc.Db	AraCyc.Name	gene_model_type	Protein.families	GO.biological.process	GO.cellular.component	GO.molecular.function	note	Derives_from	old_symbols	EC	KEGG_pathway	refseq_id	PMID

```

# Run Methylome.At with UI
Running this script will open a UI menu to run the **Methylome.At** main pipeline and the **MetaPlot** pipeline, allowing the user to change parameters as needed.
```
./Methylome.At_UI.sh
```

# Run Methylome.At manually
#### Main pipeline (**Methylome.At**):
```
./Methylome.At.sh PATH/TO/samples_table.txt
```
#### **MetaPlots** pipeline:
```
./Methylome.At_metaPlots.sh PATH/TO/samples_table.txt
```

#### Usage:
```
$ ./Methylome.At.sh --help

Usage: ./Methylome.At.sh [samples_file] [options]

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


```
$ ./Methylome.At_metaPlots.sh --help

Usage: ./Methylome.At_metaPlots.sh [samples_file] [options]

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
  --Methylome_At_path           Path to 'Methylome.At' directory [default: PATH/TO/Methylome.At]

```

## .log file
An automated log file will be created during the process. see examples from [Methylome.At](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Methylome.At_log_file.log) and [MetaPlots](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/MetaPlots_log_file.log) pipelines.


## Output Figures
*each analysis also output the data as a table

*all the results taken from [paper]()

#

#### **Total Methylation Levels**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Whole_Genome.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Heterochromatin_region.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Euchromatin_region.svg)

#

#### **Methylation Distribution**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/ChrPlot_CHG.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/ChrPlot_difference_CHG.svg)

#

#### **Gene Body and Transposable Elements Meta-plots**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/TEs_CHG_metaPlot.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Genes_CHG_metaPlot.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Genes_features_CHG_metaPlot.svg)

#

#### **DMRs direction (Gain or Loss)**

![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/CHH_gainORloss.svg)

#

#### **DMRs Distribution Mapping**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/DMRs_Density.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/ChrPlot_DMRs_CHG.svg)

#

#### **Genome Annotation for DMRs**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/CHG_genom_annotations.svg)

#

#### **Gene Ontology (GO) Analysis**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/CHG_Genes_GO.svg)

#

#### **KEGG Pathway Enrichment Test**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Promoters.KEGG.svg)

## Output Files tree

####  multiple files key:

> CNTX = CG, CHG, CHH

> FEATURE = gene, promoter, CDS, intron, 5'UTR, 3'UTR

> GO_TYPE = BP (biological process), MF (molecular function), CC (cellular component)

####  short output tree:
```
mto1_vs_wt
 ¦--conversion_rate.csv
 ¦--DMRs_CNTX_mto1_vs_wt.csv
 ¦--DMRs_Density_mto1_vs_wt.svg
 ¦--ChrPlot_CX/
 ¦--ChrPlot_DMRs/      
 ¦--gain_OR_loss/
 ¦--genome_annotation/
 ¦   ¦--CNTX/
 ¦   ¦--CX_annotation/
 ¦   ¦--CNTX_genom_annotations.svg
 ¦   °--TEs_addiotionnal_results/
 ¦--GO & KEGG analysis/
 ¦   ¦--CNTX/
 ¦   ¦   °--FEATURE/
 ¦   °--CNTX_FEATURE.svg
 ¦--metaPlots/
 ¦   ¦--TEs/Genes/Gene_features/
 ¦   ¦   ¦--metaPlot_tables/
 ¦   ¦   °--CNTX_metaPlot.svg
 °--methylation_levels/
```
* [long version of the output tree file](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/output_tree.txt)
