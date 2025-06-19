# Methylome.At

Methylome.At is a comprehensive, R-based pipeline for *Arabidopsis thaliana* that processes post-alignment WGBS data for CG, CHG and CHH DNA methylation contexts, identifies differentially methylated regions (DMRs, using [DMRcaller](https://github.com/nrzabet/DMRcaller) package) to replicates/single samples data, integrates multiple genomic resources for functional interpretation, and generates extensive visualizations and annotations to advance understanding of plant epigenetic regulation.

#

Methylome.At will produce few analysis and each analysis contains CG, CHG and CHH contexts:

* Chloroplast Conversion Rate
* PCA analysis
* Total Methylation Levels
* Methylation Distribution
* Gene Body and Transposable Elements Meta-plots
* DMRs Identification (using [DMRcaller](https://github.com/nrzabet/DMRcaller))
* DMRs Distribution Mapping
* Genome Annotation for DMRs and DMPs
* Gene Ontology (GO) Analysis
* KEGG Pathway Enrichment Test
* BigWig Files Generation for DMRs visualisation

<img
  src="https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/pipeline_scheme.png"
  alt="Methylome.At Flow"
  width="80%"
/>

#

# System Requirements

#### Using Conda environment

* Linux Environment
* [Conda](https://docs.conda.io/en/latest/miniconda.html) ([download](https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh))
* [Whiptail](https://linux.die.net/man/1/whiptail) (for UI tutorial)
* CPU: No special restrictions, but it is recommended to work with more than 10 cores for improved efficiency.

#### using local R environment

* [R](https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html) (≥ 4.4.0)
* [RCurl](https://cran.r-project.org/web/packages/RCurl/index.html)
* [textshaping](https://github.com/r-lib/textshaping) (≥ 0.4.1)
* R packages (install 'install_R_packages.R' script)

 ```text
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
RColorBrewer
circlize
 ```

#

# Installation

#### 1. Download the source code

```bash
git clone https://github.com/Yo-yerush/Methylome.At.git
cd ./Methylome.At
```

#### 2. Setup [conda](https://docs.conda.io/en/latest/miniconda.html) environment and install all the packages

* using build-in setup script:

```bash
chmod +x ./setup_env.sh
./setup_env.sh
```

###### *Check the 'setup_env' log file to verify packages installation*

###

* or manually:

```bash
conda create --name Methylome.At_env
conda activate Methylome.At_env
conda install -c conda-forge -c bioconda r-base=4.4.2 r-curl r-rcurl r-devtools zlib r-textshaping harfbuzz fribidi freetype libpng pkg-config libxml2 r-xml bioconductor-rsamtools
Rscript scripts/install_R_packages.R
chmod +x ./Methylome.At_UI.sh
chmod +x ./Methylome.At_metaPlots.sh
chmod +x ./Methylome.At.sh
chmod +x ./scripts/cgmap_to_cx.sh
chmod +x ./scripts/cx_to_cgmap.sh
```

#

# Input files

#### 1. Samples table file ([example](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/sample_table_example.txt))

###### ****tab-delimited***, no header

```text
treatment PATH/TO/CX_report.txt
```

#### 2. '**CX_report**' file is an post-alignment methylation status for every cytosine in the genome, output from [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) program

###### ****tab-delimited***, no header. See [columns definition](https://support.illumina.com/help/BaseSpace_App_MethylSeq_help/Content/Vault/Informatics/Sequencing_Analysis/Apps/swSEQ_mAPP_MethylSeq_CytosineReport.htm)

```text
Chr1     3563    +       0       6       CHG     CCG
Chr1     3564    +       5       2       CG      CGA
Chr1     3565    -       2       3       CG      CGG
Chr1     3571    -       0       5       CHH     CAA
Chr1     3577    +       1       5       CHH     CTA
```

##### Convert 'CGmap' to 'CX_report' file

```bash
./scripts/cgmap_to_cx.sh PATH/TO/input_file.CGmap PATH/TO/output_CX_report.txt
```

#### 3. Annotation and description files

* By default, Methylome.At provides *gene annotation file* in GFF3 format, *transposable elements annotation file* in text file (provide by [TAIR10](https://www.arabidopsis.org/)), and a *description file* integrated from multiple databases.
* users can alternativley use other *annotation files* (.gtf/.gff/.gff3/.csv/.txt) and *description file* (.csv/.txt)

##### Description file columns

```text
gene_id Symbol Short_description Gene_description Computational_description AraCyc.Db AraCyc.Name gene_model_type Protein.families GO.biological.process GO.cellular.component GO.molecular.function note Derives_from old_symbols EC KEGG_pathway refseq_id PMID

```

#

# Run Methylome.At with UI

Running this script will open a UI menu to run the **Methylome.At** main pipeline and the **MetaPlot** pipeline, allowing the user to change parameters as needed.

```bash
./Methylome.At_UI.sh
```

#

# Run Methylome.At manually

#### Main pipeline (**Methylome.At**)

```bash
./Methylome.At.sh PATH/TO/samples_table.txt
```

#### **MetaPlots** pipeline

```bash
./Methylome.At_metaPlots.sh PATH/TO/samples_table.txt
```

#### Usage

```text
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
  --minReadsPerCytosine         Minimum reads per cytosine [default: 4]
  --pValueThreshold             P-value threshold [default: 0.05]
  --n_cores                     Number of cores [default: 10]
  --GO_analysis                 Perform GO analysis [default: FALSE]
  --KEGG_pathways               Perform KEGG pathways analysis [default: FALSE]
  --annotation_file             Genome Annotation file [default: Methylome.At annotations file (TAIR10 based)]
  --description_file            Description file [default: Methylome.At description file]
  --TEs_file                    Transposable Elements file [default: TAIR10 'Transposable Elements' annotations]
  --Methylome_At_path           Path to Methylome.At [default: PATH/TO/Methylome.At]
```

```text
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

#

## .log file

An automated log file will be created during the process. see examples from [Methylome.At](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Methylome.At_log_file.log) and [MetaPlots](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/MetaPlots_log_file.log) pipelines.

#

## Output Tabels and Figures

*all the results taken from [paper]()

#

#### **Chloroplast Conversion Rate**

```text
| sample | conversion_rate |
|--------|-----------------|
| wt     | 99.49           |
| wt     | 99.39           |
| mto1   | 99.48           |
| mto1   | 99.50           |
| mto1   | 99.47           |

```

#

#### **PCA plots**
>
> PCA plots and data tables for CG, CHG, CHH and all contexts

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/PCA_all_contexts.svg)

#

#### **Total Methylation Levels**
>
> csv tabel for whole-genome, eu-/hetero- chromatin

```text
| type | treatment | levels       | SD          |
|-----:|:---------:|-------------:|------------:|
| CG   | wt        | 25.82221842  | 0.067418241 |
| CHG  | wt        | 6.837569844  | 0.335858464 |
| CHH  | wt        | 2.264188479  | 0.160429206 |
| CG   | mto1      | 25.94971156  | 0.203546496 |
| CHG  | mto1      | 8.239533867  | 0.617113840 |
| CHH  | mto1      | 2.540204156  | 0.119720471 |
```

<img
  src="https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Whole_Genome.svg"
  alt="fig"
  height="40%"
/>
<img
  src="https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Heterochromatin_region.svg"
  alt="fig"
  height="40%"
/>
<img
  src="https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Euchromatin_region.svg"
  alt="fig"
  height="40%"
/>

#

#### **Methylation Distribution**

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/ChrPlot_CHG.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/ChrPlot_difference_CHG.svg)

#

#### **Gene Body and Transposable Elements Meta-plots**
>
> csv tables for up.stream, gene/TE-body, down.stream region; each context; each of genes/gene-feature/TEs
> first 4 rows. position based on bin-size (20 as default).

```text
| seqnames   | start | end | width | strand | Proportion  |
|:----------:|------:|----:|------:|:------:|------------:|
| up.stream  |     1 |   1 |     1 |   *    | 0.038614587 |
| up.stream  |     2 |   2 |     1 |   *    | 0.039872008 |
| up.stream  |     3 |   3 |     1 |   *    | 0.040454043 |
| up.stream  |     4 |   4 |     1 |   *    | 0.041221039 |
```

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/TEs_CHG_metaPlot.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Genes_CHG_metaPlot.svg)
 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Genes_features_CHG_metaPlot.svg)

#

#### **DMRs direction (Gain or Loss)**

![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/CHH_gainORloss.svg)

#

#### **DMRs Distribution Mapping**
>
> csv tables for each context
> first 5 rows

```text
| seqnames |  start  |   end   | width | strand |    pValue    |    log2FC    | context | sumReadsM1 | sumReadsN1 | proportion1 | sumReadsM2 | sumReadsN2 | proportion2 | cytosinesCount | direction | regionType |
|:--------:|--------:|--------:|------:|:------:|-------------:|-------------:|:-------:|-----------:|-----------:|------------:|-----------:|-----------:|------------:|---------------:|----------:|:----------:|
| Chr1     |  56101  |  56200  |   100 | *      | 7.32E-05     | 0.966553517  | CHG     | 21         | 75         | 0.285714286 | 66         | 118        | 0.558333333 | 4              | 1         | gain       |
| Chr1     | 213701  | 213800  |   100 | *      | 0.000208476  | 3.608941364  | CHG     | 1          | 110        | 0.017857143 | 38         | 177        | 0.217877095 | 8              | 1         | gain       |
| Chr1     | 239101  | 239200  |   100 | *      | 4.06E-05     | 0.496927653  | CHG     | 42         | 73         | 0.573333333 | 88         | 108        | 0.809090909 | 4              | 1         | gain       |
| Chr1     | 432801  | 432900  |   100 | *      | 0.001693563  | 1.007546809  | CHG     | 10         | 44         | 0.239130435 | 49         | 102        | 0.480769231 | 4              | 1         | gain       |
| Chr1     | 612601  | 612700  |   100 | *      | 0.040449983  | 0.827360614  | CHG     | 18         | 66         | 0.279411765 | 58         | 117        | 0.495798319 | 5              | 1         | gain       |
```

<img
  src="https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/DMRs_Density.svg"
  alt="fig"
  width="60%"
/>

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/ChrPlot_DMRs_CHG.svg)

#

#### **Genome Annotation for DMRs**
>
> csv tables for each context; each gene-feature
> first 5 rows

```text
| gene_id   | seqnames |   start  |   end    | width | strand |  pValue   |   log2FC   | context | direction | regionType |  type  | Symbol | Short_description                                   | Gene_description | Computational_description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | AraCyc.Db | AraCyc.Name | gene_model_type | Protein.families | GO.biological.process | GO.cellular.component | GO.molecular.function | note | Derives_from | old_symbols | EC | KEGG_pathway |           refseq_id           |    PMID    | sumReadsM1 | sumReadsN1 | proportion1  | sumReadsM2 | sumReadsN2 | proportion2   | cytosinesCount |
|-----------|---------:|---------:|---------:|------:|:------:|----------:|------------:|:-------:|----------:|:----------:|:------:|:------:|:----------------------------------------------------|:----------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:---------:|:----------:|:---------------:|:-----------------:|:---------------------:|:---------------------:|:---------------------:|:----:|:-----------:|:----------:|:--:|:------------:|:---------------------------------:|:----------:|-----------:|-----------:|-------------:|-----------:|-----------:|--------------:|----------------:|
| AT1G20400 | Chr1     |  7073626 |  7073725 |   100 |   *    | 9.32E-06  | 0.6119295  | CHG     |         1 | gain       | intron |        | Protein of unknown function (DUF1204)               |                  | FUNCTIONS IN: molecular_function unknown; INVOLVED IN: biological_process unknown; LOCATED IN: cellular_component unknown; EXPRESSED IN: cultured cell; CONTAINS InterPro DOMAIN/s: Protein of unknown function DUF1204 (InterPro:IPR009596); BEST Arabidopsis thaliana protein match is: myosin heavy chain-related (TAIR:AT2G15420.1); Has 2374 Blast hits to 2027 proteins in 308 species: Archae - 22; Bacteria - 170; Metazoa - 962; Fungi - 207; Plants - 198; Viruses - 74; Other Eukaryotes - 741 (source: NCBI BLink).                                                                                                      |           |             | protein_coding  |                  |                      |                      |                      |      |             |            |    |             | NM_101891; NP_173465            | 28441590   |         41 |        106 | 0.38888889   |         62 |        104 | 0.594339623    |               4 |
| AT1G29620 | Chr1     | 10348093 | 10348192 |   100 |   *    | 4.28E-03  | 0.5191947  | CHG     |         1 | gain       | intron |        | Cytochrome C oxidase polypeptide VIB family protein |                  | Cytochrome C oxidase polypeptide VIB family protein; BEST Arabidopsis thaliana protein match is: Cytochrome C oxidase polypeptide VIB family protein (TAIR:AT1G32720.1); Has 17 Blast hits to 17 proteins in 3 species: Archae - 0; Bacteria - 0; Metazoa - 0; Fungi - 0; Plants - 17; Viruses - 0; Other Eukaryotes - 0 (source: NCBI BLink).                                                                                                                                                                                                                                                                                                                                                         |           |             | protein_coding  |                  |                      |                      |                      |      |             |            |    |             | NM_102702; NP_174255            |            |         32 |         65 | 0.49253731   |         59 |         83 | 0.705882353    |               6 |
| AT1G30060 | Chr1     | 10545193 | 10545292 |   100 |   *    | 2.96E-07  | 1.2755803  | CHG     |         1 | gain       | intron |        | COP1-interacting protein-related                    |                  | COP1-interacting protein-related; BEST Arabidopsis thaliana protein match is: COP1-interacting protein-related (TAIR:AT2G01800.1); Has 95 Blast hits to 93 proteins in 18 species: Archae - 0; Bacteria - 0; Metazoa - 6; Fungi - 0; Plants - 89; Viruses - 0; Other Eukaryotes - 0 (source: NCBI BLink).<br>acyl carrier protein metabolism<br>holo-[acyl-carrier-protein] synthase                                                                                                                                                                                                                                                                                          |           |             | protein_coding  |                  |                      |                      |                      |      |             |            |    |             | NM_102745; NP_174298            |            |         31 |        221 | 0.14349776   |        106 |        306 | 0.347402597    |               7 |
| AT1G30060 | Chr1     | 10545593 | 10545692 |   100 |   *    | 3.99E-07  | 0.9102812  | CHG     |         1 | gain       | intron |        | COP1-interacting protein-related                    |                  | COP1-interacting protein-related; BEST Arabidopsis thaliana protein match is: COP1-interacting protein-related (TAIR:AT2G01800.1); Has 95 Blast hits to 93 proteins in 18 species: Archae - 0; Bacteria - 0; Metazoa - 6; Fungi - 0; Plants - 89; Viruses - 0; Other Eukaryotes - 0 (source: NCBI BLink).<br>acyl carrier protein metabolism<br>holo-[acyl-carrier-protein] synthase                                                                                                                                                                                                                                                                                          |           |             | protein_coding  |                  |                      |                      |                      |      |             |            |    |             | NM_102745; NP_174298            |            |         24 |         88 | 0.27777778   |         70 |        134 | 0.522058824    |               4 |
```

###

> TEs superfamily frequency

```text
| Transposon_Super_Family | unique_IDs | total_DMRs | hyper_DMRs | hypo_DMRs |
|:------------------------|-----------:|-----------:|-----------:|----------:|
| LTR/Gypsy               |        927 |       1706 |       1705 |         1 |
| DNA/MuDR                |        308 |        581 |        573 |         8 |
| LTR/Copia               |        233 |        536 |        534 |         2 |
| LINE/L1                 |        172 |        344 |        342 |         2 |
```

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/CHG_genom_annotations.svg)

#

#### **Gene Ontology (GO) Analysis**
>
> csv tables for each context; each gene-feature; BP/MF/CC; gain/loss of methylation (compare to control)

```text
| GO.ID      | Term                                 | Annotated | Significant | Expected | Fisher   |
|:----------:|:-------------------------------------|----------:|------- ----:|---------:|---------:|
| GO:0036071 | N-glycan fucosylation                |         2 |           1 |        0 | 0.00047  |
| GO:0001887 | selenium compound metabolic process  |         4 |           1 |        0 | 0.00093  |
| GO:0019346 | transsulfuration                     |         4 |           1 |        0 | 0.00093  |
| GO:0009086 | methionine biosynthetic process      |        22 |           1 |     0.01 | 0.00511  |
```

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/CHG_Genes_GO.svg)

#

#### **KEGG Pathway Enrichment Test**
>
> csv tables for each context; each gene-feature; gain/loss of methylation (compare to control)

```text
| pathway.code | pathway.name                      |   p.value   | Significant | Annotated |
|:------------:|:----------------------------------|------------:|------------:|----------:|
| ath00261     | Monobactam biosynthesis           | 5.46E-06    |           3 |        14 |
| ath00232     | Caffeine metabolism               | 0.00030362  |           1 |         3 |
| ath03008     | Ribosome biogenesis in eukaryotes | 0.001228952 |           7 |        92 |
| ath04144     | Endocytosis                       | 0.001590624 |          10 |       158 |
| ath00565     | Ether lipid metabolism            | 0.001741456 |           3 |        26 |
```

 ![fig](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/Promoters.KEGG.svg)

##

## Output Files tree

#### multiple files key

> CNTX = CG, CHG, CHH

> FEATURE = gene, promoter, CDS, intron, 5'UTR, 3'UTR

> GO_TYPE = BP (biological process), MF (molecular function), CC (cellular component)

#### short output tree

```text
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
 ¦   °--TEs/Genes/Gene_features/
 ¦       ¦--metaPlot_tables/
 ¦       °--CNTX_metaPlot.svg
 ¦--PCA_plots/
 °--methylation_levels/
```

* [long version of the output tree file](https://github.com/Yo-yerush/Methylome.At/blob/main/output_example/output_tree.txt)

###

###

###

# License

This project is licensed under the [MIT License](https://github.com/Yo-yerush/Methylome.At/blob/main/LICENSE).
