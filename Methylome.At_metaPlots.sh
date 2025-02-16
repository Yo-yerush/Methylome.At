#!/bin/bash

samples_file=""

# Get the directory where the Bash script is located
Methylome_At_path=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")

# Default files path
annotation_file="$Methylome_At_path/annotation_files/Methylome.At_annotations.csv.gz"
TEs_file="$Methylome_At_path/annotation_files/TAIR10_Transposable_Elements.txt"

# Default values for optional arguments
minReadsPerCytosine=6
metaPlot_random_genes=10000
n_cores=20
bin_size_features=10
Genes_n_TEs=TRUE
Gene_features=TRUE

# Function to display help text
usage() {
  echo ""
  echo "Usage: $0 [samples_file] [options]"
  echo ""
  echo "Required argument:"
  echo "  --samples_file                Path to samples file [required]"
  echo ""
  echo "Optional arguments:"
  echo "  --Genes_n_TEs                 Analyze Genes and TEs metaPlot [logical. default: $Genes_n_TEs]"
  echo "  --Gene_features               Analyze Gene Features metaPlot [logical. default: $Gene_features]"
  echo "  --minReadsPerCytosine         Minimum reads per cytosine [default: $minReadsPerCytosine]"
  echo "  --metaPlot_random_genes       Number of random genes for metaPlots. 'all' for all the coding-genes and TEs [default: $metaPlot_random_genes]"
  echo "  --n_cores                     Number of cores [default: $n_cores]"
  echo "  --bin_size_features           Bin-size (set only for 'Gene_features' analysis!) [default: $bin_size_features]"
  echo "  --annotation_file             Genome Annotation file [default: Methylome.At annotations file (TAIR10 based)]"
  echo "  --TEs_file                    Transposable Elements file [default: TAIR10 'Transposable Elements' annotations]"
  echo "  --Methylome_At_path           Path to 'Methylome.At' directory [default: $Methylome_At_path]"
  echo ""
  exit 1
}

# Assign first positional argument to samples_file if it does not start with --
if [[ -n "$1" && "$1" != --* ]]; then
    samples_file="$1"
    shift
fi

# Check if the first argument is 'help' or a help flag
#if [[ "$1" == "-h" || "$1" == "--help" ]]; then
#  usage
#  exit 1
#fi

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --samples_file) samples_file="$2"; shift ;;
    --minReadsPerCytosine) minReadsPerCytosine="$2"; shift ;;
    --metaPlot_random_genes) metaPlot_random_genes="$2"; shift ;;
    --n_cores) n_cores="$2"; shift ;;
    --bin_size_features) bin_size_features="$2"; shift ;;
    --annotation_file) annotation_file="$2"; shift ;;
    --TEs_file) TEs_file="$2"; shift ;;
    --Methylome_At_path) Methylome_At_path="$2"; shift ;;
    --Genes_n_TEs) Genes_n_TEs="$2"; shift ;;
    --Gene_features) Gene_features="$2"; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
  shift
done

# Check for required argument
if [ -z "$samples_file" ]; then
  echo "Error: --samples_file is a required argument."
  usage
fi

# Change to the Methylome_At_path directory
cd "$Methylome_At_path" || {
    echo "Error: Cannot change directory to $Methylome_At_path"
    exit 1
}

# Extract samples unique values from the first column
is_txt=$(awk -F'\t' 'NR==1 {print NF}' "$samples_file")
is_csv=$(awk -F',' 'NR==1 {print NF}' "$samples_file")
if [ "$is_txt" -eq 2 ]; then
    samples_var=$(cut -f1 "$samples_file" | uniq)
elif [ "$is_csv" -eq 2 ]; then
    samples_var=$(cut -d',' -f1 "$samples_file" | uniq)
else
    echo "Error: Cannot read the samples file. Please provide a tab or comma separated file with two columns and no headers."
    exit 1
fi

# Use head and sed to get the first and second lines
control_s=$(echo "$samples_var" | head -n1)
treatment_s=$(echo "$samples_var" | head -n2 | tail -n1)

# Output the configuration
echo ""
echo "$treatment_s VS. $control_s"
echo "Genes and TEs metaPlot: $Genes_n_TEs"
echo "Gene Features metaPlot: $Gene_features"
echo "Min Reads Per Cytosine: $minReadsPerCytosine"
echo "MetaPlot Random Genes: $metaPlot_random_genes"
echo "Number of Cores: $n_cores"
echo "Bin-Size for 'Gene_features' analysis: $bin_size_features"
echo "Samples file: $samples_file"
echo "Annotation file: Methylome.At (TAIR10 based)"
echo "Transposable Elements file: TAIR10"
echo "Methylome_At_path: $Methylome_At_path"
echo ""

# Generate log file with a timestamp
log_file="results/${treatment_s}_vs_${control_s}_metaPlots_$(date +"%d-%m-%y").log"
echo "**  $(date +"%d-%m-%y %H:%M")" >> "$log_file"
echo "**  $treatment_s VS. $control_s" >> "$log_file"
echo "" >> "$log_file"

# Call the R script with the arguments
Rscript ./scripts/MetaPlots_run.R \
"$samples_file" \
"$Methylome_At_path" \
"$annotation_file" \
"$TEs_file" \
"$minReadsPerCytosine" \
"$metaPlot_random_genes" \
"$n_cores" \
"$bin_size_features" \
"$Genes_n_TEs" \
"$Gene_features" \
    2>> "$log_file"
    
# Output again the configurations, now to the 'log' file
echo ""
echo ""
echo ""
echo "<$treatment_s VS. $control_s> configuration:" >> "$log_file"
echo "Genes and TEs metaPlot: $Genes_n_TEs" >> "$log_file"
echo "Gene Features metaPlot: $Gene_features" >> "$log_file"
echo "Min Reads Per Cytosine: $minReadsPerCytosine" >> "$log_file"
echo "MetaPlot Random Genes: $metaPlot_random_genes" >> "$log_file"
echo "Number of Cores: $n_cores" >> "$log_file"
echo "Bin-Size for 'Gene_features' analysis: $bin_size_features" >> "$log_file"
echo "Samples file: $samples_file" >> "$log_file"
echo "Annotation file: Methylome.At (TAIR10 based)" >> "$log_file"
echo "Transposable Elements file: TAIR10" >> "$log_file"
echo "Methylome_At_path: $Methylome_At_path" >> "$log_file"
