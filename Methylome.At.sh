#!/bin/bash

samples_file=""

# Get the directory where the Bash script is located
Methylome_At_path=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")

# Default files path
annotation_file="$Methylome_At_path/annotation_files/Methylome.At_annotations.csv.gz"
description_file="$Methylome_At_path/annotation_files/Methylome.At_description_file.csv.gz"
TEs_file="$Methylome_At_path/annotation_files/TAIR10_Transposable_Elements.txt"

# Default values for optional arguments
minProportionDiff_CG=0.4
minProportionDiff_CHG=0.2
minProportionDiff_CHH=0.1
binSize=100
minCytosinesCount=4
minReadsPerCytosine=4
pValueThreshold=0.05
methyl_files_type=CX_report
img_type=svg
n_cores=10
GO_analysis=FALSE
KEGG_pathways=FALSE

# Function to display help text
usage() {
  echo ""
  echo "Usage: $0 [samples_file] [options]"
  echo ""
  echo "Required argument:"
  echo "  --samples_file                Path to samples file [required]"
  echo ""
  echo "Optional arguments:"
  echo "  --minProportionDiff_CG        Minimum proportion difference for CG [default: $minProportionDiff_CG]"
  echo "  --minProportionDiff_CHG       Minimum proportion difference for CHG [default: $minProportionDiff_CHG]"
  echo "  --minProportionDiff_CHH       Minimum proportion difference for CHH [default: $minProportionDiff_CHH]"
  echo "  --binSize                     DMRs bin size [default: $binSize]"
  echo "  --minCytosinesCount           Minimum cytosines count [default: $minCytosinesCount]"
  echo "  --minReadsPerCytosine         Minimum reads per cytosine [default: $minReadsPerCytosine]"
  echo "  --pValueThreshold             P-value threshold [default: $pValueThreshold]"
  echo "  --file_type                   Post-alignment file type - 'CX_report', 'bedMethyl' and 'CGmap' [default: '$methyl_files_type' OR determine automatically]"
  echo "  --image_type                  Output images file type [default: '$img_type']"
  echo "  --n_cores                     Number of cores [default: $n_cores]"
  echo "  --GO_analysis                 Perform GO analysis [default: $GO_analysis]"
  echo "  --KEGG_pathways               Perform KEGG pathways analysis [default: $KEGG_pathways]"
  echo "  --annotation_file             Genome Annotation file [default: Methylome.At annotations file (TAIR10 based)]"
  echo "  --description_file            Description file [default: Methylome.At description file]"
  echo "  --TEs_file                    Transposable Elements file [default: TAIR10 'Transposable Elements' annotations]"
  echo "  --Methylome_At_path           Path to Methylome.At [default: $Methylome_At_path]"
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
    --minProportionDiff_CG) minProportionDiff_CG="$2"; shift ;;
    --minProportionDiff_CHG) minProportionDiff_CHG="$2"; shift ;;
    --minProportionDiff_CHH) minProportionDiff_CHH="$2"; shift ;;
    --binSize) binSize="$2"; shift ;;
    --minCytosinesCount) minCytosinesCount="$2"; shift ;;
    --minReadsPerCytosine) minReadsPerCytosine="$2"; shift ;;
    --pValueThreshold) pValueThreshold="$2"; shift ;;
    --file_type) methyl_files_type="$2"; shift ;;
    --image_type) img_type="$2"; shift ;;
    --n_cores) n_cores="$2"; shift ;;
    --GO_analysis) GO_analysis="$2"; shift ;;
    --KEGG_pathways) KEGG_pathways="$2"; shift ;;
    --annotation_file) annotation_file="$2"; shift ;;
    --description_file) description_file="$2"; shift ;;
    --TEs_file) TEs_file="$2"; shift ;;
    --Methylome_At_path) Methylome_At_path="$2"; shift ;;
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

# ensure files unix line endings [sed -i 's/\r$//' "$samples_file"]
dos2unix "$samples_file" 2>/dev/null
dos2unix "$annotation_file" 2>/dev/null
dos2unix "$description_file" 2>/dev/null
dos2unix "$TEs_file" 2>/dev/null

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
  samples_path=$(cut -f2 "$samples_file" | head -n 1)
elif [ "$is_csv" -eq 2 ]; then
  samples_var=$(cut -d',' -f1 "$samples_file" | uniq)
  samples_path=$(cut -d',' -f2 "$samples_file" | head -n 1)
else
  echo "Error: Cannot read the samples file. Please provide a tab or comma separated file with two columns and no headers."
  exit 1
fi

# # Determine file type based on samples_path content
# if [[ "$methyl_files_type" != "CX_report" ]]; then 
#   if [[ $(basename "$samples_path" | awk -F. '{print $NF}') == "bed" ]]; then
#     methyl_files_type="bedMethyl"
#   elif [[ $(basename "$samples_path" .gz | awk -F. '{print $NF}') == "CGmap" ]]; then
#     methyl_files_type="CGmap"
#   else
#     echo "Error: Cannot read the samples file. Please provide a tab or comma separated file with two columns and no headers."
#     exit 1
#   fi
# fi

# Use head and sed to get the first and second lines
control_s=$(echo "$samples_var" | head -n1)
treatment_s=$(echo "$samples_var" | head -n2 | tail -n1)

# Output the configuration
echo ""
echo "$treatment_s VS. $control_s"
echo ""
echo "DMRs Min Proportion Diff CG: $minProportionDiff_CG"
echo "DMRs Min Proportion Diff CHG: $minProportionDiff_CHG"
echo "DMRs Min Proportion Diff CHH: $minProportionDiff_CHH"
echo "DMRs Bin size: $binSize"
echo "DMRs Min Cytosines Count: $minCytosinesCount"
echo "DMRs Min Reads Per Cytosine: $minReadsPerCytosine"
echo "DMRs P-value Threshold: $pValueThreshold"
echo ""
echo "Post-alignment file type: $methyl_files_type"
echo "Output images type: $img_type"
echo "Number of Cores: $n_cores"
echo "GO Analysis: $GO_analysis"
echo "KEGG Pathways: $KEGG_pathways"
echo ""
echo "Samples file: $samples_file"
echo "Annotation file: $annotation_file"
echo "Description file: $description_file"
echo "Transposable Elements file: $TEs_file"
echo "Methylome.At directory path: $Methylome_At_path"
echo ""
echo ""

# create results directory
mkdir -p results

# Generate log file with a timestamp
log_file="results/${treatment_s}_vs_${control_s}_$(date +"%d-%m-%y").log"
echo "**  $(date +"%d-%m-%y %H:%M")" > "$log_file"
echo "**  $treatment_s VS. $control_s" >> "$log_file"
echo "" >> "$log_file"

# Call the R script with the arguments
Rscript ./scripts/Methylome.At_run.R \
"$samples_file" \
"$Methylome_At_path" \
"$annotation_file" \
"$description_file" \
"$TEs_file" \
"$minProportionDiff_CG" \
"$minProportionDiff_CHG" \
"$minProportionDiff_CHH" \
"$binSize" \
"$minCytosinesCount" \
"$minReadsPerCytosine" \
"$pValueThreshold" \
"$methyl_files_type" \
"$img_type" \
"$n_cores" \
"$GO_analysis" \
"$KEGG_pathways" \
    2>> "$log_file"
    
# Output again the configurations, now to the 'log' file
echo "" >> "$log_file"
echo "" >> "$log_file"
echo "" >> "$log_file"
echo "configuration:" >> "$log_file"
echo "DMRs Min Proportion Diff CG: $minProportionDiff_CG" >> "$log_file"
echo "DMRs Min Proportion Diff CHG: $minProportionDiff_CHG" >> "$log_file"
echo "DMRs Min Proportion Diff CHH: $minProportionDiff_CHH" >> "$log_file"
echo "DMRs Bin size: $binSize" >> "$log_file"
echo "DMRs Min Cytosines Count: $minCytosinesCount" >> "$log_file"
echo "DMRs Min Reads Per Cytosine: $minReadsPerCytosine" >> "$log_file"
echo "DMRs P-value Threshold: $pValueThreshold" >> "$log_file"
echo "Post-alignment file type: $methyl_files_type" >> "$log_file"
echo "Output images type: $img_type" >> "$log_file"
echo "Number of Cores: $n_cores" >> "$log_file"
echo "GO Analysis: $GO_analysis" >> "$log_file"
echo "KEGG Pathways: $KEGG_pathways" >> "$log_file"
echo "Samples file: $samples_file" >> "$log_file"
echo "Annotation file: $annotation_file" >> "$log_file"
echo "Description file: $description_file" >> "$log_file"
echo "Transposable Elements file: $TEs_file" >> "$log_file"
echo "Methylome_At_path: $Methylome_At_path" >> "$log_file"