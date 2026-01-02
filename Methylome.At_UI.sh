#!/usr/bin/env bash

# Get the directory where the Bash script is located
Methylome_At_path=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
cd "$Methylome_At_path"

#####################################
# Initialize and Activate Conda Env
#####################################
if [ "$CONDA_DEFAULT_ENV" != "Methylome.At_env" ]; then
  eval "$(conda shell.bash hook)"
  conda activate Methylome.At_env
  
  if [ "$CONDA_DEFAULT_ENV" != "Methylome.At_env" ]; then
    echo "Error: Failed to activate the 'Methylome.At_env' Conda environment."
    echo "Please activate it manually using 'conda activate Methylome.At_env' and rerun the script."
    exit 1
  fi
fi

###############
# CONFIGURATION
###############
# Default parameters for run_bismark.sh:
SCRIPT_BIS_DEFAULT_genome="TAIR10"
SCRIPT_BIS_DEFAULT_ncores="8"

# Default parameters for Methylome.At.sh:
SCRIPT1_DEFAULT_minProportionDiff_CG="0.4"
SCRIPT1_DEFAULT_minProportionDiff_CHG="0.2"
SCRIPT1_DEFAULT_minProportionDiff_CHH="0.1"
SCRIPT1_DEFAULT_binSize="100"
SCRIPT1_DEFAULT_minCytosinesCount="4"
SCRIPT1_DEFAULT_minReadsPerCytosine="4"
SCRIPT1_DEFAULT_pValueThreshold="0.05"
SCRIPT1_DEFAULT_n_cores="8"
SCRIPT1_DEFAULT_GO_analysis="FALSE"
SCRIPT1_DEFAULT_KEGG_pathways="FALSE"
SCRIPT1_DEFAULT_file_type="CX_report"
SCRIPT1_DEFAULT_img_type="pdf"
SCRIPT1_DEFAULT_annotation_file="annotation_files/Methylome.At_annotations.csv.gz"
SCRIPT1_DEFAULT_description_file="annotation_files/Methylome.At_description_file.csv.gz"
SCRIPT1_DEFAULT_TEs_file="annotation_files/TAIR10_Transposable_Elements.txt"
SCRIPT1_DEFAULT_delta_H="FALSE"
SCRIPT1_DEFAULT_TEs_metaplots="FALSE"
SCRIPT1_DEFAULT_Genes_metaplots="FALSE"
SCRIPT1_DEFAULT_Gene_features_metaplots="FALSE"
SCRIPT1_DEFAULT_bin_size_features="10"
SCRIPT1_DEFAULT_metaPlot_random_genes="10000"


# Default parameters for Methylome.At_metaPlots.sh:
SCRIPT2_DEFAULT_Genes_n_TEs="TRUE"
SCRIPT2_DEFAULT_Gene_features="TRUE"
SCRIPT2_DEFAULT_minReadsPerCytosine="4"
SCRIPT2_DEFAULT_metaPlot_random_genes="10000"
SCRIPT2_DEFAULT_n_cores="8"
SCRIPT2_DEFAULT_bin_size_features="10"
SCRIPT2_DEFAULT_file_type="CX_report"
SCRIPT2_DEFAULT_img_type="pdf"
SCRIPT2_DEFAULT_annotation_file="annotation_files/Methylome.At_annotations.csv.gz"
SCRIPT2_DEFAULT_TEs_file="annotation_files/TAIR10_Transposable_Elements.txt"

# Paths to the scripts we want to run (adjust if needed)
SCRIPT_BIS_PATH="./scripts/run_bismark.sh"
SCRIPT1_PATH="./scripts/Methylome.At.sh"
SCRIPT2_PATH="./scripts/Methylome.At_metaPlots.sh"

##################
# WHIPTAIL DIALOGS
##################

# Prompt user: which scripts do you want to run?
CHOICE=$(whiptail --title "Choose scripts to run" \
  --checklist "Select which pipeline(s) to run. Use SPACE to toggle selection, ENTER to confirm, ESC to cancle." \
  18 70 4 \
  "Bismark" "Run genome alignment with Bismark" OFF \
  "Methylome.At" "Run main methylome pipeline'" ON \
  "MetaPlots" "Run 'metaPlots' without the main pipeline" OFF \
  3>&1 1>&2 2>&3)

# If user hits Cancel or ESC, exit
if [ $? -ne 0 ]; then
  echo "No scripts selected. Exiting."
  exit 1
fi

# Convert whiptail’s checklist output into an array
SELECTED_SCRIPTS=()
for item in $CHOICE; do
  # Remove surrounding quotes
  item=$(echo $item | sed 's/"//g')
  SELECTED_SCRIPTS+=("$item")
done

#############################
# Prompt for the samples_file
#############################
SAMPLES_FILE=$(whiptail --title "samples_file" --inputbox \
  "Enter the path to the samples table file:" \
  10 70 \
  "" \
  3>&1 1>&2 2>&3)
SAMPLES_FILE="${SAMPLES_FILE//\"/}"
if [ $? -ne 0 ] || [ -z "$SAMPLES_FILE" ]; then
  echo "samples_file is required. Exiting."
  exit 1
fi

########################################
# Function to edit parameters interactively
########################################

# A generic function to create a parameter menu for bismrk scripts
edit_script_bis_parameters() {
  # Parameters are expected to be set before calling this function
  while true; do
    OPTION=$(whiptail --title "'Bismark alignment' Parameters" --menu "Select a parameter to change or proceed with current settings." 25 78 16 \
      "Proceed." "Use current parameters" \
      "$SCRIPT_BIS_DEFAULT_genome" "Reference FASTA file" \
      "$SCRIPT_BIS_DEFAULT_ncores" "Number of cores" \
      3>&1 1>&2 2>&3)

    # Check if user cancelled
    [ $? -ne 0 ] && return 1

    if [ "$OPTION" = "Proceed." ]; then
      break
    elif [ "$OPTION" = "$SCRIPT_BIS_DEFAULT_genome" ]; then
      SCRIPT_BIS_DEFAULT_genome=$(whiptail --inputbox "Path to reference genome file" 10 70 "$SCRIPT_BIS_DEFAULT_genome" 3>&1 1>&2 2>&3 || echo "$SCRIPT_BIS_DEFAULT_genome")
    elif [ "$OPTION" = "$SCRIPT_BIS_DEFAULT_ncores" ]; then
      SCRIPT_BIS_DEFAULT_ncores=$(whiptail --inputbox "Number of cores" 10 70 "$SCRIPT_BIS_DEFAULT_ncores" 3>&1 1>&2 2>&3 || echo "$SCRIPT_BIS_DEFAULT_ncores")
    fi
  done
}

# A generic function to create a parameter menu for script1
edit_script1_parameters() {
    # nice aligned display helpers (monospace)
    fmt() { printf '%-10s  %-62s' "$1" "$2"; }

  # Parameters are expected to be set before calling this function
  while true; do
    OPTION=$(whiptail --title "'Methylome.At' Parameters" \
      --menu "Select a parameter to change or proceed with current settings." 30 90 22 \
      "Proceed."                "$(fmt '' 'Use current parameters')" \
      "Min diff (CG)"           "$(fmt "$SCRIPT1_minProportionDiff_CG" 'Min methylation proportion difference to call CG DMRs')" \
      "Min diff (CHG)"          "$(fmt "$SCRIPT1_minProportionDiff_CHG" 'Min methylation proportion difference to call CHG DMRs')" \
      "Min diff (CHH)"          "$(fmt "$SCRIPT1_minProportionDiff_CHH" 'Min methylation proportion difference to call CHH DMRs')" \
      "DMR bin size (bp)"       "$(fmt "$SCRIPT1_binSize" 'Bin-size window used for DMR calling')" \
      "Min cytosines / bin"     "$(fmt "$SCRIPT1_minCytosinesCount" 'Minimum cytosine count required in a bin to test it')" \
      "Min reads / cytosine"    "$(fmt "$SCRIPT1_minReadsPerCytosine" 'Minimum read coverage per cytosine to include it')" \
      "Adj. p-value cutoff"     "$(fmt "$SCRIPT1_pValueThreshold" 'Adjusted p-value threshold for calling significant DMRs')" \
      "CPU cores"               "$(fmt "$SCRIPT1_n_cores" 'Number of cores for parallel steps')" \
      "Input file format"       "$(fmt "$SCRIPT1_file_type" 'Methylation input format (CX_report/bedMethyl/CGmap)')" \
      "Figure format"           "$(fmt "$SCRIPT1_img_type" 'Output image format for plots (pdf/svg/png/tiff/...)')" \
      "GO enrichment"           "$(fmt "$SCRIPT1_GO_analysis" 'GO enrichment on DMR-associated gene-bodies/promoters')" \
      "KEGG enrichment"         "$(fmt "$SCRIPT1_KEGG_pathways" 'KEGG pathway enrichment on DMR-associated gene-bodies/promoters')" \
      "TE meta-plots"           "$(fmt "$SCRIPT1_TEs_metaplots" 'Metaplots over transposable elements (TE bodies/flanks)')" \
      "Gene-body meta-plots"    "$(fmt "$SCRIPT1_Genes_metaplots" 'Metaplots across gene bodies (TSS→TES/flank)')" \
      "Gene-feature meta-plots" "$(fmt "$SCRIPT1_Gene_features_metaplots" 'Metaplots over gene features (promoter/CDS/intron/UTR)')" \
      "Feature bin size"        "$(fmt "$SCRIPT1_bin_size_features" 'Bins-size for each gene feature region in feature meta-plots')" \
      "Random genes"            "$(fmt "$SCRIPT1_metaPlot_random_genes" 'Number of genes to sample for meta-plots (or all)')" \
      "dH analysis"                    "$(fmt "$SCRIPT1_delta_H" '')" \
      "Annotation file (gtf/gff/csv)"  "$(fmt "$SCRIPT1_annotation_file" '')" \
      "Gene descriptions (txt/csv)"    "$(fmt "$SCRIPT1_description_file" '')" \
      "TE annotation file (txt)"       "$(fmt "$SCRIPT1_TEs_file" '')" \
      3>&1 1>&2 2>&3)

    # Check if user cancelled
    # Check if user cancelled
    [ $? -ne 0 ] && return 1

    case "$OPTION" in
      "Proceed.")
        break
        ;;

      "Min diff (CG)")
        SCRIPT1_minProportionDiff_CG=$(
          whiptail --inputbox "Minimum methylation proportion difference for CG DMRs" 10 80 \
            "$SCRIPT1_minProportionDiff_CG" 3>&1 1>&2 2>&3
        ) || SCRIPT1_minProportionDiff_CG="$SCRIPT1_minProportionDiff_CG"
        ;;

      "Min diff (CHG)")
        SCRIPT1_minProportionDiff_CHG=$(
          whiptail --inputbox "Minimum methylation proportion difference for CHG DMRs" 10 80 \
            "$SCRIPT1_minProportionDiff_CHG" 3>&1 1>&2 2>&3
        ) || SCRIPT1_minProportionDiff_CHG="$SCRIPT1_minProportionDiff_CHG"
        ;;

      "Min diff (CHH)")
        SCRIPT1_minProportionDiff_CHH=$(
          whiptail --inputbox "Minimum methylation proportion difference for CHH DMRs" 10 80 \
            "$SCRIPT1_minProportionDiff_CHH" 3>&1 1>&2 2>&3
        ) || SCRIPT1_minProportionDiff_CHH="$SCRIPT1_minProportionDiff_CHH"
        ;;

      "DMR bin size (bp)")
        SCRIPT1_binSize=$(
          whiptail --inputbox "Bin-size (bp) for DMR calling" 10 80 \
            "$SCRIPT1_binSize" 3>&1 1>&2 2>&3
        ) || SCRIPT1_binSize="$SCRIPT1_binSize"
        ;;

      "Min cytosines / bin")
        SCRIPT1_minCytosinesCount=$(
          whiptail --inputbox "Minimum cytosines count per bin" 10 80 \
            "$SCRIPT1_minCytosinesCount" 3>&1 1>&2 2>&3
        ) || SCRIPT1_minCytosinesCount="$SCRIPT1_minCytosinesCount"
        ;;

      "Min reads / cytosine")
        SCRIPT1_minReadsPerCytosine=$(
          whiptail --inputbox "Minimum reads per cytosine" 10 80 \
            "$SCRIPT1_minReadsPerCytosine" 3>&1 1>&2 2>&3
        ) || SCRIPT1_minReadsPerCytosine="$SCRIPT1_minReadsPerCytosine"
        ;;

      "Adj. p-value cutoff")
        SCRIPT1_pValueThreshold=$(
          whiptail --inputbox "Adjusted p-value threshold for DMRs" 10 80 \
            "$SCRIPT1_pValueThreshold" 3>&1 1>&2 2>&3
        ) || SCRIPT1_pValueThreshold="$SCRIPT1_pValueThreshold"
        ;;

      "CPU cores")
        SCRIPT1_n_cores=$(
          whiptail --inputbox "Number of cores" 10 80 \
            "$SCRIPT1_n_cores" 3>&1 1>&2 2>&3
        ) || SCRIPT1_n_cores="$SCRIPT1_n_cores"
        ;;

      "Input file format")
        SCRIPT1_file_type=$(
          whiptail --radiolist "Select methylation file type:" 15 80 3 \
            "CX_report" "Bismark CX_report (.txt)"   $([ "$SCRIPT1_file_type" = "CX_report" ] && echo ON || echo OFF) \
            "bedMethyl" "BED methylation (.bed)"     $([ "$SCRIPT1_file_type" = "bedMethyl" ] && echo ON || echo OFF) \
            "CGmap"     "CGmap format (.CGmap)"      $([ "$SCRIPT1_file_type" = "CGmap" ] && echo ON || echo OFF) \
            3>&1 1>&2 2>&3
        ) || SCRIPT1_file_type="$SCRIPT1_file_type"
        ;;

      "Figure format")
        SCRIPT1_img_type=$(
          whiptail --radiolist "Select image format:" 18 80 6 \
            "pdf"  "Vector; publication-friendly"                     $([ "$SCRIPT1_img_type" = "pdf"  ] && echo ON || echo OFF) \
            "svg"  "Vector; editable"                                 $([ "$SCRIPT1_img_type" = "svg"  ] && echo ON || echo OFF) \
            "png"  "Raster (lossless)"                                $([ "$SCRIPT1_img_type" = "png"  ] && echo ON || echo OFF) \
            "tiff" "Raster (lossless - LZW); journal standard"        $([ "$SCRIPT1_img_type" = "tiff" ] && echo ON || echo OFF) \
            "jpeg" "Raster (lossy); small file"                       $([ "$SCRIPT1_img_type" = "jpeg" ] && echo ON || echo OFF) \
            "bmp"  "Raster (uncompressed); avoid"                     $([ "$SCRIPT1_img_type" = "bmp"  ] && echo ON || echo OFF) \
            3>&1 1>&2 2>&3
        ) || SCRIPT1_img_type="$SCRIPT1_img_type"
        ;;

      "GO enrichment")
        SCRIPT1_GO_analysis=$(
          whiptail --radiolist "Perform GO enrichment analysis?" 12 80 2 \
            "TRUE"  "Yes"  $([ "$SCRIPT1_GO_analysis" = "TRUE"  ] && echo ON || echo OFF) \
            "FALSE" "No"   $([ "$SCRIPT1_GO_analysis" = "FALSE" ] && echo ON || echo OFF) \
            3>&1 1>&2 2>&3
        ) || SCRIPT1_GO_analysis="$SCRIPT1_GO_analysis"
        ;;

      "KEGG enrichment")
        SCRIPT1_KEGG_pathways=$(
          whiptail --radiolist "Perform KEGG pathway enrichment analysis?" 12 80 2 \
            "TRUE"  "Yes"  $([ "$SCRIPT1_KEGG_pathways" = "TRUE"  ] && echo ON || echo OFF) \
            "FALSE" "No"   $([ "$SCRIPT1_KEGG_pathways" = "FALSE" ] && echo ON || echo OFF) \
            3>&1 1>&2 2>&3
        ) || SCRIPT1_KEGG_pathways="$SCRIPT1_KEGG_pathways"
        ;;

      "TE meta-plots")
        SCRIPT1_TEs_metaplots=$(
          whiptail --radiolist "Generate metaplots for Transposable Elements (TEs)?" 12 80 2 \
            "TRUE"  "Yes"  $([ "$SCRIPT1_TEs_metaplots" = "TRUE"  ] && echo ON || echo OFF) \
            "FALSE" "No"   $([ "$SCRIPT1_TEs_metaplots" = "FALSE" ] && echo ON || echo OFF) \
            3>&1 1>&2 2>&3
        ) || SCRIPT1_TEs_metaplots="$SCRIPT1_TEs_metaplots"
        ;;

      "Gene-body meta-plots")
        SCRIPT1_Genes_metaplots=$(
          whiptail --radiolist "Generate metaplots across gene bodies?" 12 80 2 \
            "TRUE"  "Yes"  $([ "$SCRIPT1_Genes_metaplots" = "TRUE"  ] && echo ON || echo OFF) \
            "FALSE" "No"   $([ "$SCRIPT1_Genes_metaplots" = "FALSE" ] && echo ON || echo OFF) \
            3>&1 1>&2 2>&3
        ) || SCRIPT1_Genes_metaplots="$SCRIPT1_Genes_metaplots"
        ;;

      "Gene-feature meta-plots")
        SCRIPT1_Gene_features_metaplots=$(
          whiptail --radiolist "Generate metaplots over gene features (promoter/CDS/UTRs/introns)?" 12 90 2 \
            "TRUE"  "Yes"  $([ "$SCRIPT1_Gene_features_metaplots" = "TRUE"  ] && echo ON || echo OFF) \
            "FALSE" "No"   $([ "$SCRIPT1_Gene_features_metaplots" = "FALSE" ] && echo ON || echo OFF) \
            3>&1 1>&2 2>&3
        ) || SCRIPT1_Gene_features_metaplots="$SCRIPT1_Gene_features_metaplots"
        ;;

      "Feature bin size")
        SCRIPT1_bin_size_features=$(
          whiptail --inputbox "Bin size for feature metaplots" 10 80 \
            "$SCRIPT1_bin_size_features" 3>&1 1>&2 2>&3
        ) || SCRIPT1_bin_size_features="$SCRIPT1_bin_size_features"
        ;;

      "Random genes")
        SCRIPT1_metaPlot_random_genes=$(
          whiptail --inputbox "Number of random genes for metaplots (or 'all')" 10 80 \
            "$SCRIPT1_metaPlot_random_genes" 3>&1 1>&2 2>&3
        ) || SCRIPT1_metaPlot_random_genes="$SCRIPT1_metaPlot_random_genes"
        ;;

      "dH analysis")
        SCRIPT1_delta_H=$(
          whiptail --radiolist "Enable dH analysis?" 12 80 2 \
            "TRUE"  "Yes"  $([ "$SCRIPT1_delta_H" = "TRUE"  ] && echo ON || echo OFF) \
            "FALSE" "No"   $([ "$SCRIPT1_delta_H" = "FALSE" ] && echo ON || echo OFF) \
            3>&1 1>&2 2>&3
        ) || SCRIPT1_delta_H="$SCRIPT1_delta_H"
        ;;

      "Annotation file (gtf/gff/csv)")
        SCRIPT1_annotation_file=$(
          whiptail --inputbox "Path to annotation file (GFF/GTF/CSV as required)" 10 90 \
            "$SCRIPT1_annotation_file" 3>&1 1>&2 2>&3
        ) || SCRIPT1_annotation_file="$SCRIPT1_annotation_file"
        ;;

      "Gene descriptions (txt/csv)")
        SCRIPT1_description_file=$(
          whiptail --inputbox "Path to gene description file" 10 90 \
            "$SCRIPT1_description_file" 3>&1 1>&2 2>&3
        ) || SCRIPT1_description_file="$SCRIPT1_description_file"
        ;;

      "TE annotation file (txt)")
        SCRIPT1_TEs_file=$(
          whiptail --inputbox "Path to TE annotation file" 10 90 \
            "$SCRIPT1_TEs_file" 3>&1 1>&2 2>&3
        ) || SCRIPT1_TEs_file="$SCRIPT1_TEs_file"
        ;;

      *)
        # Unknown / spacer items (if you add headers later)
        ;;
    esac
  done
}



###################
# Gather run_bismark.sh
###################
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Bismark" ]]; then
    # Initialize parameters with defaults
    SCRIPT_BIS_genome="$SCRIPT_BIS_DEFAULT_genome"
    SCRIPT_BIS_ncores="$SCRIPT_BIS_DEFAULT_ncores"

    # Directly go to the parameters selection menu
    edit_script_bis_parameters || exit 1
fi

###################
# Gather Methylome.At.sh
###################
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Methylome.At" ]]; then
    # Initialize parameters with defaults
    SCRIPT1_minProportionDiff_CG="$SCRIPT1_DEFAULT_minProportionDiff_CG"
    SCRIPT1_minProportionDiff_CHG="$SCRIPT1_DEFAULT_minProportionDiff_CHG"
    SCRIPT1_minProportionDiff_CHH="$SCRIPT1_DEFAULT_minProportionDiff_CHH"
    SCRIPT1_binSize="$SCRIPT1_DEFAULT_binSize"
    SCRIPT1_minCytosinesCount="$SCRIPT1_DEFAULT_minCytosinesCount"
    SCRIPT1_minReadsPerCytosine="$SCRIPT1_DEFAULT_minReadsPerCytosine"
    SCRIPT1_pValueThreshold="$SCRIPT1_DEFAULT_pValueThreshold"
    SCRIPT1_n_cores="$SCRIPT1_DEFAULT_n_cores"
    SCRIPT1_GO_analysis="$SCRIPT1_DEFAULT_GO_analysis"
    SCRIPT1_KEGG_pathways="$SCRIPT1_DEFAULT_KEGG_pathways"
    SCRIPT1_file_type="$SCRIPT1_DEFAULT_file_type"
    SCRIPT1_img_type="$SCRIPT1_DEFAULT_img_type"
    SCRIPT1_annotation_file="$SCRIPT1_DEFAULT_annotation_file"
    SCRIPT1_description_file="$SCRIPT1_DEFAULT_description_file"
    SCRIPT1_TEs_file="$SCRIPT1_DEFAULT_TEs_file"
    SCRIPT1_TEs_metaplots="$SCRIPT1_DEFAULT_TEs_metaplots"
    SCRIPT1_Genes_metaplots="$SCRIPT1_DEFAULT_Genes_metaplots"
    SCRIPT1_Gene_features_metaplots="$SCRIPT1_DEFAULT_Gene_features_metaplots"
    SCRIPT1_bin_size_features="$SCRIPT1_DEFAULT_bin_size_features"
    SCRIPT1_metaPlot_random_genes="$SCRIPT1_DEFAULT_metaPlot_random_genes"
    SCRIPT1_delta_H="$SCRIPT1_DEFAULT_delta_H"

    # Directly go to the parameters selection menu
    edit_script1_parameters || exit 1
fi

###################
# Gather Methylome.At_metaPlots.sh
###################
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "MetaPlots" ]]; then
    # Initialize parameters with defaults
    SCRIPT2_Genes_n_TEs="$SCRIPT2_DEFAULT_Genes_n_TEs"
    SCRIPT2_Gene_features="$SCRIPT2_DEFAULT_Gene_features"
    SCRIPT2_minReadsPerCytosine="$SCRIPT2_DEFAULT_minReadsPerCytosine"
    SCRIPT2_metaPlot_random_genes="$SCRIPT2_DEFAULT_metaPlot_random_genes"
    SCRIPT2_n_cores="$SCRIPT2_DEFAULT_n_cores"
    SCRIPT2_bin_size_features="$SCRIPT2_DEFAULT_bin_size_features"
    SCRIPT2_file_type="$SCRIPT2_DEFAULT_file_type"
    SCRIPT2_img_type="$SCRIPT2_DEFAULT_img_type"
    SCRIPT2_annotation_file="$SCRIPT2_DEFAULT_annotation_file"
    SCRIPT2_TEs_file="$SCRIPT2_DEFAULT_TEs_file"
    
    # Directly go to the parameters selection menu
    edit_script2_parameters || exit 1
fi

###########
# RUN SCRIPTS
###########

# Construct a message listing the chosen scripts
chosen_message=""
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Bismark" ]]; then
    chosen_message+="'Bismark' "
fi
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Methylome.At" ]]; then
    chosen_message+="'Methylome.At' "
fi
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "MetaPlots" ]]; then
    if [ -n "$chosen_message" ]; then
        chosen_message+="and "
    fi
    chosen_message+="'MetaPlots' "
fi

# Trim trailing space
chosen_message=$(echo "$chosen_message" | sed 's/[[:space:]]*$//')

# If somehow no scripts are chosen (shouldn't happen due to earlier checks), handle gracefully
if [ -z "$chosen_message" ]; then
    chosen_message="No scripts selected"
fi

# Display yes/no dialog
if (whiptail --title "All done!" --yesno "You have chosen to run: $chosen_message.\n\nWould you like to proceed?" 12 70); then

  cd "$Methylome_At_path"

  # Bismark pipeline for 'cx_report' files
  if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Bismark" ]]; then
    echo "Running run_bismark.sh..."
    SAMPLES_FILE_CX=$(bash "$SCRIPT_BIS_PATH" -s "$SAMPLES_FILE" -g "$SCRIPT_BIS_genome" -n "$SCRIPT_BIS_ncores" -o "${Methylome_At_path}/bismark_CX_reports" --cx --mat | tail -n1)
  else
    SAMPLES_FILE_CX="$SAMPLES_FILE"
  fi

  cd "$Methylome_At_path"

  # Methylome.At pipeline invocation
  if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Methylome.At" ]]; then
    echo "Running Methylome.At.sh..."
    bash "$SCRIPT1_PATH" \
      --samples_file "$SAMPLES_FILE_CX" \
      --minProportionDiff_CG "$SCRIPT1_minProportionDiff_CG" \
      --minProportionDiff_CHG "$SCRIPT1_minProportionDiff_CHG" \
      --minProportionDiff_CHH "$SCRIPT1_minProportionDiff_CHH" \
      --binSize "$SCRIPT1_binSize" \
      --minCytosinesCount "$SCRIPT1_minCytosinesCount" \
      --minReadsPerCytosine "$SCRIPT1_minReadsPerCytosine" \
      --pValueThreshold "$SCRIPT1_pValueThreshold" \
      --n_cores "$SCRIPT1_n_cores" \
      --GO_analysis "$SCRIPT1_GO_analysis" \
      --KEGG_pathways "$SCRIPT1_KEGG_pathways" \
      --file_type "$SCRIPT1_file_type" \
      --image_type "$SCRIPT1_img_type" \
      --annotation_file "$SCRIPT1_annotation_file" \
      --description_file "$SCRIPT1_description_file" \
      --TEs_file "$SCRIPT1_TEs_file" \
      --MP_TEs "$SCRIPT1_TEs_metaplots" \
      --MP_Genes "$SCRIPT1_Genes_metaplots" \
      --MP_Gene_features "$SCRIPT1_Gene_features_metaplots" \
      --MP_features_bin_size "$SCRIPT1_bin_size_features" \
      --metaPlot_random "$SCRIPT1_metaPlot_random_genes" \
      --dH "$SCRIPT1_delta_H"
  fi

  cd "$Methylome_At_path"

  # MetaPlots pipeline
  if [[ " ${SELECTED_SCRIPTS[*]} " =~ "MetaPlots" ]]; then
    echo "Running Methylome.At_metaPlots.sh..."
    bash "$SCRIPT2_PATH" \
      --samples_file "$SAMPLES_FILE_CX" \
      --Genes_n_TEs "$SCRIPT2_Genes_n_TEs" \
      --Gene_features "$SCRIPT2_Gene_features" \
      --minReadsPerCytosine "$SCRIPT2_minReadsPerCytosine" \
      --metaPlot_random_genes "$SCRIPT2_metaPlot_random_genes" \
      --n_cores "$SCRIPT2_n_cores" \
      --bin_size_features "$SCRIPT2_bin_size_features" \
      --file_type "$SCRIPT2_file_type" \
      --image_type "$SCRIPT2_img_type" \
      --annotation_file "$SCRIPT2_annotation_file" \
      --TEs_file "$SCRIPT2_TEs_file"
  fi
else

  echo "You chose not to run the scripts. Exiting."
  exit 0
fi