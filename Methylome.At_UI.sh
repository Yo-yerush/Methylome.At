#!/usr/bin/env bash
#
# A Bash script that uses whiptail to configure and run 'Methylome.At.sh' and 'Methylome.At_metaPlots.sh'
# Make sure you have `whiptail` installed (commonly in the `newt` package).
# Adjust the paths to 'Methylome.At.sh' and 'Methylome.At_metaPlots.sh' if needed.

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
# Default parameters for Methylome.At.sh:
SCRIPT1_DEFAULT_minProportionDiff_CG="0.4"
SCRIPT1_DEFAULT_minProportionDiff_CHG="0.2"
SCRIPT1_DEFAULT_minProportionDiff_CHH="0.1"
SCRIPT1_DEFAULT_binSize="100"
SCRIPT1_DEFAULT_minCytosinesCount="4"
SCRIPT1_DEFAULT_minReadsPerCytosine="6"
SCRIPT1_DEFAULT_pValueThreshold="0.05"
SCRIPT1_DEFAULT_n_cores="10"
SCRIPT1_DEFAULT_GO_analysis="TRUE"
SCRIPT1_DEFAULT_KEGG_pathways="TRUE"
SCRIPT1_DEFAULT_annotation_file="annotation_files/Methylome.At_annotations.csv.gz"
SCRIPT1_DEFAULT_description_file="annotation_files/Methylome.At_description_file.csv.gz"
SCRIPT1_DEFAULT_TEs_file="annotation_files/TAIR10_Transposable_Elements.txt"

# Default parameters for Methylome.At_metaPlots.sh:
SCRIPT2_DEFAULT_Genes_n_TEs="TRUE"
SCRIPT2_DEFAULT_Gene_features="TRUE"
SCRIPT2_DEFAULT_minReadsPerCytosine="6"
SCRIPT2_DEFAULT_metaPlot_random_genes="10000"
SCRIPT2_DEFAULT_n_cores="20"
SCRIPT2_DEFAULT_bin_size_features="10"
SCRIPT2_DEFAULT_annotation_file="annotation_files/Methylome.At_annotations.csv.gz"
SCRIPT2_DEFAULT_TEs_file="annotation_files/TAIR10_Transposable_Elements.txt"

# Paths to the scripts we want to run (adjust if needed)
SCRIPT1_PATH="./Methylome.At.sh"
SCRIPT2_PATH="./Methylome.At_metaPlots.sh"

##################
# WHIPTAIL DIALOGS
##################

# Prompt user: which scripts do you want to run?
CHOICE=$(whiptail --title "Choose scripts to run" \
  --checklist "Select which pipeline(s) to run. Use SPACE to toggle selection, ENTER to confirm, ESC to cancle." \
  18 70 4 \
  "Methylome.At" "Run main methylome pipeline'" ON \
  "metaPlots" "Run 'metaPlot' pipeline" OFF \
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
  "Enter the path to the samples file (required):" \
  10 70 \
  "" \
  3>&1 1>&2 2>&3)
if [ $? -ne 0 ] || [ -z "$SAMPLES_FILE" ]; then
  echo "samples_file is required. Exiting."
  exit 1
fi

########################################
# Function to edit parameters interactively
########################################

# A generic function to create a parameter menu for script1
edit_script1_parameters() {
  # Parameters stored in variables
  # They are expected to be set before calling this function
  while true; do
    OPTION=$(whiptail --title "'Methylome.At' Parameters" --menu "Select a parameter to change or proceed with current settings." 25 78 15 \
      "Use these parameters" "Proceed with current settings." \
      "minProportionDiff_CG: $SCRIPT1_minProportionDiff_CG" "Minimum proportion difference for CG" \
      "minProportionDiff_CHG: $SCRIPT1_minProportionDiff_CHG" "Minimum proportion difference for CHG" \
      "minProportionDiff_CHH: $SCRIPT1_minProportionDiff_CHH" "Minimum proportion difference for CHH" \
      "binSize: $SCRIPT1_binSize" "DMRs bin size" \
      "minCytosinesCount: $SCRIPT1_minCytosinesCount" "Minimum cytosines count" \
      "minReadsPerCytosine: $SCRIPT1_minReadsPerCytosine" "Minimum reads per cytosine" \
      "pValueThreshold: $SCRIPT1_pValueThreshold" "P-value threshold" \
      "n_cores: $SCRIPT1_n_cores" "Number of cores" \
      "GO_analysis: $SCRIPT1_GO_analysis" "Perform GO analysis (TRUE/FALSE)" \
      "KEGG_pathways: $SCRIPT1_KEGG_pathways" "Perform KEGG pathways analysis (TRUE/FALSE)" \
      "annotation_file: $SCRIPT1_annotation_file" "Path to annotation file" \
      "description_file: $SCRIPT1_description_file" "Path to description file" \
      "TEs_file: $SCRIPT1_TEs_file" "Path to Transposable Elements file" \
      3>&1 1>&2 2>&3)

    # Check if user cancelled
    [ $? -ne 0 ] && return 1

case "$OPTION" in
  "Use these parameters")
    # User is done editing
    break
    ;;
  minProportionDiff_CG*) 
    SCRIPT1_minProportionDiff_CG=$(whiptail --inputbox "Minimum proportion difference for CG" 10 70 "$SCRIPT1_minProportionDiff_CG" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_minProportionDiff_CG")
    ;;
  minProportionDiff_CHG*) 
    SCRIPT1_minProportionDiff_CHG=$(whiptail --inputbox "Minimum proportion difference for CHG" 10 70 "$SCRIPT1_minProportionDiff_CHG" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_minProportionDiff_CHG")
    ;;
  minProportionDiff_CHH*)
    SCRIPT1_minProportionDiff_CHH=$(whiptail --inputbox "Minimum proportion difference for CHH" 10 70 "$SCRIPT1_minProportionDiff_CHH" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_minProportionDiff_CHH")
    ;;
  binSize*)
    SCRIPT1_binSize=$(whiptail --inputbox "DMRs bin size" 10 70 "$SCRIPT1_binSize" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_binSize")
    ;;
  minCytosinesCount*)
    SCRIPT1_minCytosinesCount=$(whiptail --inputbox "Minimum cytosines count" 10 70 "$SCRIPT1_minCytosinesCount" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_minCytosinesCount")
    ;;
  minReadsPerCytosine*)
    SCRIPT1_minReadsPerCytosine=$(whiptail --inputbox "Minimum reads per cytosine" 10 70 "$SCRIPT1_minReadsPerCytosine" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_minReadsPerCytosine")
    ;;
  pValueThreshold*)
    SCRIPT1_pValueThreshold=$(whiptail --inputbox "P-value threshold" 10 70 "$SCRIPT1_pValueThreshold" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_pValueThreshold")
    ;;
  n_cores*)
    SCRIPT1_n_cores=$(whiptail --inputbox "Number of cores" 10 70 "$SCRIPT1_n_cores" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_n_cores")
    ;;
  GO_analysis*)
    SCRIPT1_GO_analysis=$(whiptail --radiolist "Perform GO analysis?" 12 70 2 \
      "TRUE" "Perform GO analysis" ON \
      "FALSE" "Skip GO analysis" OFF \
      3>&1 1>&2 2>&3 || echo "$SCRIPT1_GO_analysis")
    ;;
  KEGG_pathways*)
    SCRIPT1_KEGG_pathways=$(whiptail --radiolist "Perform KEGG pathways analysis?" 12 70 2 \
      "TRUE" "Perform KEGG pathways analysis" ON \
      "FALSE" "Skip KEGG pathways analysis" OFF \
      3>&1 1>&2 2>&3 || echo "$SCRIPT1_KEGG_pathways")
    ;;
  annotation_file*)
    SCRIPT1_annotation_file=$(whiptail --inputbox "Path to annotation file" 10 70 "$SCRIPT1_annotation_file" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_annotation_file")
    ;;
  description_file*)
    SCRIPT1_description_file=$(whiptail --inputbox "Path to description file" 10 70 "$SCRIPT1_description_file" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_description_file")
    ;;
  TEs_file*)
    SCRIPT1_TEs_file=$(whiptail --inputbox "Path to Transposable Elements file" 10 70 "$SCRIPT1_TEs_file" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_TEs_file")
    ;;
esac
  done
}

# A generic function to create a parameter menu for script2
edit_script2_parameters() {
  while true; do
    OPTION=$(whiptail --title "'MetaPlots' Parameters" --menu "Select a parameter to change or proceed with current settings." 25 78 12 \
      "Use these parameters" "Proceed with current settings." \
      "Genes_n_TEs ($SCRIPT2_Genes_n_TEs)" "Analyze Genes and TEs metaPlot (TRUE/FALSE)" \
      "Gene_features ($SCRIPT2_Gene_features)" "Analyze Gene Features metaPlot (TRUE/FALSE)" \
      "minReadsPerCytosine ($SCRIPT2_minReadsPerCytosine)" "Minimum reads per cytosine" \
      "metaPlot_random_genes ($SCRIPT2_metaPlot_random_genes)" "Number of random genes ('all' = all)" \
      "n_cores ($SCRIPT2_n_cores)" "Number of cores" \
      "bin_size_features ($SCRIPT2_bin_size_features)" "Bin-size for Gene_features analysis" \
      "annotation_file ($SCRIPT2_annotation_file)" "Path to genome annotation file" \
      "TEs_file ($SCRIPT2_TEs_file)" "Path to Transposable Elements file" \
      3>&1 1>&2 2>&3)

    [ $? -ne 0 ] && return 1

    case "$OPTION" in
  "Use these parameters")
    break
    ;;
  Genes_n_TEs*)
    SCRIPT2_Genes_n_TEs=$(whiptail --radiolist "Analyze Genes and TEs metaPlot?" 12 70 2 \
      "TRUE" "Yes" ON \
      "FALSE" "No" OFF \
      3>&1 1>&2 2>&3 || echo "$SCRIPT2_Genes_n_TEs")
    ;;
  Gene_features*)
    SCRIPT2_Gene_features=$(whiptail --radiolist "Analyze Gene Features metaPlot?" 12 70 2 \
      "TRUE" "Yes" ON \
      "FALSE" "No" OFF \
      3>&1 1>&2 2>&3 || echo "$SCRIPT2_Gene_features")
    ;;
  minReadsPerCytosine*)
    SCRIPT2_minReadsPerCytosine=$(whiptail --inputbox "Minimum reads per cytosine" 10 70 "$SCRIPT2_minReadsPerCytosine" 3>&1 1>&2 2>&3 || echo "$SCRIPT2_minReadsPerCytosine")
    ;;
  metaPlot_random_genes*)
    SCRIPT2_metaPlot_random_genes=$(whiptail --inputbox "Number of random genes ('all' for all)" 10 70 "$SCRIPT2_metaPlot_random_genes" 3>&1 1>&2 2>&3 || echo "$SCRIPT2_metaPlot_random_genes")
    ;;
  n_cores*)
    SCRIPT2_n_cores=$(whiptail --inputbox "Number of cores" 10 70 "$SCRIPT2_n_cores" 3>&1 1>&2 2>&3 || echo "$SCRIPT2_n_cores")
    ;;
  bin_size_features*)
    SCRIPT2_bin_size_features=$(whiptail --inputbox "Bin-size for Gene_features analysis" 10 70 "$SCRIPT2_bin_size_features" 3>&1 1>&2 2>&3 || echo "$SCRIPT2_bin_size_features")
    ;;
  annotation_file*)
    SCRIPT2_annotation_file=$(whiptail --inputbox "Path to genome annotation file" 10 70 "$SCRIPT2_annotation_file" 3>&1 1>&2 2>&3 || echo "$SCRIPT2_annotation_file")
    ;;
  TEs_file*)
    SCRIPT2_TEs_file=$(whiptail --inputbox "Path to Transposable Elements file" 10 70 "$SCRIPT2_TEs_file" 3>&1 1>&2 2>&3 || echo "$SCRIPT2_TEs_file")
    ;;
esac
  done
}

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
    SCRIPT1_annotation_file="$SCRIPT1_DEFAULT_annotation_file"
    SCRIPT1_description_file="$SCRIPT1_DEFAULT_description_file"
    SCRIPT1_TEs_file="$SCRIPT1_DEFAULT_TEs_file"

    # Directly go to the parameters selection menu
    edit_script1_parameters || exit 1
fi

###################
# Gather Methylome.At_metaPlots.sh
###################
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "metaPlots" ]]; then
    # Initialize parameters with defaults
    SCRIPT2_Genes_n_TEs="$SCRIPT2_DEFAULT_Genes_n_TEs"
    SCRIPT2_Gene_features="$SCRIPT2_DEFAULT_Gene_features"
    SCRIPT2_minReadsPerCytosine="$SCRIPT2_DEFAULT_minReadsPerCytosine"
    SCRIPT2_metaPlot_random_genes="$SCRIPT2_DEFAULT_metaPlot_random_genes"
    SCRIPT2_n_cores="$SCRIPT2_DEFAULT_n_cores"
    SCRIPT2_bin_size_features="$SCRIPT2_DEFAULT_bin_size_features"
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
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Methylome.At" ]]; then
    chosen_message+="'Methylome.At' "
fi
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "metaPlots" ]]; then
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
  # User selected Yes, proceed with running the scripts
  # Methylome.At.sh invocation
  if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Methylome.At" ]]; then
    echo "Running Methylome.At.sh..."
    bash "$SCRIPT1_PATH" \
      --samples_file "$SAMPLES_FILE" \
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
      --annotation_file "$SCRIPT1_annotation_file" \
      --description_file "$SCRIPT1_description_file" \
      --TEs_file "$SCRIPT1_TEs_file"
  fi

  # Methylome.At_metaPlots.sh invocation
  if [[ " ${SELECTED_SCRIPTS[*]} " =~ "metaPlots" ]]; then
    echo "Running Methylome.At_metaPlots.sh..."
    bash "$SCRIPT2_PATH" \
      --samples_file "$SAMPLES_FILE" \
      --Genes_n_TEs "$SCRIPT2_Genes_n_TEs" \
      --Gene_features "$SCRIPT2_Gene_features" \
      --minReadsPerCytosine "$SCRIPT2_minReadsPerCytosine" \
      --metaPlot_random_genes "$SCRIPT2_metaPlot_random_genes" \
      --n_cores "$SCRIPT2_n_cores" \
      --bin_size_features "$SCRIPT2_bin_size_features" \
      --annotation_file "$SCRIPT2_annotation_file" \
      --TEs_file "$SCRIPT2_TEs_file"
  fi

  echo "All done!"
else
  # User selected No, just exit without running
  echo "You chose not to run the scripts. Exiting."
  exit 0
fi
