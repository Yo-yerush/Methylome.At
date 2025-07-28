#!/bin/bash

# Initialize conda
eval "$(conda shell.bash hook)"

# Determine the directory of 'Methylome.At'
script_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
cd "$script_dir"

# Generate log file with a timestamp
log_file="setup_env.log"
echo "**  $(date +"%d-%m-%y %H:%M")" > "$log_file"
echo "" >> "$log_file"

# Create and activate the conda environment
packages=("r-curl" "r-rcurl" "zlib" "r-textshaping" "harfbuzz" "fribidi" "freetype" "libpng" "pkg-config" "libxml2" "r-xml" "bioconductor-rsamtools") # "r-devtools"
echo "Creating Conda environment..." >> "$log_file"
if conda create --name Methylome.At_env -c conda-forge -c bioconda r-base=4.4.3 "${packages[@]}" -y && conda activate Methylome.At_env; then
    echo "Conda environment created and activated successfully." >> "$log_file"
    echo "Conda environment name: $CONDA_DEFAULT_ENV" >> "$log_file"
else
    echo "Failed to create and activate Conda environment." >> "$log_file"
    exit 1
fi

echo "" >> "$log_file"

# Check if packages are installed from conda-forge
echo "Verifying that packages are installed from conda-forge..." >> "$log_file"
error=0
for pkg in "r-base" "${packages[@]}"; do
    # Get the package info
    pkg_info=$(conda list "$pkg")
    # Extract the channel from the package info
    channel=$(echo "$pkg_info" | awk '/^'"$pkg"'[[:space:]]/ {print $NF}')
    if [ "$channel" = "conda-forge" ] || [ "$channel" = "bioconda" ]; then
        echo "* installed $pkg: yes"
    else
        echo "* installed $pkg: no"
    fi
done

echo "" >> "$log_file"

# Run the R script to install additional R packages
echo "Install R packages..." >> "$log_file"
Rscript scripts/install_R_packages.R 2>> "$log_file"

# Permission
chmod +x ./Methylome.At.sh
chmod +x ./Methylome.At_metaPlots.sh
chmod +x ./Methylome.At_UI.sh
chmod +x scripts/*.sh
