#!/bin/bash

# Initialize conda
eval "$(conda shell.bash hook)"

# Create and activate the conda environment
echo "Creating Conda environment..."
if conda create --name Methylome.At_env -c conda-forge r-base r-curl r-rcurl r-textshaping -y && conda activate Methylome.At_env; then
    echo "Conda environment created and activated successfully."
else
    echo "Failed to create and activate Conda environment."
    exit 1
fi

# Check if R packages are installed from conda-forge
echo "Verifying that packages are installed from conda-forge..."
packages=("r-base" "r-curl" "r-rcurl" "r-textshaping")
error=0
for pkg in "${packages[@]}"; do
    # Get the package info
    pkg_info=$(conda list "$pkg")
    # Extract the channel from the package info
    channel=$(echo "$pkg_info" | awk '/^'"$pkg"'[[:space:]]/ {print $NF}')
    if [ "$channel" == "conda-forge" ]; then
        echo "$pkg is installed from conda-forge."
    else
        echo "Error: $pkg is not installed from conda-forge. Installed from $channel."
        exit 1
    fi
done

# Run the R script to install additional R packages
echo "Install R packages..."
Rscript scripts/install_R_packages.R

# Determine the directory of 'Methylome.At'
$script_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
cd "$script_dir"