#!/bin/bash

# Initialize conda
eval "$(conda shell.bash hook)"

# Generate log file with a timestamp
log_file="setup_env.log"
echo "**  $(date +"%d-%m-%y %H:%M")" > "$log_file"
echo "" >> "$log_file"

# Check R version
#R_version=$(R --version | head -n1 | awk '{print $3}')
#if [[ "$R_version" < "4.4" ]]; then
#echo "R version is '<4.4.0', install '4.4.2'" >> "$log_file"
#    packages=("r-base=4.4.2" "r-curl" "r-rcurl" "r-devtools" "zlib" "r-textshaping" "harfbuzz" "fribidi" "freetype" "libpng" "pkg-config")
#else
#echo "R is already istalled (version '$R_version')" >> "$log_file"
#    packages=("r-curl" "r-rcurl" "r-devtools" "zlib" "r-textshaping" "harfbuzz" "fribidi" "freetype" "libpng" "pkg-config") # "r-textshaping"
#fi
#echo "" >> "$log_file"

packages=("r-curl" "r-rcurl" "r-devtools" "zlib" "r-textshaping" "harfbuzz" "fribidi" "freetype" "libpng" "pkg-config")

# Create and activate the conda environment
echo "Creating Conda environment..." >> "$log_file"
if conda create --name Methylome.At_env -c conda-forge r-base=4.4.2 "${packages[@]}" -y && conda activate Methylome.At_env; then
    echo "Conda environment created and activated successfully." >> "$log_file"
    echo "Conda environment name: $CONDA_DEFAULT_ENV" >> "$log_file"
else
    echo "Failed to create and activate Conda environment." >> "$log_file"
    exit 1
fi

echo "" >> "$log_file"

# Check if R packages are installed from conda-forge
echo "Verifying that packages are installed from conda-forge..." >> "$log_file"
error=0
for pkg in "r-base" "${packages[@]}"; do
    # Get the package info
    pkg_info=$(conda list "$pkg")
    # Extract the channel from the package info
    channel=$(echo "$pkg_info" | awk '/^'"$pkg"'[[:space:]]/ {print $NF}')
    if [ "$channel" == "conda-forge" ]; then
        echo "* installed $pkg: yes" >> "$log_file"
    else
        echo "* installed $pkg: no" >> "$log_file"
        #exit 1
    fi
done

echo "" >> "$log_file"

# Run the R script to install additional R packages
echo "Install R packages..." >> "$log_file"
Rscript scripts/install_R_packages.R 2>> "$log_file"

# Determine the directory of 'Methylome.At'
$script_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
cd "$script_dir"

# Permission
chmod +x ./Methylome.At.sh
chmod +x ./Methylome.At_metaPlots.sh
chmod +x ./Methylome.At_UI.sh
chmod +x scripts/cgmap_to_cx.sh