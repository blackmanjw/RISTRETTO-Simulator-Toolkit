#!/bin/zsh

# Directory containing the FITS files
DATA_DIR="pychelle_output/onaxis"

# Loop over each .fits file in the directory
for fits_file in "$DATA_DIR"/*.fits; do
    # Get the filename without the directory and extension
    filename=$(basename "$fits_file" .fits)

    echo "Processing $filename ..."

    # Run analysis.py for both indices
    python analysis.py "onaxis/$filename" 1
    python analysis.py "onaxis/$filename" 2
done
