#!/bin/zsh
# Run add-halpha.py for multiple planet spectra and H-alpha templates

set -e

# List of planet spectra
planets=(
  "2MJ1612b_BT-Settl-CIFIST-1200K-3.5logg_140000.txt"
  "PDS70b_BT-Settl-CIFIST-1400K-4logg_140000.txt"
  "PDS70c_BT-Settl-CIFIST-1300K-4logg_140000.txt"
  "WISPIT2b_BT-Settl-CIFIST-1400K-4logg_140000.txt"
)

# List of H-alpha profiles
halpha_files=(
  "Ha_60_14.dat"
  "Ha_60_12.dat"
  "Ha_80_14.dat"
  "Ha_80_12.dat"
)

echo "Starting add-halpha runs..."
for planet in "${planets[@]}"; do
  for ha in "${halpha_files[@]}"; do
    echo "Processing $planet with $ha..."
    python add-halpha.py "$planet" "$ha"
  done
done

echo "All add-halpha runs complete!"
