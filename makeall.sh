#!/bin/zsh
# Run all makecsv scripts in batch mode (Zsh version)

# Exit immediately if a command fails
set -e

echo "Starting batch processing..."

# Stars
python makecsv.py 2MJ1612_star_BT-Settl-CIFIST-3900K-4logg_140000.txt
python makecsv.py PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt
python makecsv.py WISPIT2_star_BT-Settl-CIFIST-4400K-4logg_140000.txt

python makecsv_batch.py 2MJ1612_star_BT-Settl-CIFIST-3900K-4logg_140000.txt --noshow
python makecsv_batch.py PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt --noshow
python makecsv_batch.py WISPIT2_star_BT-Settl-CIFIST-4400K-4logg_140000.txt --noshow

# Planets (repeat for multiple Ha profiles)
halphas=(
  "Ha_60_12"
  "Ha_60_14"
  "Ha_80_12"
  "Ha_80_14"
)

for ha in "${halphas[@]}"; do
  echo "Processing planets for ${ha}..."

  python makecsv.py 2MJ1612b_BT-Settl-CIFIST-1200K-3.5logg_140000_${ha}.txt
  python makecsv.py PDS70b_BT-Settl-CIFIST-1400K-4logg_140000_${ha}.txt
  python makecsv.py PDS70c_BT-Settl-CIFIST-1300K-4logg_140000_${ha}.txt
  python makecsv.py WISPIT2b_BT-Settl-CIFIST-1400K-4logg_140000_${ha}.txt

done

echo "All processing complete!"
