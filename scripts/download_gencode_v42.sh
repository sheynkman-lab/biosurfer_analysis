#!/bin/sh
#Author: Mayank Murali
#Project: Biosurfer

#Script to download GENCODE v42 files from Zenodo (https://zenodo.org/record/7297008)

echo "================================================================="
echo "  Downloading GENCODE v42 data ..."
echo "================================================================="

cd data
mkdir A_gencode_v42
cd A_gencode_v42

wget https://zenodo.org/records/10822882/files/biosurfer_gencode_v42_data.zip
unzip biosurfer_gencode_v42_data.zip 
rm -rf __MACOSX biosurfer_gencode_v42_data.zip
