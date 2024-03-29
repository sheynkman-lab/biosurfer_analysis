#!/bin/sh
#Author: Mayank Murali
#Project: Biosurfer

#Script to download GENCODE toy files from Zenodo (https://zenodo.org/record/7297008)

echo "================================================================="
echo "  Downloading GENCODE toy data ..."
echo "================================================================="

cd data
mkdir A_gencode_toy
cd A_gencode_toy

wget https://zenodo.org/records/10822882/files/biosurfer_gencode_toy_data.zip
unzip biosurfer_gencode_toy_data.zip 
rm -rf __MACOSX biosurfer_gencode_toy_data.zip
