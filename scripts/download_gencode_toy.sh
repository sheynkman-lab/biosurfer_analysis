#!/bin/sh
#Author: Mayank Murali
#Project: Biosurfer

#Script to download GENCODE toy files from Zenodo (https://zenodo.org/record/7182809)

echo "================================================================="
echo "  Downloading GENCODE toy data ..."
echo "================================================================="

cd data
mkdir 3_gencode_toy
cd 3_gencode_toy

wget https://zenodo.org/record/7182809/files/biosurfer_gencode_toy_data.zip
unzip biosurfer_gencode_toy_data.zip 
rm -rf __MACOSX biosurfer_gencode_toy_data.zip
