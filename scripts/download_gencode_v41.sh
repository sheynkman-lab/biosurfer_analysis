#!/bin/sh
#Author: Mayank Murali
#Project: Biosurfer

#Script to download GENCODE v41 files from Zenodo (https://zenodo.org/record/7182809)

echo "================================================================="
echo "  Downloading GENCODE v41 data ..."
echo "================================================================="

cd data
mkdir A_gencode_v41
cd A_gencode_v41

wget https://zenodo.org/record/7182809/files/biosurfer_gencode_v41_data.zip
unzip biosurfer_gencode_v41_data.zip 
rm -rf __MACOSX biosurfer_gencode_v41_data.zip
