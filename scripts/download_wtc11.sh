#!/bin/sh
#Author: Mayank Murali
#Project: Biosurfer

#Script to download PacBio WTC11 from Zenodo (https://zenodo.org/record/7182809)

echo "================================================================="
echo "  Downloading PacBio WTC11 data ..."
echo "================================================================="

cd data
mkdir A_wtc11
cd A_wtc11

wget https://zenodo.org/record/7182809/files/biosurfer_wtc11_data.zip
unzip biosurfer_wtc11_data.zip 
rm -rf __MACOSX biosurfer_wtc11_data.zip