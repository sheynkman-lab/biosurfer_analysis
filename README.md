# Biosurfer Analysis

## Analysis accompanying the manuscript *"Biosurfer: Connecting genomic, transcriptomic, and proteomic information layers to track mechanisms of protein isoform variation"*

This repository contains steps to run the biosurfer analysis, which reproduces the results, summary plots, and figures for the Biosurfer manuscript ([bioRxiv]()).


**Contents**
1. [Download and install Biosurfer](#download-and-install-biosurfer)
2. [Download input data](#download-input-data)
3. [Run Biosurfer modules](#run-biosurfer-modules)
    i.  [Load database](#load-database)
    ii.  [Run hybrid alignment](#run-hybrid-alignment)
    iii.  [Visualize protein isoforms](#visualize-protein-isoforms)
4. [Global characterization of altered protein regions in the human annotation (GENCODE)](#post-processing)
    i. [Altered protein regions across the human proteome](#genome-wide-summary)
    ii. [Analysis of alternative splicing events that alter the N-terminus of proteins](#n-term)
    iii. [Characterization of splicing patterns underlying internal protein region differences](#internal-region)
    iv. [Analyzing splicing patterns for C-terminal alterations](#c-term)
    


<a id="download-and-install-biosurfer"></a>
## 1. Download and install Biosurfer


#### Create the conda environment for Biosurfer
```
conda create --name biosurfer-install --channel conda-forge python=3 pip 
```
#### Activate the conda environment:
```
conda activate biosurfer-install

conda install --channel conda-forge graph-tool
```

#### Install Biosurfer (stable version): Clone the repository
```
git clone -b dev --single-branch https://github.com/sheynkman-lab/biosurfer.git
```    
#### Run setup 
```
pip install --editable biosurfer
```
---

<a id="download-input-data"></a>
## 2. Download input data

```
for source in gencode_toy gencode_v41 wtc11
do
    bash "./scripts/download_$source.sh"
done
```
---

<a id="run-biosurfer-modules"></a>
## 3. Run Biosurfer modules

<a id="load-database"></a>
### i. Load database
    
#### gencode_toy
```
biosurfer load_db \
    --source=GENCODE \
    --gtf A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.gtf \
    --tx_fasta A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.transcripts.fa \
    --tl_fasta A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.translations.fa \
    -d gencode_toy
```
#### gencode_v42
```
biosurfer load_db \
    --source=GENCODE \
    --gtf A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.basic.annotation.gtf \
    --tx_fasta A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.pc_transcripts.fa \
    --tl_fasta A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.pc_translations.fa \
    -d gencode_v42
```
#### wtc11
```
biosurfer load_db \
    --source=GENCODE \
    --gtf A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.basic.annotation.gtf \
    --tx_fasta A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.pc_transcripts.fa \
    --tl_fasta A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.pc_translations.fa \
    -d wtc11
```
```    
    biosurfer load_db \
    --source=PacBio \
    --gtf A_wtc11/biosurfer_wtc11_data/wtc11_with_cds.gtf \
    --tx_fasta A_wtc11/biosurfer_wtc11_data/wtc11_corrected.fasta \
    --tl_fasta A_wtc11/biosurfer_wtc11_data/wtc11_orf_refined.fasta \
    --sqanti A_wtc11/biosurfer_wtc11_data/wtc11_classification.txt \
    -d wtc11
```
---

<a id="run-hybrid-alignment"></a>
### ii. Run hybdrid alignment

#### gencode_toy
```
mkdir B_hybrid_aln_results_toy
biosurfer hybrid_alignment \
    -d gencode_toy \
    -o B_hybrid_aln_results_toy \
    --gencode
```
#### gencode_v42
```    
mkdir B_hybrid_aln_gencode_v42
biosurfer hybrid_alignment \
    -d gencode_v42 \
    -o B_hybrid_aln_gencode_v42 \
    --gencode
```
#### wtc11
```
mkdir B_hybrid_aln_wtc11
biosurfer hybrid_alignment \
    -d wtc11 \
    -o B_hybrid_aln_wtc11
```
---

<a id="visualize-protein-isoforms"></a>
### iii. Visualize protein isoforms

```
bash ./scripts/isoform_plotting.sh
```
---

<a id="post-processing"></a>
## 4. Global characterization of altered protein regions in the human annotation (GENCODE)


#### Install required libraries 
```
pip install ipykernel xlsxwriter openpyxl plotly
```

<a id="genome-wide-summary"></a>
### i. Altered protein regions across the human proteome
```
python3 ./scripts/genome_wide_summary.py
```

<a id="n-term"></a>
### ii. Analysis of alternative splicing events that alter the N-terminus of proteins
```
python3 ./scripts/n_termini_summary.py
```

<a id="internal-region"></a>
### iii. Characterization of splicing patterns underlying internal protein region differences
```
python3 ./scripts/internal_summary.py
```

<a id="c-term"></a>
### iv. Analyzing splicing patterns for C-terminal alterations
```
python3 ./scripts/c_termini_summary.py
```


