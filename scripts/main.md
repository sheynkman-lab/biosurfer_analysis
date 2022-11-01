# Biosurfer Analysis

#### Steps to run the biosurfer analysis that reproduces results, summary plots and figures for Biosurfer manuscript.

## Step 1 - Create the conda environment for biosurfer

```
conda create --name biosurfer-install --channel conda-forge python=3 pip 
```

### Then activate the environment:
```
conda activate biosurfer-install

conda install --channel conda-forge graph-tool
```


---

## Step 2 - install biosurfer (the stable dev-version)

### Clone the repository
```
git clone -b dev --single-branch https://github.com/sheynkman-lab/biosurfer.git
```    
### Run setup 
```
pip install --editable biosurfer
```


---

## Step 3 - download input data
```
for source in gencode_toy gencode_v41 wtc11
do
    bash "./scripts/download_$source.sh"
done
```


---

## Step 4 - Create toy dataset on biosurfer
    
### gencode_toy
```
biosurfer load_db \
    --source=GENCODE \
    --gtf A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.gtf \
    --tx_fasta A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.transcripts.fa \
    --tl_fasta A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.translations.fa \
    -d gencode_toy
```
### gencode_v42
```
biosurfer load_db \
    --source=GENCODE \
    --gtf A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.basic.annotation.gtf \
    --tx_fasta A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.pc_transcripts.fa \
    --tl_fasta A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.pc_translations.fa \
    -d gencode_v42

```
### wtc11
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

## Step 5 - Run biosurfer hybrid alignment on the toy dataset
### gencode_toy
```
mkdir B_hybrid_aln_results_toy
biosurfer hybrid_alignment \
    -d gencode_toy \
    -o B_hybrid_aln_results_toy \
    --gencode
```
### gencode_v42
```    
mkdir B_hybrid_aln_gencode_v42
biosurfer hybrid_alignment \
    -d gencode_v42 \
    -o B_hybrid_aln_gencode_v42 \
    --gencode
```
### wtc11
```
mkdir B_hybrid_aln_wtc11
biosurfer hybrid_alignment \
    -d wtc11 \
    -o B_hybrid_aln_wtc11
```

---


## Step 6 - Post-process hybrid alignment output table / genome wide analysis


### Install required libraries 
```
pip install ipykernel xlsxwriter openpyxl plotly
```

### Run Python script
```
python3 ./scripts/genome_wide_summary.py
```

## Step 7 - N-termini summary 

### Run Python script
```
python3 ./scripts/n_termini_summary.py
```

## Step 8 - C-termini summary

### Run Python script
```
python3 ./scripts/c_termini_summary.py
```

## Step 9 - Internal region summary

### Run Python script
```
python3 ./scripts/internal_summary.py
```
## Step 10 - Isoform plot visualization using Biosurfer
```
bash ./scripts/isoform_plotting.sh
```
