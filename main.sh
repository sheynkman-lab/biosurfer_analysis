# Steps to run the biosurfer analysis scripts

# Step 1 - create the conda environment for biosurfer

# bash ./scripts/create_conda_env.sh
conda create --name biosurfer-install --channel conda-forge python=3 pip 
conda install graph-tool

# Then activate the environment:
conda activate biosurfer-install

# Step 2 - install biosurfer (the stable dev-version)
bash ./scripts/install_biosurfer.sh

# Step 3 - download GENCODE toy data
bash ./scripts/download_gencode_toy.sh 

# Step 4 - Create toy dataset on biosurfer
cd biosurfer # Have to cd in to biosurfer directory to run biosurfer commands

biosurfer load_db \
 --source=GENCODE \
 --gtf ../A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.gtf \
 --tx_fasta ../A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.transcripts.fa \
 --tl_fasta ../A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.translations.fa \
 -d gencode_toy

# Step 5 - Run biosurfer hybrid alignment on the toy dataset
mkdir ../B_hybrid_aln_results_toy
biosurfer hybrid_alignment -d gencode_toy -o ../B_hybrid_aln_results_toy 

# Step 6 - Create example visualization of isoforms using biosurfer
biosurfer plot -d gencode_toy --gene CRYBG2


## Steps from here will be based on Jupyter notebook (ipynb) for plotting

#Step 7 - Post-process hybrid alignment output table
cd .. # cd out of biosurfer directory
mkdir C_toy_plots

#Install required libraries 
pip install xlsxwriter
pip install openpyxl
pip install plotly

# If terminal command based
ipython
run ./scripts/hybrid_aln_table_processing.ipynb

#Else on ipynb through website (pip install notebook)
jupyter notebook ./scripts/hybrid_aln_table_processing.ipynb




