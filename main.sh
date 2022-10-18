# Steps to run the biosurfer analysis scripts

# Step 1 - create the conda environment for biosurfer

# bash ./scripts/create_conda_env.sh
conda create --name biosurfer-install --channel conda-forge python=3 pip 

# Then activate the environment:
conda activate biosurfer-install

conda install --channel conda-forge graph-tool

# Step 2 - install biosurfer (the stable dev-version)
#bash ./scripts/install_biosurfer.sh

    # Clone the repository
    git clone -b dev --single-branch https://github.com/sheynkman-lab/biosurfer.git
        
    # Run setup 
    pip install --editable biosurfer

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
biosurfer hybrid_alignment -d gencode_toy -o ../B_hybrid_aln_results_toy --gencode


## Steps from here until step 10 will be based on Jupyter notebook (ipynb) for plotting

#Step 6 - Post-process hybrid alignment output table / genome wide analysis
cd .. # cd out of biosurfer directory
mkdir C_toy_plots

#Install required libraries 
pip install xlsxwriter openpyxl plotly

    # If terminal command based
    ipython
    run ./scripts/genome_wide_summary.ipynb

    #Else on ipynb through website (pip install notebook)
    jupyter notebook ./scripts/genome_wide_summary.ipynb


# Step 7 - N-termini summary 

mkdir D_nterm_plots

    # If terminal command based
    ipython
    run ./scripts/n_termini_summary.ipynb

    #Else on ipynb through website (pip install notebook)
    jupyter notebook ./scripts/n_termini_summary.ipynb

# Step 8 - C-termini summary
mkdir E_cterm_plots

    # If terminal command based
    ipython
    run ./scripts/c_termini_summary.ipynb

    #Else on ipynb through website (pip install notebook)
    jupyter notebook ./scripts/c_termini_summary.ipynb


# Step 9 - Internal region summary
mkdir F_internal_region_plots

    # If terminal command based
    ipython
    run ./scripts/internal_region_summary.ipynb

    #Else on ipynb through website (pip install notebook)
    jupyter notebook ./scripts/internal_region_summary.ipynb

# Step 10 - Isoform plot visualization using Biosurfer
bash ./scripts/isoform_plotting.sh

