# Steps to run the biosurfer analysis scripts

# Step 1 - create the conda environment for biosurfer

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

# Step 3 - download input data
for source in gencode_toy gencode_v41 wtc11
do
    bash "./scripts/download_$source.sh"
done

# Step 4 - Create toy dataset on biosurfer
    # gencode_toy
    biosurfer load_db \
    --source=GENCODE \
    --gtf A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.gtf \
    --tx_fasta A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.transcripts.fa \
    --tl_fasta A_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.translations.fa \
    -d gencode_toy

    # gencode_v41
    biosurfer load_db --source=GENCODE --gtf A_gencode_v41/biosurfer_gencode_v41_data/gencode.v41.annotation.gtf --tx_fasta A_gencode_v41/biosurfer_gencode_v41_data/gencode.v41.pc_transcripts.fa --tl_fasta A_gencode_v41/biosurfer_gencode_v41_data/gencode.v41.pc_translations.fa -d gencode_v41

    #gencode_v42
    biosurfer load_db --source=GENCODE --gtf A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.basic.annotation.gtf --tx_fasta A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.pc_transcripts.fa --tl_fasta A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.pc_translations.fa -d gencode_v42


    # wtc11 with APPRIS from gencode_41
    biosurfer load_db --source=GENCODE --gtf A_gencode_v41/biosurfer_gencode_v41_data/gencode.v41.annotation.gtf --tx_fasta A_gencode_v41/biosurfer_gencode_v41_data/gencode.v41.pc_transcripts.fa --tl_fasta A_gencode_v41/biosurfer_gencode_v41_data/gencode.v41.pc_translations.fa -d wtc11_db
    biosurfer load_db --source=PacBio --gtf A_wtc11/biosurfer_wtc11_data/wtc11_with_cds.gtf --tx_fasta A_wtc11/biosurfer_wtc11_data/wtc11_corrected.fasta --tl_fasta A_wtc11/biosurfer_wtc11_data/wtc11_orf_refined.fasta --sqanti A_wtc11/biosurfer_wtc11_data/wtc11_classification.txt -d wtc11_db

    # wtc11 with APPRIS from gencode_42
    biosurfer load_db --source=GENCODE --gtf A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.basic.annotation.gtf --tx_fasta A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.pc_transcripts.fa --tl_fasta A_gencode_v42/biosurfer_gencode_v42_data/gencode.v42.pc_translations.fa -d wtc11_db_v42
    biosurfer load_db --source=PacBio --gtf A_wtc11/biosurfer_wtc11_data/wtc11_with_cds.gtf --tx_fasta A_wtc11/biosurfer_wtc11_data/wtc11_corrected.fasta --tl_fasta A_wtc11/biosurfer_wtc11_data/wtc11_orf_refined.fasta --sqanti A_wtc11/biosurfer_wtc11_data/wtc11_classification.txt -d wtc11_db_v42
    biosurfer load_db --source=PacBio --gtf A_wtc11/biosurfer_wtc11_data/wtc11_with_cds.gtf --tx_fasta A_wtc11/biosurfer_wtc11_data/wtc11_corrected.fasta --tl_fasta A_wtc11/biosurfer_wtc11_data/wtc11_orf_refined.fasta --sqanti A_wtc11/biosurfer_wtc11_data/wtc11_classification.txt -d wtc11_db

# Step 5 - Run biosurfer hybrid alignment on the toy dataset
    # gencode_toy
    mkdir B_hybrid_aln_results_toy
    biosurfer hybrid_alignment -d gencode_toy -o B_hybrid_aln_results_toy --gencode

    # gencode_v41
    mkdir B_hybrid_aln_gencode_v41
    biosurfer hybrid_alignment -d gencode_v41 -o B_hybrid_aln_gencode_v41 --gencode

    # gencode_v42
    mkdir B_hybrid_aln_gencode_v42
    biosurfer hybrid_alignment -d gencode_v42 -o B_hybrid_aln_gencode_v42 --gencode

    # wtc11
    mkdir B_hybrid_aln_wtc11
    biosurfer hybrid_alignment -d wtc11 -o B_hybrid_aln_wtc11

    # wtc11_v42
    mkdir B_hybrid_aln_wtc11_v42
    biosurfer hybrid_alignment -d wtc11_db_v42 -o B_hybrid_aln_wtc11_v42 

## Steps from here until step 10 will be based on Jupyter notebook (ipynb) for plotting

#Step 6 - Post-process hybrid alignment output table / genome wide summary plots

# Run this py script on terminal
python3 python3 ./scripts/genome_wide_summary.py


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

