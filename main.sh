# Steps to run the biosurfer analysis scripts

# Step 1 - create the conda environment for biosurfer

# bash ./scripts/create_conda_env.sh
conda create --name bio

# Then activate the environment:
conda activate bio_env

# Step 2.1 - install pip and graph-tool dependency through conda
conda install -c conda-forge graph-tool
coda install pip

# Step 2 - install biosurfer
bash ./scripts/install_biosurfer.sh


# Step 3 - download GENCODE toy data

bash ./scripts/download_gencode_toy.sh 

# Step 4 - run biosurfer on toy dataset

cd biosurfer

# create the database
biosurfer load_db \
 --source=GENCODE \
 --gtf ../3_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.gtf \
 --tx_fasta ../3_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.transcripts.fa \
 --tl_fasta ../3_gencode_toy/biosurfer_gencode_toy_data/gencode.v38.toy.translations.fa \
 -d gencode_toy

# Step 5 - run hybrid alignment on the toy dataset
mkdir ../5_hybrid_aln_results_toy
biosurfer hybrid_alignment -d gencode_toy -o ../5_hybrid_aln_results_toy --summary

# Step 6 - visualize the hybrid alignment results
# Note that the results will be written out in `5_hybrid_aln_results_toy`
biosurfer plot -d gencode_toy --gene CRYBG2

