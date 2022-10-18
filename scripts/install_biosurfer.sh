%%writefile create_conda_env.sh
#!/usr/bin/env bash
#Author: Mayank Murali
#Project: Biosurfer

#Script to download and install Biosufer from GitHub

# Clone the repository
git clone -b dev --single-branch https://github.com/sheynkman-lab/biosurfer.git
    
# Move to the folder
cd biosurfer
    
# Run setup 
pip install --editable .

