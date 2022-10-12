%%writefile create_conda_env.sh
#!/usr/bin/env bash
#Author: Mayank Murali
#Project: Biosurfer

#Script to create a new Conda environment

# User input to create a new conda env
# read -p "Create new conda env (y/n)?" CONT

# if [ "$CONT" == "n" ]; then
#   echo "exit";
# else

# # Prompt user for conda env name
#   echo "Creating new conda environment, choose name"
#   read input_variable
#   echo "Name $input_variable was chosen";

  # Create environment.yml or not
  # read -p "Create 'enviroment.yml', will overwrite if exist (y/n)?"
  #   if [ "$CONT" == "y" ]; then
  #     # yes: create enviroment.yml
echo "# BASH: conda env create

# source activate phd
name: bio_env
dependencies:
- pip: ">environment.yml    
    
  #list name of packages
conda env create
    # else
    #     echo "installing base packages"
    #     conda create --name bio_env\
    #     python=3 jupyter notebook numpy rpy2\
    #     pandas scipy numpy scikit-learn seaborn 
    # fi
echo "to exit: source deactivate"
# fi

# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate $input_variable
# conda init
# conda activate bio_env