#!/usr/bin/bash
#
# Script to bootstrap the enrivonment for members of the barkanlab pirg
#

# Add shared scripts to user path
echo "PATH=$PATH:/projects/barkanlab/shared/bin" >> ~/.bash_profile
echo "export PATH" >> ~/.bash_profile

# Load needed module
echo "module load python3" >> ~/.bash_profile

# Move to project folder
echo "cd /projects/barkanlab/$USER" >> ~/.bash_profile

# Setup permissions and link to the illumina data location
chmod -R g+rwxs /projects/barkanlab/$USER
ln -s /projects/barkanlab/shared/illumina_data /projects/barkanlab/$USER/illumina_data

# Load new environment
source ~/.bash_profile
