#!/bin/bash

# Update and install essential tools
sudo apt-get update -y
sudo apt-get install -y build-essential wget unzip git python3 python3-pip openjdk-11-jdk

# Install BWA
sudo apt-get install -y bwa

# Install SAMtools
sudo apt-get install -y samtools

# Install GATK
GATK_VERSION="4.3.0.0"
wget https://github.com/broadinstitute/gatk/releases/download/$GATK_VERSION/gatk-$GATK_VERSION.zip
unzip gatk-$GATK_VERSION.zip
sudo mv gatk-$GATK_VERSION /usr/local/gatk
export PATH="/usr/local/gatk:$PATH"

# Install CNVkit
pip3 install cnvkit

# Install FastQC
sudo apt-get install -y fastqc

# Install MultiQC
pip3 install multiqc

# Install Trimmomatic
TRIMMOMATIC_VERSION="0.39"
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/$TRIMMOMATIC_VERSION/Trimmomatic-$TRIMMOMATIC_VERSION.zip
unzip Trimmomatic-$TRIMMOMATIC_VERSION.zip
sudo mv Trimmomatic-$TRIMMOMATIC_VERSION /usr/local/Trimmomatic

# Add Trimmomatic to PATH
echo 'export PATH="/usr/local/Trimmomatic:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Install Snakemake
pip3 install snakemake

# Install ANNOVAR (Requires manual registration to download)
# Please download from https://annovar.openbioinformatics.org/en/latest/user-guide/download/
# After downloading, move it to /usr/local/annovar and add to PATH
# sudo mv annovar /usr/local/annovar
# echo 'export PATH="/usr/local/annovar:$PATH"' >> ~/.bashrc
# source ~/.bashrc

echo "All tools installed successfully!"
