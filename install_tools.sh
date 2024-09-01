#!/bin/bash

# Update and install essential tools
sudo apt-get update -y
sudo apt-get install -y build-essential wget unzip git python3 python3-pip openjdk-11-jdk curl

# Install BWA
echo "Installing BWA..."
sudo apt-get install -y bwa

# Install SAMtools
echo "Installing SAMtools..."
sudo apt-get install -y samtools

# Install GATK
echo "Installing GATK..."
GATK_VERSION="4.6.0.0"
wget https://github.com/broadinstitute/gatk/releases/download/$GATK_VERSION/gatk-$GATK_VERSION.zip
unzip gatk-$GATK_VERSION.zip
sudo mv gatk-$GATK_VERSION /usr/local/gatk
echo 'export PATH="/usr/local/gatk:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Install CNVkit
echo "Installing CNVkit..."
pip3 install cnvkit

# Install FastQC
echo "Installing FastQC..."
sudo apt-get install -y fastqc

# Install MultiQC
echo "Installing MultiQC..."
pip3 install multiqc

# Install Trimmomatic
echo "Installing Trimmomatic..."
TRIMMOMATIC_VERSION="0.4"
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/$TRIMMOMATIC_VERSION/Trimmomatic-$TRIMMOMATIC_VERSION.zip
unzip Trimmomatic-$TRIMMOMATIC_VERSION.zip
sudo mv Trimmomatic-$TRIMMOMATIC_VERSION /usr/local/Trimmomatic
echo 'export PATH="/usr/local/Trimmomatic:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Install Snakemake
echo "Installing Snakemake..."
pip3 install snakemake

# Install HTSlib (SAMtools dependency)
echo "Installing HTSlib..."
HTSLIB_VERSION="1.2"
wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2
tar -xvjf htslib-$HTSLIB_VERSION.tar.bz2
cd htslib-$HTSLIB_VERSION
./configure
make
sudo make install
cd ..

# Install ANNOVAR (Requires manual download)
echo "Installing ANNOVAR (requires manual download)..."
# Note: ANNOVAR requires a manual download from the official site
# Please download it from https://annovar.openbioinformatics.org/en/latest/user-guide/download/
# After downloading, move it to /usr/local/annovar and add to PATH
# sudo mv annovar /usr/local/annovar
# echo 'export PATH="/usr/local/annovar:$PATH"' >> ~/.bashrc
# source ~/.bashrc

# Clean up
rm -rf gatk-$GATK_VERSION.zip Trimmomatic-$TRIMMOMATIC_VERSION.zip htslib-$HTSLIB_VERSION.tar.bz2 htslib-$HTSLIB_VERSION

echo "All tools installed successfully!"
