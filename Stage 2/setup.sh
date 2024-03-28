#!/usr/bin/env bash

# I opted for Homebrew since I encountered difficulty 
finding suitable Conda channels with Homebrew 

# Install wget for downloading the sequence files
brew install wget

# Install FastQC for quality control
brew install FastQC

# Install MultiQC
pip install multiqc

# Install Trimmomatic for trimming and filtering
# Create a directory to store Trimmomatic
mkdir -p ~/bin/trimmomatic

# Navigate to the directory
cd ~/bin/trimmomatic

# Download Trimmomatic JAR file
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip

# Unzip the downloaded file
unzip Trimmomatic-0.39.zip

# Remove the zip file (optional)
rm Trimmomatic-0.39.zip

# Return to home directory
cd

# Install bwa for genome mapping
brew install bwa

# Install samtools for working with SAM/BAM files
brew install samtools

# Install freebayes for variant calling
brew install freebayes
