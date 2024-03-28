#!/usr/bin/env bash

# I opted for Homebrew since I encountered difficulty 
finding suitable Conda channels with Homebrew 

# Install wget for downloading the sequence files
brew install wget

# Install FastQC for quality control
brew install FastQC

# Install Trimmomatic for trimming
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
mv Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar

# Install bwa for genome mapping
brew install bwa

# Install samtools for working with SAM/BAM files
brew install samtools

# Install freebayes for variant calling
brew install freebayes
