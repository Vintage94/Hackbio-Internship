#!/usr/bin/env bash

# Make a new directory for Project 3 and change to Project 3 directory
mkdir Project3 && cd "$_"

# Download the forward sequence
wget "https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz?download=1"

# Download the reverse sequence
wget "https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz?download=1"

# Rename the files
mv ERR8774458_1.fastq.gz?download=1 ERR8774458_1.fastq.gz
mv ERR8774458_2.fastq.gz?download=1 ERR8774458_2.fastq.gz

# Make directory to output the results of quality control
mkdir QC_Reports

# Quality control
fastqc ERR8774458_1.fastq.gz ERR8774458_2.fastq.gz -o QC_Reports

#Use MultiQC to summarize the QC results
multiqc .

# Move back to the Project3 directory
cd ..

# Trimming
java -jar ~/bin/trimmomatic/trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 ERR8774458_1.fastq.gz ERR8774458_2.fastq.gz paired1.fastq.gz unpaired1.fastq.gz paired2.fastq.gz unpaired2.fastq.gz SLIDINGWINDOW:3:20 MINLEN:50


# Mapping
# Download the reference genome
wget "https://zenodo.org/records/10886725/files/Reference.fasta?download=1"

# Variant calling 


