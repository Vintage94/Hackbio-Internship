#!/usr/bin/env bash

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

# Trimming


# Mapping


# Variant calling 


