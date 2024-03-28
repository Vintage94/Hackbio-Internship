#!/usr/bin/env bash

# Download the forward sequence
wget "https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz?download=1"

# Download the reverse sequence
wget "https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz?download=1"

# Rename the files
mv ERR8774458_1.fastq.gz?download=1 ERR8774458_1.fastq.gz
mv ERR8774458_2.fastq.gz?download=1 ERR8774458_2.fastq.gz

# Quality control


# Trimming


# Mapping


# Variant calling 


