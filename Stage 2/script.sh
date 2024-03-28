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

# Use MultiQC to summarize the QC results
multiqc .

# Move back to the Project3 directory
cd ..

# Trimming and Filtering
java -jar ~/bin/trimmomatic/trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 ERR8774458_1.fastq.gz ERR8774458_2.fastq.gz paired1.fastq.gz unpaired1.fastq.gz paired2.fastq.gz unpaired2.fastq.gz SLIDINGWINDOW:3:20 MINLEN:50

# Decompress the files
gunzip paired1.fastq.gz
gunzip paired2.fastq.gz

# Download the reference genome
wget "https://zenodo.org/records/10886725/files/Reference.fasta?download=1"

# Rename the reference genome
mv Reference.fasta?download=1 Reference.fasta

# Make a new directory 'Mapping' to output our mapping results
mkdir Mapping

# Index the reference genome
bwa index Reference.fasta

# Map the preprocessed sequences to the reference genome
bwa mem Reference.fasta paired1.fastq paired2.fastq > Mapping/aligned.sam

# Change back to the Project3 directory
cd..

# Convert the sam file to a bam file
samtools view -@ 20 -S -b Mapping/aligned.sam > Mapping/aligned.bam

# Sort bam file based on the order the reads were mapped
samtools sort -@ 32 -o Mapping/aligned.sorted.bam Mapping/aligned.bam

# Index the sorted bam file
samtools index Mapping/aligned.sorted.bam

# Variant calling 
# Index the reference genome using samtools
samtools faidx Reference.fasta

# Variant calling with freebayes
freebayes -f Reference.fasta Mapping/aligned.sorted.bam > variants.vcf

# Compress the vcf file
bgzip variants.vcf

# Decopress the vcf file
tabix variants.vcf.gz


