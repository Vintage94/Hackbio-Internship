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

# Trimming and Filtering
java -jar ~/bin/trimmomatic/trimmomatic-0.39/trimmomatic-0.39.jar PE -@ 8 "$name"_R1.fastq.gz "$name"_R2.fastq.gz \
    "$name"_R1_trimmed.fastq.gz "$name"_R1_unpaired.fastq.gz \
    "$name"_R2_trimmed.fastq.gz "$name"_R2_unpaired.fastq.gz \
    TRAILING:20 MINLEN:50

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

# Convert the vcf file to a .csv file for further analysis and visualisation
echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER" > output.csv
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' variants.vcf.gz >> output.csv



