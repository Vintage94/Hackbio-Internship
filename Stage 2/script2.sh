#!/usr/bin/env bash
# Error handling: Exit script if any command fails or encounters unset variables
set -eu

# Define the URL for the reference genome and the desired filename
reference_url="https://example.com/reference_genome.fasta"
reference_name="reference_genome.fasta"

# Download the reference genome if it's not already present
if [ ! -f "$reference_name" ]; then
    echo "Downloading reference genome..."
    wget -O "$reference_name" "$reference_url"
    echo "Downloaded reference genome."
fi

# Define an array of sample URLs
sample_urls=(
  "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz"
  "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz"
  "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz"
  "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz"
  "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz"
  "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz"
  "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz"
  "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz"
  "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz"
  "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz"
)

# Index the reference genome for mapping (if not already indexed)
if [ ! -f "$reference_name".bwt ]; then
    bwa index "$reference_name"
    echo "Indexed reference genome for mapping."
fi

# Loop through each sample URL
for url in "${sample_urls[@]}"; do
  # Extract the sample name from the URL
  name=$(basename "$url" | cut -d '_' -f 1)
  echo "Processing sample: $name"
  
  # Download the dataset
  wget -O "$name.fastq.gz" "$url"
  echo "Downloaded dataset for $name."

  # Create directories for quality control and mapping
  mkdir -p "$name"/QC_Reports
  mkdir -p "$name"/Mapping
  echo "Created directories for quality control and mapping for $name."
  
  # Quality control
  fastqc "$name.fastq.gz" -o "$name"/QC_Reports
  echo "Performed quality control for $name."

  # Summarize QC results
  multiqc "$name"/QC_Reports
  echo "Summarized QC results for $name."
  
  # Trim using Trimomatic
  java -jar trimmomatic-0.39.jar PE -threads 8 "$name"_R1.fastq.gz "$name"_R2.fastq.gz \
    "$name"_R1_trimmed.fastq.gz "$name"_R1_unpaired.fastq.gz \
    "$name"_R2_trimmed.fastq.gz "$name"_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  echo "Trimmed reads for $name."

  # Mapping
  bwa mem "$reference_name" "$name"_R1_trimmed.fastq.gz "$name"_R2_trimmed.fastq.gz > "$name"/Mapping/"$name".sam
  samtools view -@ 20 -S -b "$name"/Mapping/"$name".sam > "$name"/Mapping/"$name".bam
  samtools sort -@ 32 -o "$name"/Mapping/"$name".sorted.bam "$name"/Mapping/"$name".bam
  samtools index "$name"/Mapping/"$name".sorted.bam
  echo "Performed mapping for $name."

  # Index the reference genome for variant calling (if not already indexed)
  if [ ! -f "$reference_name".fai ]; then
      samtools faidx "$reference_name"
      echo "Indexed reference genome for variant calling."
  fi

  # Variant calling using freebayes
  freebayes -f "$reference_name" "$name"/Mapping/"$name".sorted.bam > "$name"/"$name".vcf
  echo "Performed variant calling for $name."

  # Compress and index the VCF file
  bgzip "$name"/"$name".vcf
  tabix "$name"/"$name".vcf.gz
  echo "Compressed and indexed VCF file for $name."

  # Convert VCF to CSV
  echo -e "Sample\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER" > "$name"/"$name".csv
  bcftools query -l "$name"/"$name".vcf.gz | sed 's/$/'$'\t''/' > "$name"/sample_column.txt
  paste "$name"/sample_column.txt <(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' "$name"/"$name".vcf.gz) >> "$name"/"$name".csv
  echo "Converted VCF to CSV for $name."

  # Move the CSV file to the main directory
  mv "$name"/"$name".csv .

done

# Merge CSV files for all samples
echo "Merging CSV files..."
cat *.csv > merged.csv
echo "Merged CSV files for all samples."

echo "Analysis completed successfully!"
