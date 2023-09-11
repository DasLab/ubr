#!/bin/bash
#SBATCH --job-name=ubr_split
#SBATCH --output=ubr_split.out
#SBATCH --error=ubr_split.err
#SBATCH --partition=biochem,owners
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G

FASTQ_DIR=../data_full

ubr_split.py -q 1000000 -s pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa  -b RTBbarcodes_PK50_RNA.fasta -1 $FASTQ_DIR/Sample1_S1_L001_R1_001.fastq.gz -2 $FASTQ_DIR/Sample1_S1_L001_R2_001.fastq.gz

echo "DONE"
