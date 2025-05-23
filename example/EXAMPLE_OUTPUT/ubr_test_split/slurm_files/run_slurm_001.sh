#!/bin/bash
#SBATCH --job-name=ubr_run
#SBATCH --output=ubr_run.o%j
#SBATCH --error=ubr_run.e%j
#SBATCH --partition=biochem,owners
#SBATCH --time=48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=48G

ubr_run.py -s UBR/001/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/001/RTBbarcodes_PK50_RNA.fasta -1 UBR/001/Sample1_S1_L001_R1_001_400k.part_001.fastq.gz -2 UBR/001/Sample1_S1_L001_R2_001_400k.part_001.fastq.gz -O UBR/001 > UBR/001/ubr_run.out 2> UBR/001/ubr_run.err &
ubr_run.py -s UBR/002/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/002/RTBbarcodes_PK50_RNA.fasta -1 UBR/002/Sample1_S1_L001_R1_001_400k.part_002.fastq.gz -2 UBR/002/Sample1_S1_L001_R2_001_400k.part_002.fastq.gz -O UBR/002 > UBR/002/ubr_run.out 2> UBR/002/ubr_run.err &
ubr_run.py -s UBR/003/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/003/RTBbarcodes_PK50_RNA.fasta -1 UBR/003/Sample1_S1_L001_R1_001_400k.part_003.fastq.gz -2 UBR/003/Sample1_S1_L001_R2_001_400k.part_003.fastq.gz -O UBR/003 > UBR/003/ubr_run.out 2> UBR/003/ubr_run.err &
ubr_run.py -s UBR/004/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/004/RTBbarcodes_PK50_RNA.fasta -1 UBR/004/Sample1_S1_L001_R1_001_400k.part_004.fastq.gz -2 UBR/004/Sample1_S1_L001_R2_001_400k.part_004.fastq.gz -O UBR/004 > UBR/004/ubr_run.out 2> UBR/004/ubr_run.err &

wait
echo "DONE"
