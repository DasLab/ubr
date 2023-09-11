#!/bin/bash
#SBATCH --job-name=ubr_run
#SBATCH --output=ubr_run.o%j
#SBATCH --error=ubr_run.e%j
#SBATCH --partition=biochem,owners
#SBATCH --time=48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=48G

ubr_run.py -s UBR/001/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/001/RTBbarcodes_PK50_RNA.fasta -1 UBR/001/Sample1_S1_L001_R1_001.part_001.fastq.gz -2 UBR/001/Sample1_S1_L001_R2_001.part_001.fastq.gz -O UBR/001 > UBR/001/ubr_run.out 2> UBR/001/ubr_run.err &
ubr_run.py -s UBR/002/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/002/RTBbarcodes_PK50_RNA.fasta -1 UBR/002/Sample1_S1_L001_R1_001.part_002.fastq.gz -2 UBR/002/Sample1_S1_L001_R2_001.part_002.fastq.gz -O UBR/002 > UBR/002/ubr_run.out 2> UBR/002/ubr_run.err &
ubr_run.py -s UBR/003/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/003/RTBbarcodes_PK50_RNA.fasta -1 UBR/003/Sample1_S1_L001_R1_001.part_003.fastq.gz -2 UBR/003/Sample1_S1_L001_R2_001.part_003.fastq.gz -O UBR/003 > UBR/003/ubr_run.out 2> UBR/003/ubr_run.err &
ubr_run.py -s UBR/004/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/004/RTBbarcodes_PK50_RNA.fasta -1 UBR/004/Sample1_S1_L001_R1_001.part_004.fastq.gz -2 UBR/004/Sample1_S1_L001_R2_001.part_004.fastq.gz -O UBR/004 > UBR/004/ubr_run.out 2> UBR/004/ubr_run.err &
ubr_run.py -s UBR/005/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/005/RTBbarcodes_PK50_RNA.fasta -1 UBR/005/Sample1_S1_L001_R1_001.part_005.fastq.gz -2 UBR/005/Sample1_S1_L001_R2_001.part_005.fastq.gz -O UBR/005 > UBR/005/ubr_run.out 2> UBR/005/ubr_run.err &
ubr_run.py -s UBR/006/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/006/RTBbarcodes_PK50_RNA.fasta -1 UBR/006/Sample1_S1_L001_R1_001.part_006.fastq.gz -2 UBR/006/Sample1_S1_L001_R2_001.part_006.fastq.gz -O UBR/006 > UBR/006/ubr_run.out 2> UBR/006/ubr_run.err &
ubr_run.py -s UBR/007/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/007/RTBbarcodes_PK50_RNA.fasta -1 UBR/007/Sample1_S1_L001_R1_001.part_007.fastq.gz -2 UBR/007/Sample1_S1_L001_R2_001.part_007.fastq.gz -O UBR/007 > UBR/007/ubr_run.out 2> UBR/007/ubr_run.err &
ubr_run.py -s UBR/008/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/008/RTBbarcodes_PK50_RNA.fasta -1 UBR/008/Sample1_S1_L001_R1_001.part_008.fastq.gz -2 UBR/008/Sample1_S1_L001_R2_001.part_008.fastq.gz -O UBR/008 > UBR/008/ubr_run.out 2> UBR/008/ubr_run.err &
ubr_run.py -s UBR/009/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/009/RTBbarcodes_PK50_RNA.fasta -1 UBR/009/Sample1_S1_L001_R1_001.part_009.fastq.gz -2 UBR/009/Sample1_S1_L001_R2_001.part_009.fastq.gz -O UBR/009 > UBR/009/ubr_run.out 2> UBR/009/ubr_run.err &
ubr_run.py -s UBR/010/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/010/RTBbarcodes_PK50_RNA.fasta -1 UBR/010/Sample1_S1_L001_R1_001.part_010.fastq.gz -2 UBR/010/Sample1_S1_L001_R2_001.part_010.fastq.gz -O UBR/010 > UBR/010/ubr_run.out 2> UBR/010/ubr_run.err &
ubr_run.py -s UBR/011/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/011/RTBbarcodes_PK50_RNA.fasta -1 UBR/011/Sample1_S1_L001_R1_001.part_011.fastq.gz -2 UBR/011/Sample1_S1_L001_R2_001.part_011.fastq.gz -O UBR/011 > UBR/011/ubr_run.out 2> UBR/011/ubr_run.err &
ubr_run.py -s UBR/012/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/012/RTBbarcodes_PK50_RNA.fasta -1 UBR/012/Sample1_S1_L001_R1_001.part_012.fastq.gz -2 UBR/012/Sample1_S1_L001_R2_001.part_012.fastq.gz -O UBR/012 > UBR/012/ubr_run.out 2> UBR/012/ubr_run.err &
ubr_run.py -s UBR/013/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/013/RTBbarcodes_PK50_RNA.fasta -1 UBR/013/Sample1_S1_L001_R1_001.part_013.fastq.gz -2 UBR/013/Sample1_S1_L001_R2_001.part_013.fastq.gz -O UBR/013 > UBR/013/ubr_run.out 2> UBR/013/ubr_run.err &
ubr_run.py -s UBR/014/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/014/RTBbarcodes_PK50_RNA.fasta -1 UBR/014/Sample1_S1_L001_R1_001.part_014.fastq.gz -2 UBR/014/Sample1_S1_L001_R2_001.part_014.fastq.gz -O UBR/014 > UBR/014/ubr_run.out 2> UBR/014/ubr_run.err &
ubr_run.py -s UBR/015/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/015/RTBbarcodes_PK50_RNA.fasta -1 UBR/015/Sample1_S1_L001_R1_001.part_015.fastq.gz -2 UBR/015/Sample1_S1_L001_R2_001.part_015.fastq.gz -O UBR/015 > UBR/015/ubr_run.out 2> UBR/015/ubr_run.err &
ubr_run.py -s UBR/016/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/016/RTBbarcodes_PK50_RNA.fasta -1 UBR/016/Sample1_S1_L001_R1_001.part_016.fastq.gz -2 UBR/016/Sample1_S1_L001_R2_001.part_016.fastq.gz -O UBR/016 > UBR/016/ubr_run.out 2> UBR/016/ubr_run.err &
ubr_run.py -s UBR/017/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/017/RTBbarcodes_PK50_RNA.fasta -1 UBR/017/Sample1_S1_L001_R1_001.part_017.fastq.gz -2 UBR/017/Sample1_S1_L001_R2_001.part_017.fastq.gz -O UBR/017 > UBR/017/ubr_run.out 2> UBR/017/ubr_run.err &
ubr_run.py -s UBR/018/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/018/RTBbarcodes_PK50_RNA.fasta -1 UBR/018/Sample1_S1_L001_R1_001.part_018.fastq.gz -2 UBR/018/Sample1_S1_L001_R2_001.part_018.fastq.gz -O UBR/018 > UBR/018/ubr_run.out 2> UBR/018/ubr_run.err &

wait
echo "DONE"
