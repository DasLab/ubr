# Ubr (pronounced 'uber')

Ubr is a pipeline to process RNA chemical mapping with unique 3' barcodes for each RNA and mutational profiling

U = Ultraplex (for demultiplexing)
B = Bowtie2 (for alignment)
R = RNAFramework (for assigning mutations to chemical modification events)

(c) R. Das, Howard Hughes Medical Institute and Stanford University, 2023

Note: structure modeling & scoring scripts have mostly moved to https://github.com/eternagame/OpenKnotScore. 

## Requirements

You need *Python3*, and the following packages:

- `ultraplex`, available via `pip3 install ultraplex`.
- `bowtie2`, available via [conda](https://anaconda.org/bioconda/bowtie2) or [direct download](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.1/) 
- `RNAFramework`, available from GitHub, with install directions [here](https://rnaframework-docs.readthedocs.io/en/latest/#installation).
- `bbmerge.sh`, available in [bbmap](https://sourceforge.net/projects/bbmap/) -- requires `java`.
- `seqkit`, available via [conda](https://anaconda.org/bioconda/seqkit) or [direct download](https://bioinf.shenwei.me/seqkit/download/)
- `samtools`, available for [download or github](http://www.htslib.org/)
- `MATLAB`, for final data processing and visualizing output. (Ideally these last scripts should be ported to Python/matplotlib to increase accessibility.)

## Example

### Step 1. Process Illumina FASTQ data to get counts (command-line)

The goal of this first of the two steps is to go from raw Illumina reads in FASTQ files to 'counts' files that count up each mutation. 

Following works through an example run that holds information on 2729 sequences from the Eterna OpenKnot Pilot "PK50" library. This RNA library was prepared from DNA pools ordered from both CustomArray/Genscript and Twist. The two RNA libraries were probed with the SHAPE reagent 1M7 and also, as controls for background subtraction, with mock treatment with DMSO.

There are two `FASTA`-formatted with the relvevant information.

 The four experimental conditions are encoded by 12-nt barcodes introduced during reverse transcription, which are encoded in the file `RTBbarcodes_PK50_RNA.fasta`:

```
>RTB008_Twist_PK50_1M7
GGTATATGTACA

>RTB010_CustomArray_PK50_1M7
ATGCATGCACAG

>RTB012_Twist_PK50_nomod
GCAAATGTGCTA

>RTB014_CustomArray_PK50_nomod
CTTTCCCACACT
```

And the 2729 sequences are in `pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa`:

```
>11459092	 I'm knot giving up	 Merida
GGGAACGACUCGAGUAGAGUCGAAAACAUUCCCAAAUUCCACCUUGGUGAUGGCACCCGGAGAGGAGCCAUCACCACACAAAUUUCGAUUUGUGAAAAGAAACAACAACAACAAC

>11459090	 Astros-pnote-11	 Astromon
GGGAACGACUCGAGUAGAGUCGAAAACGGUGCCAGUAUAGCGUAGCGCGGCGGUACACCGCCGAGCUGAGCUGGAAGAAGUAGUUCGCUACUUCAAAAGAAACAACAACAACAAC
...
```

The data are in MiSeq output FASTQ files for Read 1 and Read 2. (Note that in our standard run, we encode experimental conditions not in the i5/i7 regions but within the cDNA itself with the RTB primers.) 

For the example in this repo, 400,000 reads have been extracted from the middle of this file:

```
Sample1_S1_L001_R1_001_400k.fastq.gz
Sample1_S1_L001_R2_001_400k.fastq.gz
```

the full sequencing and example input/output are available at the Google Drive link described below.

There are three ways to run the first step of `UBR`, as (A) single processor run (good for testing), (B) multiple processor run (good for MiSeq-scale data), and (C) cluster run (with SLURM, necessary for large runs, e.g., NovaSeq).

#### Option 1A. Single processor run (good for testing)

These data files are available in `ubr_test`. Run the following command:

```
FASTQ_DIR=../data

ubr_run.py -s pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa  -b RTBbarcodes_PK50_RNA.fasta -1 $FASTQ_DIR/Sample1_S1_L001_R1_001_400k.fastq.gz -2 $FASTQ_DIR/Sample1_S1_L001_R2_001_400k.fastq.gz
```

(This command is also available in the `README_RUN` file, which you can run by typing `source README_RUN`.)

You should get:

```
Read in 2729 sequences from pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa.
Read in 4 primer barcodes from RTBbarcodes_PK50_RNA.fasta.
bbmerge.sh in=../data/Sample1_S1_L001_R1_001_400k.fastq.gz in2=../data/Sample1_S1_L001_R2_001_400k.fastq.gz out=../data/Sample1_S1_L001_R1_001_400k_MERGED.assembled.fastq.gz pigz=f unpigz=f > 0_merge_pairs.out 2> 0_merge_pairs.err

ultraplex -i ../data/Sample1_S1_L001_R1_001_400k_MERGED.assembled.fastq.gz -b primer_barcodes.csv  -d 1_ultraplex/  --dont_build_reference --ignore_no_match --threads 1 > 1_ultraplex/1_ultraplex.out 2> 1_ultraplex/1_ultraplex.err

bowtie2-build seq.fasta 2_bowtie2/bowtie-build/seq --threads 1 > 2_bowtie2/bowtie-build/bt2.out 2> 2_bowtie2/bowtie-build/bt2.err
bowtie2 --end-to-end --sensitive --maxins=800 --ignore-quals --no-unal --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 -x 2_bowtie2/bowtie-build/seq  -U 1_ultraplex/ultraplex_demux_5bc_GGTATATGTACA.fastq.gz -S 2_bowtie2/RTB008_Twist_PK50_1M7/bowtie2.sam --threads 1  > 2_bowtie2/RTB008_Twist_PK50_1M7/bowtie2.out 2> 2_bowtie2/RTB008_Twist_PK50_1M7/bowtie2.err
bowtie2 --end-to-end --sensitive --maxins=800 --ignore-quals --no-unal --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 -x 2_bowtie2/bowtie-build/seq  -U 1_ultraplex/ultraplex_demux_5bc_ATGCATGCACAG.fastq.gz -S 2_bowtie2/RTB010_CustomArray_PK50_1M7/bowtie2.sam --threads 1  > 2_bowtie2/RTB010_CustomArray_PK50_1M7/bowtie2.out 2> 2_bowtie2/RTB010_CustomArray_PK50_1M7/bowtie2.err
bowtie2 --end-to-end --sensitive --maxins=800 --ignore-quals --no-unal --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 -x 2_bowtie2/bowtie-build/seq  -U 1_ultraplex/ultraplex_demux_5bc_GCAAATGTGCTA.fastq.gz -S 2_bowtie2/RTB012_Twist_PK50_nomod/bowtie2.sam --threads 1  > 2_bowtie2/RTB012_Twist_PK50_nomod/bowtie2.out 2> 2_bowtie2/RTB012_Twist_PK50_nomod/bowtie2.err
bowtie2 --end-to-end --sensitive --maxins=800 --ignore-quals --no-unal --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 -x 2_bowtie2/bowtie-build/seq  -U 1_ultraplex/ultraplex_demux_5bc_CTTTCCCACACT.fastq.gz -S 2_bowtie2/RTB014_CustomArray_PK50_nomod/bowtie2.sam --threads 1  > 2_bowtie2/RTB014_CustomArray_PK50_nomod/bowtie2.out 2> 2_bowtie2/RTB014_CustomArray_PK50_nomod/bowtie2.err

rf-count --processors 1 -wt 1 -fast -f seq.fasta -m -cc -rd -ni -ds 1  -orc  -o 3_rf_count/RTB008_Twist_PK50_1M7 2_bowtie2/RTB008_Twist_PK50_1M7//bowtie2.sam >> 3_rf_count/rf-count_RTB008_Twist_PK50_1M7.out 2>> 3_rf_count/rf-count_RTB008_Twist_PK50_1M7.err
gzip 3_rf_count/RTB008_Twist_PK50_1M7/raw_counts/bowtie2.txt
rf-count --processors 1 -wt 1 -fast -f seq.fasta -m -cc -rd -ni -ds 1  -orc  -o 3_rf_count/RTB010_CustomArray_PK50_1M7 2_bowtie2/RTB010_CustomArray_PK50_1M7//bowtie2.sam >> 3_rf_count/rf-count_RTB010_CustomArray_PK50_1M7.out 2>> 3_rf_count/rf-count_RTB010_CustomArray_PK50_1M7.err
gzip 3_rf_count/RTB010_CustomArray_PK50_1M7/raw_counts/bowtie2.txt
rf-count --processors 1 -wt 1 -fast -f seq.fasta -m -cc -rd -ni -ds 1  -orc  -o 3_rf_count/RTB012_Twist_PK50_nomod 2_bowtie2/RTB012_Twist_PK50_nomod//bowtie2.sam >> 3_rf_count/rf-count_RTB012_Twist_PK50_nomod.out 2>> 3_rf_count/rf-count_RTB012_Twist_PK50_nomod.err
gzip 3_rf_count/RTB012_Twist_PK50_nomod/raw_counts/bowtie2.txt
rf-count --processors 1 -wt 1 -fast -f seq.fasta -m -cc -rd -ni -ds 1  -orc  -o 3_rf_count/RTB014_CustomArray_PK50_nomod 2_bowtie2/RTB014_CustomArray_PK50_nomod//bowtie2.sam >> 3_rf_count/rf-count_RTB014_CustomArray_PK50_nomod.out 2>> 3_rf_count/rf-count_RTB014_CustomArray_PK50_nomod.err
gzip 3_rf_count/RTB014_CustomArray_PK50_nomod/raw_counts/bowtie2.txt

rf-rctools view 3_rf_count/RTB008_Twist_PK50_1M7/bowtie2.rc > 4_rctools/RTB008_Twist_PK50_1M7/rf_count.csv && gzip 4_rctools/RTB008_Twist_PK50_1M7/rf_count.csv
rf-rctools view 3_rf_count/RTB010_CustomArray_PK50_1M7/bowtie2.rc > 4_rctools/RTB010_CustomArray_PK50_1M7/rf_count.csv && gzip 4_rctools/RTB010_CustomArray_PK50_1M7/rf_count.csv
rf-rctools view 3_rf_count/RTB012_Twist_PK50_nomod/bowtie2.rc > 4_rctools/RTB012_Twist_PK50_nomod/rf_count.csv && gzip 4_rctools/RTB012_Twist_PK50_nomod/rf_count.csv
rf-rctools view 3_rf_count/RTB014_CustomArray_PK50_nomod/bowtie2.rc > 4_rctools/RTB014_CustomArray_PK50_nomod/rf_count.csv && gzip 4_rctools/RTB014_CustomArray_PK50_nomod/rf_count.csv

Created: RTB008_Twist_PK50_1M7.muts.txt.gz and RTB008_Twist_PK50_1M7.coverage.txt.gz for 2729 sequences with total coverage 44495
Created: RTB010_CustomArray_PK50_1M7.muts.txt.gz and RTB010_CustomArray_PK50_1M7.coverage.txt.gz for 2729 sequences with total coverage 39888
Created: RTB012_Twist_PK50_nomod.muts.txt.gz and RTB012_Twist_PK50_nomod.coverage.txt.gz for 2729 sequences with total coverage 34820
Created: RTB014_CustomArray_PK50_nomod.muts.txt.gz and RTB014_CustomArray_PK50_nomod.coverage.txt.gz for 2729 sequences with total coverage 32928

Created raw_counts/RTB008_Twist_PK50_1M7.AC.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.AG.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.AT.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.CA.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.CG.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.CT.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.GA.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.GC.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.GT.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.TA.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.TC.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.TG.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.ins.txt.gz,raw_counts/RTB008_Twist_PK50_1M7.del.txt.gz for 2729 sequences (found 2319 designs)
Created raw_counts/RTB010_CustomArray_PK50_1M7.AC.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.AG.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.AT.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.CA.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.CG.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.CT.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.GA.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.GC.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.GT.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.TA.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.TC.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.TG.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.ins.txt.gz,raw_counts/RTB010_CustomArray_PK50_1M7.del.txt.gz for 2729 sequences (found 2559 designs)
Created raw_counts/RTB012_Twist_PK50_nomod.AC.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.AG.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.AT.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.CA.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.CG.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.CT.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.GA.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.GC.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.GT.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.TA.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.TC.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.TG.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.ins.txt.gz,raw_counts/RTB012_Twist_PK50_nomod.del.txt.gz for 2729 sequences (found 2275 designs)
Created raw_counts/RTB014_CustomArray_PK50_nomod.AC.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.AG.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.AT.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.CA.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.CG.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.CT.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.GA.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.GC.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.GT.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.TA.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.TC.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.TG.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.ins.txt.gz,raw_counts/RTB014_CustomArray_PK50_nomod.del.txt.gz for 2729 sequences (found 2549 designs)

Timings:
0_merge_pairs 00:00:03
1_ultraplex   00:00:07
2_bowtie2     00:01:43
3_rf_count    00:00:34
4_rctools     00:00:02
Final merge   00:00:01

Total time: 00:02:32
```

The files you need are the `.txt.gz` files that show up in the directory as well as in the `raw_counts/` subdirectory. 

These are simple comma-separated text files which contain one line per sequence (2729 sequences) with integers as follows:

 * `RTB012...`, `RTB014..` tags are primer names, given in the primer file (`RTBbarcodes_PK50_RNA.fasta`).
 * The files represent:
 	* `muts` - number of mutations detected at each position, including deletions (assigned to the right most position of ambiguous homonucleotide stretches).
 	* `coverage` - number of sequencing reads detected at each position
	* `ins`, `del` - number of insertions and deletions detected at each position
	*  `AC`, `AT`, - number of each kind of mutation detected at each position.

Example output is available in the directory `example/EXAMPLE_OUTPUT/ubr_test/`.

These files are what you need for the next step, (see **Step 2** below). 

#### Option 1B. Multiple processors on one node. (Fine for MiSeq)

The above data processing step is efficient for, say, 400,000 lines but becomes slow when processing 15-30M reads for a MiSeq. At that scale, it's more efficient to pre-split the FASTQ into several partitions on which `ubr_run.py` can be run in parallel. 

For preparing that split and the separate run directories, we can run `ubr_split.py`.

An example is given in `example/ubr_test_split`. In that directory, run:

```
FASTQ_DIR=../data

ubr_split.py -q 100000 -s pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa  -b RTBbarcodes_PK50_RNA.fasta -1 $FASTQ_DIR/Sample1_S1_L001_R1_001_400k.fastq.gz -2 $FASTQ_DIR/Sample1_S1_L001_R2_001_400k.fastq.gz
```

(This command-line is also available in `README_SPLIT`, which you can run with `source README_SPLIT`.)

Output should look like:

```
Read in 2729 sequences from pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa.
Read in 4 primer barcodes from RTBbarcodes_PK50_RNA.fasta.

seqkit split2 -s 100000 -1 ../data/Sample1_S1_L001_R1_001_400k.fastq.gz -2 ../data/Sample1_S1_L001_R2_001_400k.fastq.gz -O UBR --threads 12
[INFO] flag -1/--read1 and -2/--read2 given, ignore: -
[INFO] split seqs from ../data/Sample1_S1_L001_R1_001_400k.fastq.gz and ../data/Sample1_S1_L001_R2_001_400k.fastq.gz
[INFO] split into 100000 seqs per file
[INFO] write 100000 sequences to file: UBR/Sample1_S1_L001_R1_001_400k.part_001.fastq.gz
[INFO] write 100000 sequences to file: UBR/Sample1_S1_L001_R2_001_400k.part_001.fastq.gz
[INFO] write 100000 sequences to file: UBR/Sample1_S1_L001_R2_001_400k.part_002.fastq.gz
[INFO] write 100000 sequences to file: UBR/Sample1_S1_L001_R1_001_400k.part_002.fastq.gz
[INFO] write 100000 sequences to file: UBR/Sample1_S1_L001_R1_001_400k.part_003.fastq.gz
[INFO] write 100000 sequences to file: UBR/Sample1_S1_L001_R2_001_400k.part_003.fastq.gz
[INFO] write 100000 sequences to file: UBR/Sample1_S1_L001_R1_001_400k.part_004.fastq.gz
[INFO] write 100000 sequences to file: UBR/Sample1_S1_L001_R2_001_400k.part_004.fastq.gz

Timings:
seqkit split2 00:00:00
file mv/cp    00:00:00

Total time: 00:00:00


All 4 commands can be run with:
 source all_commands.sh

Or you can go into each subdirectory of UBR/ and run:
 source ubr_run.sh

Or to queue up 1 slurm jobs on sherlock you can run:
 source sbatch_commands.sh
```

You can look at the command-lines of the four `ubr_run.py` jobs that need to be queued up by typing `cat all_commands.sh`. You should see:

```
ubr_run.py -s UBR/001/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/001/RTBbarcodes_PK50_RNA.fasta -1 UBR/001/Sample1_S1_L001_R1_001_400k.part_001.fastq.gz -2 UBR/001/Sample1_S1_L001_R2_001_400k.part_001.fastq.gz -O UBR/001 > UBR/001/ubr_run.out 2> UBR/001/ubr_run.err &
ubr_run.py -s UBR/002/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/002/RTBbarcodes_PK50_RNA.fasta -1 UBR/002/Sample1_S1_L001_R1_001_400k.part_002.fastq.gz -2 UBR/002/Sample1_S1_L001_R2_001_400k.part_002.fastq.gz -O UBR/002 > UBR/002/ubr_run.out 2> UBR/002/ubr_run.err &
ubr_run.py -s UBR/003/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/003/RTBbarcodes_PK50_RNA.fasta -1 UBR/003/Sample1_S1_L001_R1_001_400k.part_003.fastq.gz -2 UBR/003/Sample1_S1_L001_R2_001_400k.part_003.fastq.gz -O UBR/003 > UBR/003/ubr_run.out 2> UBR/003/ubr_run.err &
ubr_run.py -s UBR/004/pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa -b UBR/004/RTBbarcodes_PK50_RNA.fasta -1 UBR/004/Sample1_S1_L001_R1_001_400k.part_004.fastq.gz -2 UBR/004/Sample1_S1_L001_R2_001_400k.part_004.fastq.gz -O UBR/004 > UBR/004/ubr_run.out 2> UBR/004/ubr_run.err &
```

If you're runnig on a computer that can run at least 4 threads, then go ahead and type:

```
source all_commands.sh 
```

You'll see four processes queue up, which you can monitor with, e.g., `top`. You'll see files show up in the subdirectories `UBR/001`,...`UBR/004`.

Or you can look at the log files:

```
tail -f UBR/001/ubr_run.out
```
When jbs are done, the end of the file will report timings, e.g.:

```
Timings:
0_merge_pairs 00:00:03
1_ultraplex   00:00:02
2_bowtie2     00:00:28
3_rf_count    00:00:15
4_rctools     00:00:02
Final merge   00:00:01

Total time: 00:00:53
```

To collate counts from the four runs on the four splits, now type:

```
ubr_merge.py UBR
```

Output will be:

```
['UBR/*/']

Will merge:
RTB008_Twist_PK50_1M7.coverage.txt.gz
RTB010_CustomArray_PK50_1M7.coverage.txt.gz
RTB012_Twist_PK50_nomod.coverage.txt.gz
...
RTB014_CustomArray_PK50_nomod.ins.txt.gz
RTB008_Twist_PK50_1M7.del.txt.gz
RTB010_CustomArray_PK50_1M7.del.txt.gz
RTB012_Twist_PK50_nomod.del.txt.gz
RTB014_CustomArray_PK50_nomod.del.txt.gz
Compiled    44495 total counts for   2729 sequences from      4 of      4 files into: RTB008_Twist_PK50_1M7.coverage.txt.gz
Compiled    39888 total counts for   2729 sequences from      4 of      4 files into: RTB010_CustomArray_PK50_1M7.coverage.txt.gz
Compiled    34820 total counts for   2729 sequences from      4 of      4 files into: RTB012_Twist_PK50_nomod.coverage.txt.gz
Compiled    32928 total counts for   2729 sequences from      4 of      4 files into: RTB014_CustomArray_PK50_nomod.coverage.txt.gz
Compiled     4514 total counts for   2729 sequences from      4 of      4 files into: RTB008_Twist_PK50_1M7.muts.txt.gz
Compiled     5044 total counts for   2729 sequences from      4 of      4 files into: RTB010_CustomArray_PK50_1M7.muts.txt.gz
Compiled     2634 total counts for   2729 sequences from      4 of      4 files into: RTB012_Twist_PK50_nomod.muts.txt.gz
...
Compiled     2889 total counts for   2729 sequences from      4 of      4 files into: raw_counts/RTB012_Twist_PK50_nomod.ins.txt.gz
Compiled     3309 total counts for   2729 sequences from      4 of      4 files into: raw_counts/RTB014_CustomArray_PK50_nomod.ins.txt.gz
Compiled     3200 total counts for   2729 sequences from      4 of      4 files into: raw_counts/RTB008_Twist_PK50_1M7.del.txt.gz
Compiled     3792 total counts for   2729 sequences from      4 of      4 files into: raw_counts/RTB010_CustomArray_PK50_1M7.del.txt.gz
Compiled      885 total counts for   2729 sequences from      4 of      4 files into: raw_counts/RTB012_Twist_PK50_nomod.del.txt.gz
Compiled     2356 total counts for   2729 sequences from      4 of      4 files into: raw_counts/RTB014_CustomArray_PK50_nomod.del.txt.gz

Timings:
Read in data: 00:00:02
Output  data: 00:00:02

Total time: 00:00:05
```

Again, note that the above example is for 400,000 reads, split into 4 partitions of 100,000 reads. 

In actuality, a MiSeq run will produce 15-30M reads, which you'll want to partition into, say,  8 splits of 2M-4M reads each, by providing an option like `-q 2000000` to `ubr_split.py` above. It would then take roughly 10-30 mins to run each split on one thread on a typical desktop.  

Example output is available in the directory `example/EXAMPLE_OUTPUT/ubr_test_split/`.

#### Option 1C. Cluster runs (SLURM; needed for NovaSeq scale runs)

For NovaSeq runs (1B-8B reads), you'll likely need to do a cluster run. 

The data for a full example directory, even for a Miseq, is too big to store in a GitHub repo, but the example data are available in this [Google Drive folder](). 

Download those data, and put the files in the `example/data_full/` directory.

Then, go into `example/ubr_full`.  The `ubr_split.py` command is provided in `run_ubr_split.sh`, which contains:

```
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
```

This is set up with appropriate parameters for SLURM to run on Stanford's Sherlock cluster, asking for a node from the `biochem,owners` partition -- you may need to edit that line to point to the partitions on which you have access. 

You can queue up the `ubr_split` command if you're on a SLURM cluster, by typing:

```
sbatch run_ubr_split.sh
```

That can take a while. For this example, it actually doesn't take long to run the split (only 5 minutes), so you can also just directly run the command with `sh run_ubr_split.sh`. However, for big jobs (e.g., NovaSeq), the split can take 24 hrs. It's best to not hog a login node to do the split, or to risk a timeout, so in those cases, use `sbatch` to queue up the job.

The splits are 1M reads per job and in this case lead to 18 partitions, since the starting Miseq run has XXX reads.

After the `ubr_split.py` completes (look for `DONE` to show up in the out file!), you can run the jobs on slurm with:

```
source sbatch_commands.sh
```

This runs one or more batches of jobs in which 24 `ubr_run.py` jobs will be run on each node.  Monitor as above, looking inside `ubr_run.out` files in the `UBR/*/` subdirectories.

When done, you can run `ubr_merge.py UBR/`. That can take a super long time if you have a lot of subdirectories, in which case type:

```
ubr_merge.py UBR -s
```

This will trigger the definition of a bunch of parallel jobs which will collate each of the `*.txt.gz` counts files in separate threads. Example output below:

```
```

Go ahead and queue up the `ubr_merge.py` jobs with `sbatch ubr_merge_commands.sh`.

In the end, you should have the `.txt.gz` files in the `ubr_full` example directory and the `raw_counts/` subdirectory.

The core data that need to be saved -- ideally transferred from the cluster to a laptop/desktop, which should be synced to Drive/Box/DropBox -- are the files in the directory without the massive `UBR` subdirectory. That's where you can run step 2.

Example output (*not* including the UBR subdirectories) is in `example/EXAMPLE_OUTPUT/ubr_full/`. Full example output is available at the [Google Drive link]().


### Step 2. Process counts to get & visualize reactivities (MATLAB)



