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

## Example

*TODO: fill out*

