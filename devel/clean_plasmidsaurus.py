#!/usr/bin/env python3

import argparse
import os
import subprocess
import glob

def run_cutadapt(input_file, output_file, adapt, preadapt):
    cmd = [
        "cutadapt",
        "-a", adapt,
        "-o", output_file,
        input_file,
        "--trimmed-only",
        "-g", preadapt,
        "--revcomp"
    ]
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser(description="Clean plasmidsaurus FASTQ files using cutadapt.")
    parser.add_argument("input_files", nargs="+", help="Input FASTQ/FASTQ.gz file(s)")
    parser.add_argument("-a", "--adapt", default="AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG", help="Adapter sequence")
    parser.add_argument("-g", "--preadapt", default="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", help="Pre-adapter sequence")
    parser.add_argument("-o", "--output", default="out.fastq.gz", help="Output FASTQ file")

    args = parser.parse_args()

    temp_files = []

    for input_file in args.input_files:
        base_name = os.path.basename(input_file)
        temp_file1 = f"{base_name}.0.gz"
        temp_file2 = f"{base_name}.1.gz"

        run_cutadapt(input_file, temp_file1, args.adapt, args.preadapt)
        run_cutadapt(temp_file1, temp_file2, args.adapt, args.preadapt)

        temp_files.extend([temp_file1, temp_file2])

    # Concatenate all .1.gz files
    with open(args.output, 'wb') as outfile:
        for filename in glob.glob('*.1.gz'):
            with open(filename, 'rb') as readfile:
                outfile.write(readfile.read())

    # Remove intermediate files
    for temp_file in temp_files:
        os.remove(temp_file)

if __name__ == "__main__":
    main()
