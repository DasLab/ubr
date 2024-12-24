#!/usr/bin/env python3
import argparse
import os
import time
import shutil
import gzip
import glob
import random

parser = argparse.ArgumentParser(
                    prog = 'prepare_fasta_from_bottleneck_sam.py',
                    description = 'looks for specific sequences that are seen repeatedly in SAM file',
                    epilog = 'Make sure to run samtools sort and samtools markdup --duplicate_count on SAM file')

parser.add_argument('sam_file', help='SAM file from DNA sequencing run of bottlenecked library')
parser.add_argument('-s','--sequence',default='GGGAACGACTCGAGTAGAGTCGAAAAACACTGNNNNNNNCAGTGAGCATGNNNNCTACGTTCGCGTAGNNNNCATGCAAAAGAAACAACAACAACAAC',help='sequence to match')
parser.add_argument('-o','--output_fasta',default='unique_sequences.fasta',help='Name of output file')

args = parser.parse_args()
ref_seq = args.sequence
rna_length = len(ref_seq)
sam_file = args.sam_file

assert( shutil.which( 'samtools' ) )

n_pos = []
match_pos = []
for i in range(len(ref_seq)):
    if ref_seq[i] == 'N': n_pos.append(i)
    else: match_pos.append(i)
assert( len(n_pos) > 0 )

lines = os.popen( 'samtools view %s' % sam_file ).readlines()
headers = []
rna_seqs = []
all_counts = 0
out_counts = 0
for line in lines:
    cols = line.strip().split()
    dc = cols[-1]
    assert( dc[:5]=='dc:i:' )
    count = dc[5:]
    cigar = cols[5]
    seq = cols[9]
    all_counts += int( count )

    # Check for exact sequence match
    if cigar != '%dM' % rna_length: continue # later can rescue indels
    matches = True
    if seq.find('N')>0: continue
    for i in match_pos:
        if seq[i] != ref_seq[i]:
            matches= False
            break
    if not matches: continue

    out_counts += int( count )
    rna_seq = seq.replace('T','U')
    header = '>'
    for (i,pos) in enumerate(n_pos):
        if (i>0 and pos > n_pos[i-1]+1): header += '_'
        header += rna_seq[pos]

    headers.append(header)
    rna_seqs.append(rna_seq)
    #print( header,seq )



fid = open(args.output_fasta,'w')
for (header,seq) in zip(headers,rna_seqs):
      fid.write('%s\n%s\n\n' % (header,rna_seq))
fid.close()
print('Found %d of %d counts in sequences to output.' % (out_counts,all_counts) )
print('Outputted %d of %d sequences to %s' % (len(rna_seqs),len(lines), args.output_fasta))





