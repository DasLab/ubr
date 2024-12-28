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
parser.add_argument('-o','--output_fasta',default='',help='Name of output file')
parser.add_argument('--markdup',action='store_true',help='Run samtools sort and markdup')
parser.add_argument('--min_count',type=int,default=1,help='Minimum number of counts to output to fasta file')
parser.add_argument('-v','--verbose',action='store_true',help='Verbose output')

args = parser.parse_args()
ref_seq = args.sequence
rna_length = len(ref_seq)
sam_file = args.sam_file

assert( shutil.which( 'samtools' ) )

if args.markdup:
    assert( sam_file.find('.sam') > -1 )
    sort_file = sam_file.replace('.sam','.sorted.sam')
    if not os.path.isfile( sort_file ):
        command = 'samtools sort %s -o %s' % (sam_file, sort_file)
        print( command )
        os.system( command )
    assert( os.path.isfile( sort_file ) )
    markdup_file = sam_file.replace('.sam','.sorted.markdup.sam')
    if not os.path.isfile( markdup_file ):
        command = 'samtools markdup -r -s -f %s %s %s --duplicate-count' % (markdup_file+'.txt',sort_file,markdup_file)
        print( command )
        os.system( command )
    assert( os.path.isfile( markdup_file ) )
    sam_file = markdup_file

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
pass_length = 0
pass_min_count = 0
pass_no_n = 0
out_counts = 0
ref_printed = False
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
    pass_length += 1

    if int(count) < args.min_count: continue
    pass_min_count += 1

    matches = True
    if seq.find('N')>0: continue
    pass_no_n += 1

    mismatches = []
    mismatch_pos = []
    for i in match_pos:
        if seq[i] != ref_seq[i]:
            matches= False
            mismatch_pos.append(i)

    if len(mismatch_pos) > 0:
        if args.verbose:
            if not ref_printed:
                print(ref_seq, 'ref')
                ref_printed = True
            for q in range(len(ref_seq)):
                if q in mismatch_pos: print('X',end='')
                else: print(' ',end='')
            print()
            print(seq, cols[0])
        continue

    out_counts += int( count )
    rna_seq = seq.replace('T','U')
    header = '>'
    for (i,pos) in enumerate(n_pos):
        if (i>0 and pos > n_pos[i-1]+1): header += '_'
        header += rna_seq[pos]

    headers.append(header)
    rna_seqs.append(rna_seq)
    #print( header,seq )


output_fasta = args.output_fasta
if len(output_fasta) == 0: output_fasta = sam_file+'.fasta'

fid = open(output_fasta,'w')
for (header,rna_seq) in zip(headers,rna_seqs):
      fid.write('%s\n%s\n\n' % (header,rna_seq))
fid.close()
if args.verbose:
    print('all: %d, pass_length: %d, pass_min_count: %d, pass_no_n: %d, pass_match: %d' % (len(lines),pass_length,pass_min_count,pass_no_n,len(rna_seqs)) )
    print('Found %8d of %8d counts in sequences to output.' % (out_counts,all_counts) )
print('Outputted %8d of %8d sequences to %s' % (len(rna_seqs),len(lines),output_fasta))





