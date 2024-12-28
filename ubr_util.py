#!/usr/bin/env python3
import argparse
import os

parser = argparse.ArgumentParser(
                    prog = 'ubr_utile.py',
                    description = 'Utility scripts for UBR' )

parser.add_argument('-s','--sequences_fasta', required=True, help='FASTA of RNA sequences')
parser.add_argument('-t','--threads',default=1, type=int, help='Number of threads for Bowtie2 and RNAFramework')
parser.add_argument('--build_bowtie2_index', action='store_true', help='Prepare bowtie-build directory from FASTA file')

def check_sequence(sequence):
    # check sequences are RNA/DNA
    for c in sequence:
        if c not in 'ACGTUN': return False
    return True

def check_dup(mylist,mytag, force = False):
    if (len(set(mylist)) != len(mylist)):
        print('%s not unique! %d are duplicates.' % (mytag,len(mylist)-len(set(mylist))) )
        myset = set()
        for (i,x) in enumerate(mylist):
            if x in myset: print('%6d %s' % (i+1,x) )
            myset.add(x)
        if not force: exit(0)
    return

def read_fasta( fasta_file, force = False ):
    lines = open( fasta_file ).readlines()
    sequences = []
    headers = []
    header = None
    sequence = ''
    for line in lines:
        if len(line) > 0 and line[0] == '>':
            if header is not None:
                headers.append(header)
                sequences.append(sequence)
            sequence = ''
            header = line[1:].strip('\n')
            continue
        sequence = sequence + line.strip('\n')
    if header is not None:
        headers.append(header)
        sequences.append(sequence)
    assert( len(sequences) == len(headers ) )
    check_dup( sequences,'Sequences', force )
    check_dup( headers,'Headers', force )
    return (sequences,headers)

def create_seq_fasta( sequences, headers, wd = './' ):
    # Need to make sure we have DNA versions
    seq_file = wd + 'seq.fasta'
    fid = open( seq_file, 'w' )
    for (sequence,header) in zip(sequences, headers): fid.write('>%s\n%s\n' % (header,sequence.replace('U','T')) )
    fid.close()
    return seq_file

def build_bowtie2_index( sequences_fasta ):
    (sequences,headers) = read_fasta( sequences_fasta )
    seq_file = create_seq_fasta( sequences, headers )

    bowtie_build_dir = 'bowtie-build'
    bt2_prefix = '%s/seq'% (bowtie_build_dir)
    bt2_file = '%s.1.bt2' % bt2_prefix

    os.makedirs(bowtie_build_dir,exist_ok = True)
    command = 'bowtie2-build %s %s --threads %d > %s/bt2.out 2> %s/bt2.err' % (seq_file,bt2_prefix,args.threads,bowtie_build_dir,bowtie_build_dir)
    print(command)
    os.system( command )
    print('Created: %s' % bowtie_build_dir )

if __name__ == "__main__":
    args = parser.parse_args()

    if args.build_bowtie2_index:
        build_bowtie2_index( args.sequences_fasta )
