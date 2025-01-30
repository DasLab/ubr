#! /usr/bin/env python3
import argparse
from os import popen

parser = argparse.ArgumentParser(
                    prog = 'sam_to_2d.py',
                    description = 'Processes SAM file to create 1D and 2D mut/del profiles',
                    epilog = 'Extremely basic processing of CIGAR string in bowtie2 SAM file to produce alignments and then coding.')

parser.add_argument('--sam', nargs='*',default='RTB010_CustomArray_PK50_1M7_BindTheFivePrimeEnd.sam.txt',help='SAM-formatted text file, Bowtie2')
parser.add_argument('-s','--sequences_fasta', required=True, default='pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa',help='Fasta file with DNA/RNA sequences used for Bowtie2')
parser.add_argument('--mut_cutoff',type=int,default=10,help='Filter for maximum number of mut/del in read' )
parser.add_argument('-o','--out_file_tag',type=str,default='',help='Tag for outfile [default: derive from first samfile]' )

args = parser.parse_args()

def read_fasta( fasta_file ):
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
    return (sequences,headers)

# figure out reference sequence
def figure_out_sequence(id):
    found_sequence = 0
    for (sequence,header) in zip(sequences,headers):
        if header.find(id) > -1:
            found_sequence = 1
            break
    assert(found_sequence)
    full_sequence = sequence.replace('U','T')
    #print(header,full_sequence,len(full_sequence))
    return (full_sequence,id)

(sequences,headers) = read_fasta( args.sequences_fasta )

out_file_tag = args.out_file_tag
if len( out_file_tag ) == 0:
    if len( args.sam ) == 1:
        out_file_tag = args.sam[0].replace('.sam.txt','').replace('.sam','').replace('.bam','')
    else:
        print('Supply out_file_tag with -o since you have multiple sam files.')
        exit(0)

counts_2d = {}
count = 0
count_filter = 0
for samfile in args.sam:
    print( 'Doing SAM file %s' % samfile )
    if  samfile.find('.sam.txt')>-1:
        fid = open(samfile)
    else:
        assert( samfile[-4:]=='.sam' or  samfile[-4:]=='.bam' )
        fid = popen( 'samtools view %s' % samfile )

    line = fid.readline()
    while line[0]=='@': line = fid.readline()

    cols = line.split()
    (full_sequence,id) = figure_out_sequence(cols[2])
    Nres = len(full_sequence)

    if len( counts_2d ) == 0:
        for i in ['mut','del']:
            counts_2d[i] = {}
            for j in ['mut','del']:
                counts_2d[i][j] = []
                for n in range(Nres):
                    counts_2d[i][j].append( [0]*Nres )

    while line:
        cols = line.split()
        start_pos = int(cols[3])
        cigar = cols[5]
        signed_tmpl_len = int(cols[8])
        seq = cols[9]
        if ( id != cols[2] ): (full_sequence,id) = figure_out_sequence(cols[2])

        ###################################
        # create alignment output for seq.
        ###################################
        # parse cigar string, e.g., "22M1I23M1D69M" => "22M 1I 23M 1D 69M"
        cpos = 0 # cigar position
        spos = 0 # sequence position
        seqa = ''
        seqa += '.'*(start_pos-1)
        for k,s in enumerate(cigar):
            num = cigar[cpos:(k+1)]
            if not num.isnumeric():
                nres = int(cigar[cpos:k])
                indelcode = cigar[k]
                if indelcode == 'M':
                    seqa += seq[spos:(spos+nres)]
                    spos += nres
                elif indelcode == 'D':
                    seqa += '-'*nres
                    spos += 0
                elif indelcode == 'I':
                    spos += nres
                cpos = k+1 # advance to next chunk of cigar string
        assert(len(seqa) <= len(full_sequence))
        #TODO actually take into account signed_tmpl_len, rather than assuming aligned to beginning/end of sequence.
        if signed_tmpl_len < 0:
            seqa = '.'*(len(full_sequence)-len(seqa)) + seqa
        else:
            seqa = seqa + '.'*(len(full_sequence)-len(seqa))
        assert(len(seqa)==len(full_sequence))

        pos = {}
        pos['mut'] = []
        pos['del'] = []
        for n,nt_pair in enumerate(zip(full_sequence,seqa)):
            if nt_pair[0] not in 'ACGT-': continue
            if nt_pair[1] not in 'ACGT-': continue
            if nt_pair[1] == '-': pos['del'].append(n)
            elif nt_pair[0] != nt_pair[1]: pos['mut'].append(n)

        count += 1
        if count % 1000 == 0: print('Processed %d lines...' % count)

        total_mutdel = len(pos['mut']) + len(pos['del'])
        if total_mutdel <= args.mut_cutoff:
            count_filter += 1
            for tag1 in ['mut','del']:
                for tag2 in ['mut','del']:
                    for n1 in pos[tag1]:
                        for n2 in pos[tag2]:
                            counts_2d[tag1][tag2][n1][n2] += 1

        line = fid.readline()

    fid.close()
    print('Processed %d from %d lines\n' % (count_filter,count) )

for tag1 in ['mut','del']:
    for tag2 in ['mut','del']:
        outfile = out_file_tag + '.counts_%s_%s.txt' % (tag1,tag2)
        fid_out = open( outfile, 'w')
        for n1 in range(Nres):
            for n2 in range(Nres):
                fid_out.write('%d' % counts_2d[tag1][tag2][n1][n2])
                if n2 < Nres-1: fid_out.write(',')
            fid_out.write('\n')
        print('Wrote %dx%d to %s' % (Nres,Nres,outfile) )
        fid_out.close()

# 1D profiles
for tag in ['mut','del']:
    outfile = out_file_tag + '.counts_%s.txt' % (tag)
    fid_out = open( outfile, 'w')
    for n in range(Nres):
        fid_out.write('%d' % counts_2d[tag][tag][n][n])
        if n < Nres-1: fid_out.write(',')
    fid_out.write('\n')
    print('Wrote %d to %s' % (Nres,outfile) )
    fid_out.close()

# combine mut/del
counts_2d_all = []
for n1 in range( Nres ): counts_2d_all.append( [0]*Nres )
outfile = out_file_tag + '.counts_%s_%s.txt' % ('mutdel','mutdel')
fid_out = open( outfile, 'w')
for n1 in range(Nres):
    for n2 in range(Nres):
        for tag1 in ['mut','del']:
            for tag2 in ['mut','del']:
                counts_2d_all[n1][n2] += counts_2d[tag1][tag2][n1][n2]
        fid_out.write('%d' % counts_2d_all[n1][n2])
        if n2 < Nres-1: fid_out.write(',')
    fid_out.write('\n')
print('\nWrote %dx%d to %s' % (Nres,Nres,outfile) )
fid_out.close()

outfile = out_file_tag + '.counts_%s.txt' % ('mutdel')
fid_out = open( outfile, 'w')
for n in range(Nres):
    fid_out.write('%d' % counts_2d_all[n][n])
    if n < Nres-1: fid_out.write(',')
fid_out.write('\n')
print('Wrote %d to %s' % (Nres,outfile) )
fid_out.close()
