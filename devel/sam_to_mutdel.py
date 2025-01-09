#! /usr/bin/env python3
import argparse
parser = argparse.ArgumentParser(
                    prog = 'sam_to_mutdel.py',
                    description = 'Processes SAM file to create simple mut/del file with 0,1,...16 codes',
                    epilog = 'Extremely basic processing of CIGAR string in bowtie2 SAM file to produce alignments and then coding.')

parser.add_argument('--sam', required=True,default='RTB010_CustomArray_PK50_1M7_BindTheFivePrimeEnd.sam.txt',help='SAM-formatted text file, Bowtie2')
parser.add_argument('-s','--sequences_fasta', required=True, default='pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa',help='Fasta file with DNA/RNA sequences used for Bowtie2')

args = parser.parse_args()
samfile  = args.sam
assert( samfile.find('.sam')>-1 )
outfile = samfile.replace('.sam','.mutdel')
outfile_align = samfile.replace('.sam','.aln')

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

(sequences,headers) = read_fasta( args.sequences_fasta )

fid_out = open(outfile,'w')
fid_aln = open(outfile_align,'w')
fid = open(samfile)
line = fid.readline()
count = 0

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

cols = line.split()
(full_sequence,id) = figure_out_sequence(cols[2])

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
    fid_aln.write( '%s\n' % seqa )
    assert(len(seqa)==len(full_sequence))

    nt_mut_codes = {'AA':0,'AC':1,'AG':2,'AT':3,'CA':4,'CC':0,'CG':5,'CT':6,
                    'GA':7,'GC':8,'GG':0,'GT':9,'TA':10,'TC':11,'TG':12,'TT':0,
                    'A-':13,'C-':14,'G-':15,'T-':16}
    codes_out = ['0']*len(full_sequence)
    for n,nt_pair in enumerate(zip(full_sequence,seqa)):
        if nt_pair[0] not in 'ACGT-': continue
        if nt_pair[1] not in 'ACGT-': continue
        nt_mut = nt_pair[0]+nt_pair[1]
        codes_out[n] = str(nt_mut_codes[nt_mut])

    fid_out.write(' '.join(codes_out))
    fid_out.write('\n')

    line = fid.readline()
    count += 1
    if count % 1000 == 0: print('Processed %d lines...' % count)

fid.close()
fid_out.close()
fid.close()
fid_aln.close()
print('Wrote %d lines into %s\n' % (count,outfile) )

print('Wrote %d lines into %s\n' % (count,outfile_align) )


