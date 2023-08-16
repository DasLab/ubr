#!/usr/bin/env python3
import argparse
import os
import time
import shutil

parser = argparse.ArgumentParser(
                    prog = 'ubr_run.py',
                    description = 'Gets mutation counts from FASTQs',
                    epilog = 'Runs Ultraplex, Bowtie2, RNA-framework with .csv output.\nRead 1 is assumed to be primer barcode, then reverse complement of RNA sequence.')

parser.add_argument('-s','--sequences_fasta', required=True, help='FASTA of RNA sequences')
parser.add_argument('-b','--primer_barcodes_fasta', required=True, help='FASTA of primer barcodes, first nucleotides of Read 1, prepended to cDNA')
parser.add_argument('-1','--read1_fastq', required=True, help='FASTQ (can be gzipped) of Illumina run')
parser.add_argument('-2','--read2_fastq', help='FASTQ (can be gzipped) of Read 2')
parser.add_argument('-ow','--overwrite',action = 'store_true', help='overwrite all previous files')
parser.add_argument('-O','--outdir',default='',help='output directory for all files')
parser.add_argument('-t','--threads',default=1, type=int, help='Number of threads for Bowtie2 and RNAFramework')
parser.add_argument('-mpb','--merge_pairs_bbmerge',action = 'store_true',help='Merge paired reads with bbmerge.sh')
parser.add_argument('-mpp','--merge_pairs_pear',action = 'store_true',help='Merge paired reads with PEAR')

parser.add_argument('-nm','--no_mixed',action = 'store_true',help='No mixed reads in Bowtie2')
parser.add_argument('-sm','--score_min',help='minimum score for Bowtie2')
parser.add_argument('-mq','--map_quality',default=10,type=int,help='minimum Bowtie2 MAPQ to consider read')

parser.add_argument('-lc','--length_cutoff',action = 'store_true',help='Use length cutoff of 0.92 length for RNAframework')
parser.add_argument('-nlc','--no_length_cutoff',action = 'store_true',help=argparse.SUPPRESS)
parser.add_argument('-norc','--no_output_raw_counts',action = 'store_true',help='do not output raw counts from RNAframework')
parser.add_argument('-orc','--output_raw_counts',action = 'store_true',help=argparse.SUPPRESS)
parser.add_argument('-me','--max_edit_distance',default=0.0,type=float,help='max edit distance for RNAFramework (0.15)')


args = parser.parse_args()
if args.no_length_cutoff: print( '--no_length_cutoff is on by default now! Flag will be deprecated later.' )
if args.output_raw_counts:  print( '--output_raw_counts is on by default now! Flag Will be deprecated later.' )
if (args.merge_pairs_bbmerge and args.merge_pairs_pear):
    print( '\nSpecify either --merge_pairs_bbmerge or --merge_pair_pear, not both')
    exit()

time_start = time.time()

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

# Check executables!
assert( shutil.which( 'ultraplex' ) )
assert( shutil.which( 'bowtie2-build' ) )
assert( shutil.which( 'bowtie2' ) )
assert( shutil.which( 'rf-count' ) )
assert( shutil.which( 'samtools' ) )
if args.merge_pairs_bbmerge:
    assert( shutil.which('bbmerge.sh') )
    assert( shutil.which('java') )
if args.merge_pairs_pear: assert( shutil.which('pear') )

assert( os.path.isfile( args.sequences_fasta ) )
assert( os.path.isfile( args.primer_barcodes_fasta ) )
assert( os.path.isfile( args.read1_fastq ) )
if args.read2_fastq: assert( os.path.isfile( args.read2_fastq ) )
read1_fastq = args.read1_fastq
read2_fastq = args.read2_fastq

# Read in FASTA files (sequences, and primer barcodes)
(sequences,headers) = read_fasta( args.sequences_fasta )
print( 'Read in %d sequences from %s.' % (len(sequences),args.sequences_fasta) )
(primer_barcodes,primer_names) = read_fasta( args.primer_barcodes_fasta )
print( 'Read in %d primer barcodes from %s.' % (len(primer_barcodes),args.primer_barcodes_fasta) )

# check sequences are RNA/DNA
def check_sequence(sequence):
    for c in sequence:
        if c not in 'ACGTU': return False
    return True

for sequence in sequences:
    if not check_sequence(sequence): exit('problem with sequence in sequences file: %s' % sequence )
for sequence in primer_barcodes:
    if not check_sequence(sequence): exit('problem with sequence in primer barcode file: %s' % sequence )

wd = args.outdir
if len(wd)>0:
    if wd[-1] != '/': wd += '/'
    if not os.path.isdir( wd ): os.makedirs( wd, exist_ok = True )

# Merge

if args.merge_pairs_bbmerge or args.merge_pairs_pear:
    out_prefix = args.read1_fastq.replace('.fq','').replace('.fastq','').replace('.gz','') + '_MERGED'
    merged_fastq = out_prefix+'.assembled.fastq'
    if os.path.isfile( merged_fastq ):
        print('Merged file already exists, skipping merge:',merged_fastq)
    else:
        if args.merge_pairs_bbmerge:
            command = 'bbmerge.sh in=%s in2=%s out=%s > %s0_merge_pairs.out 2> %s0_merge_pairs.err' % (args.read1_fastq, args.read2_fastq, merged_fastq, wd, wd)
        else:
            assert( args.merge_pairs_pear )
            command = 'pear -f %s -r %s -o  %s > %s0_merge_pairs.out 2> %s0_merge_pairs.err' % (args.read1_fastq, args.read2_fastq, out_prefix, wd, wd)
        print(command)
        os.system( command )
    assert(os.path.isfile( merged_fastq ) )
    read1_fastq = merged_fastq
    read2_fastq = None
time_merge_pairs = time.time()

# Ultraplex -- round 1 to demultiplex with respect to reverse transcription primers (e.g., RTB barcodes).
primer_barcodes_csv_file = wd + 'primer_barcodes.csv'
fid = open( primer_barcodes_csv_file, 'w' )
for barcode in primer_barcodes:  fid.write( barcode+'\n' )
fid.close()

outdir = wd + '1_ultraplex/'
os.makedirs(outdir,exist_ok = True)

# Did we already run ultraplex?
print()
any_ultraplex_out_files = False
for primer_barcode,primer_name in zip(primer_barcodes,primer_names):
    if read2_fastq:
        i1 = wd + '1_ultraplex/ultraplex_demux_5bc_%s_Fwd.fastq.gz'  % primer_barcode
        i2 = wd + '1_ultraplex/ultraplex_demux_5bc_%s_Rev.fastq.gz'  % primer_barcode
    else:
        i1 = wd + '1_ultraplex/ultraplex_demux_5bc_%s.fastq.gz'  % primer_barcode
        i2 = ''

    if os.path.isfile(i1) or os.path.isfile(i2):
        any_ultraplex_out_files = True
        break

if not any_ultraplex_out_files or args.overwrite:
    extra_flags = ''
    fastq_flags = ' -i %s' % (read1_fastq)
    if read2_fastq:  fastq_flags += ' -i2 %s' % (read2_fastq)
    command = 'ultraplex%s%s -b %s  -d %s  --dont_build_reference --ignore_no_match --threads %d > %s1_ultraplex/1_ultraplex.out 2> %s1_ultraplex/1_ultraplex.err' % (fastq_flags, extra_flags, primer_barcodes_csv_file, outdir, args.threads, wd, wd)
    print(command)
    os.system( command )
else:
    print('\nSkipping %s' % outdir)
time_ultraplex1 = time.time()


# Bowtie2
print()
bowtie_build_dir = wd + '2_bowtie2/bowtie-build'
os.makedirs(bowtie_build_dir,exist_ok = True)
bt2_prefix = '%s/seq'% (bowtie_build_dir)
bt2_file = '%s.1.bt2' % bt2_prefix

# Need to make sure we have DNA versions
seq_file = wd + 'seq.fasta'
fid = open( seq_file, 'w' )
for (sequence,header) in zip(sequences, headers): fid.write('>%s\n%s\n' % (header,sequence.replace('U','T')) )
fid.close()

if not os.path.isfile(bt2_file) or args.overwrite:
    command = 'bowtie2-build %s %s --threads %d > %s/bt2.out 2> %s/bt2.err' % (seq_file,bt2_prefix,args.threads,bowtie_build_dir,bowtie_build_dir)
    print(command)
    os.system( command )

for primer_barcode,primer_name in zip(primer_barcodes,primer_names):
    if read2_fastq:
        i1 = wd + '1_ultraplex/ultraplex_demux_5bc_%s_Fwd.fastq.gz'  % primer_barcode
        i2 = wd + '1_ultraplex/ultraplex_demux_5bc_%s_Rev.fastq.gz'  % primer_barcode
        if not os.path.isfile(i1): continue
        if not os.path.isfile(i2): continue
        fastq_flags = ' -1 %s -2 %s' % (i1,i2)
    else:
        i1 = wd + '1_ultraplex/ultraplex_demux_5bc_%s.fastq.gz'  % primer_barcode
        if not os.path.isfile(i1): continue
        fastq_flags = ' -U %s' % i1

    outdir = wd + '2_bowtie2/%s' % primer_name
    os.makedirs(outdir,exist_ok = True)

    # bowtie2-align
    sam_file = '%s/bowtie2.sam' % outdir
    extra_flags = ''
    if args.no_mixed: extra_flags += ' --no-mixed'
    if args.score_min != None:  extra_flags += ' --score-min %s' % args.score_min
    if not os.path.isfile( sam_file ) or args.overwrite:
        # previously used --local --sensitive-local, but bowtie2 would punt on aligning 3' ends and misalign reads to some short parasite replicons.
        command = 'bowtie2 --end-to-end --sensitive --maxins=800 --ignore-quals --no-unal --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 -x %s %s -S %s --threads %d %s > %s/bowtie2.out 2> %s/bowtie2.err' % (bt2_prefix,fastq_flags,sam_file,args.threads,extra_flags,outdir,outdir)
        print(command)
        os.system( command )
    else:
        print( 'Skipping bowtie2-align for %s' % outdir )

time_bowtie2 = time.time()

# Some deprecated options
if args.length_cutoff:
    min_seq_length = min( map( lambda x:len(x), sequences ) )
    MIN_READ_LENGTH = int(min_seq_length * 0.92)
    print( '\nrf-count will throw out reads with length smaller than 92%% of minimal sequence length: %d' % MIN_READ_LENGTH )
else:
    MIN_READ_LENGTH = 1

# rf-count to assign muts/dels/inserts
print()
rf_count_dir = wd + '3_rf_count'
os.makedirs(rf_count_dir,exist_ok = True)
#rf_count_outfile = rf_count_dir+'/rf-count.out'
#rf_count_errfile = rf_count_dir+'/rf-count.err'
for (n,primer_name) in enumerate(primer_names):
    bowtie_align_dir = wd + '2_bowtie2/%s/' % (primer_name)
    sam_file = '%s/bowtie2.sam' % bowtie_align_dir
    if not os.path.isfile( sam_file ): continue

    outdir = '%s/%s' % (rf_count_dir,primer_name)
    if args.overwrite or not os.path.isfile(outdir+'/bowtie2.rc'):
        # need to remove dir and start job; rf-count's overwrite option does crazy things.
        if os.path.isdir(outdir): shutil.rmtree(outdir)
        rf_count_outfile = rf_count_dir+'/rf-count_%s.out' % primer_name
        rf_count_errfile = rf_count_dir+'/rf-count_%s.err' % primer_name
        if os.path.isfile(rf_count_outfile): os.remove( rf_count_outfile )
        if os.path.isfile(rf_count_errfile): os.remove( rf_count_errfile )
        extra_flags = ''
        if not args.no_output_raw_counts: extra_flags = ' -orc '
        if args.map_quality != 10:  extra_flags += ' --map-quality %d' % args.map_quality
        if args.max_edit_distance > 0:  extra_flags += ' --max-edit-distance %f' % args.max_edit_distance

        command = 'rf-count --processors %d -wt 1 -fast -f %s -m -cc -rd -ni -ds %d %s -o %s %s >> %s 2>> %s' % (args.threads, seq_file, MIN_READ_LENGTH, extra_flags, outdir, sam_file, rf_count_outfile, rf_count_errfile)
        print(command)
        os.system( command )
    else:
        print( 'Skipping rf-count into: %s' % outdir )

time_rf_count = time.time()

# rf-rctools view to create an easy to parse .csv
print()
for primer_name in primer_names:
    rc_file = wd + '3_rf_count/%s/bowtie2.rc' % (primer_name)

    outdir = wd + '4_rctools/%s' % (primer_name)
    rf_count_file = outdir + '/rf_count.csv'
    if args.overwrite or not os.path.isfile(rf_count_file):
        os.makedirs(outdir,exist_ok = True)
        command = 'rf-rctools view %s > %s' % (rc_file, rf_count_file)
        print(command)
        os.system( command )
        assert( len(open(rf_count_file).readlines()) == 5 * len(sequences) )
    else:
        print( 'Skipping rf-rctools into: %s' % outdir )
time_rctools = time.time()

# Compile information into final .txt files.
print()
Npos = 0
for primer_name in primer_names:
    outfile_counts = wd + '%s.muts.txt' % primer_name
    outfile_coverage = wd + '%s.coverage.txt' % primer_name
    fid_counts = open(outfile_counts,'w')
    fid_coverage = open(outfile_coverage,'w')

    N = 0
    infile = wd + '4_rctools/%s/rf_count.csv' % (primer_name)
    total_coverage = 0
    if os.path.isfile( infile ):
        lines = open( infile ).readlines()
        N = int(len(lines)/5)
        all_muts = []
        all_coverage = []
        tally_coverage = (N <= 50000) # takes long time for long libraries
        for n in range( N ):
            muts     = lines[5*n+2].strip('\n')
            coverage = lines[5*n+3].strip('\n')
            fid_counts.write( muts+'\n' )
            fid_coverage.write( coverage + '\n' )
            if n == 0: Npos = len(coverage.split(','))
            if tally_coverage: total_coverage += max(map( lambda x:int(x), coverage.split(',') ))
    else:
        print( 'WARNING! Could not find %s' % infile )
    fid_counts.close()
    fid_coverage.close()
    total_coverage_string = ''
    if tally_coverage: total_coverage_string = ' with total coverage %d' % total_coverage
    print( 'Created: %s and %s for %d sequences%s' % (outfile_counts,outfile_coverage,N,total_coverage_string) )

# Compile information on mutation-type read counts (if available)
print()
design_name_idx = {}
for (idx,header) in enumerate(headers): design_name_idx[ header.strip().split()[0] ] = idx
assert( len(design_name_idx) == N )

mut_types = ['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del']
outdir = wd+'raw_counts/'
os.makedirs(outdir,exist_ok = True)
for primer_name in primer_names:
    infile = wd + '3_rf_count/%s/raw_counts/bowtie2.txt' % (primer_name)
    if os.path.isfile( infile ):
        lines = open( infile ).readlines()
        N = int(len(lines)/16)
        outfiles_raw_counts = []
        for (k,mut_type) in enumerate(mut_types):
            raw_count_lines = [','.join( ['0']*Npos )]*len(design_name_idx)
            for n in range( N ):
                design_name = lines[16*n].strip('\n')
                assert( design_name in design_name_idx )
                idx = design_name_idx[ design_name ]
                cols = lines[16*n+1+k].strip('\n').split()
                assert( cols[0] == mut_type)
                raw_count_lines[idx] = cols[1]
            outfile_raw_counts = outdir + '%s.%s.txt' % (primer_name,mut_type)
            outfiles_raw_counts.append(outfile_raw_counts)
            fid_raw_counts = open(outfile_raw_counts,'w')
            for line in raw_count_lines: fid_raw_counts.write( line+'\n' )
            fid_raw_counts.close()
        print( 'Created %s for %d sequences (found %d designs)' % (','.join(outfiles_raw_counts),len(design_name_idx),N) )
    else:
        print( 'WARNING! Could not find %s' % infile )

time_end=time.time()

print( '\nTimings:')
print( '0_merge_pairs ' + time.strftime("%H:%M:%S",time.gmtime(time_merge_pairs-time_start) ) )
print( '1_ultraplex   ' + time.strftime("%H:%M:%S",time.gmtime(time_ultraplex1-time_merge_pairs) ) )
print( '2_bowtie2     ' + time.strftime("%H:%M:%S",time.gmtime(time_bowtie2-time_ultraplex1) ) )
print( '3_rf_count    ' + time.strftime("%H:%M:%S",time.gmtime(time_rf_count-time_bowtie2) ) )
print( '4_rctools     ' + time.strftime("%H:%M:%S",time.gmtime(time_rctools-time_rf_count) ) )
print( 'Final merge   ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_rctools) ) )

print( '\nTotal time: ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_start) ) )

