#!/usr/bin/env python3
import argparse
import os
import time
import shutil
import gzip
import glob
import random
import ubr_check_stats
from ubr_util import read_fasta,check_sequence,create_seq_fasta

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
parser.add_argument('-nmp','--no_merge_pairs',action = 'store_true',help='do not merge paired end reads before Bowtie2' )
parser.add_argument('--cmuts',action = 'store_true',help='use cmuts instead of RNAFramework' )
parser.add_argument('--precomputed_bowtie_build_dir',default='',help=argparse.SUPPRESS) # precomputed bowtie-build directory, which otherwise takes forever to generate on the fly for >1M seqs
parser.add_argument('-f','--force',action='store_true',help='override some warnings' )
parser.add_argument('--dedup',action='store_true',help='dedup if N''s in barcode' )

# Deprecated/secret
parser.add_argument('--excise_barcode',default=0,type=int,help=argparse.SUPPRESS) # 'remove this many nucleotides from merged FASTQ and sequences' )
parser.add_argument('-nm','--no_mixed',action = 'store_true',help=argparse.SUPPRESS )# 'No mixed reads in Bowtie2')
parser.add_argument('-sm','--score_min',help=argparse.SUPPRESS )#'minimum score for Bowtie2')
parser.add_argument('-mq','--map_quality',default=10,type=int,help=argparse.SUPPRESS )#help='minimum Bowtie2 MAPQ to consider read')
parser.add_argument('-lc','--length_cutoff',action = 'store_true',help=argparse.SUPPRESS )#help='Use length cutoff of 0.92 length for RNAFramework')
parser.add_argument('-norc','--no_output_raw_counts',action = 'store_true',help=argparse.SUPPRESS )#help='do not output raw counts from RNAFramework')
parser.add_argument('-me','--max_edit_distance',default=0.0,type=float,help=argparse.SUPPRESS )#help='max edit distance for RNAFramework (0.15)')
parser.add_argument('-mpp','--merge_pairs_pear',action = 'store_true',help=argparse.SUPPRESS)
parser.add_argument('--force_merge_pairs',action = 'store_true',help=argparse.SUPPRESS) # force merge pairs (don't bother to check for overlap)
parser.add_argument('--skip_ultraplex',action = 'store_true',help=argparse.SUPPRESS)
parser.add_argument('--cutadapt',action = 'store_true',help=argparse.SUPPRESS) # force cutadapt trimming of Read2 side for pre-demuxed Ultima
parser.add_argument('--no_collapse',action = 'store_true',help=argparse.SUPPRESS) # no collapse option in rf-count
parser.add_argument('--use_tmp_dir',action = 'store_true',help=argparse.SUPPRESS) # For cmuts, run job in /tmp/ to try to reduce disk i/o
parser.add_argument('--ultima',action='store_true',help=argparse.SUPPRESS) # recognize Ultima adapter in ultraplex
parser.add_argument('--minimap2',action='store_true',help=argparse.SUPPRESS) # use minimap2 instead of bowtie2

args = parser.parse_args()
assert( not( args.no_merge_pairs and args.merge_pairs_pear ) )

time_start = time.time()

# Check executables!
if not args.skip_ultraplex: assert( shutil.which( 'ultraplex' ) )
if args.minimap2:
    assert( shutil.which( 'minimap2' ) )
else:
    assert( shutil.which( 'bowtie2-build' ) )
    assert( shutil.which( 'bowtie2' ) )
if args.cmuts: assert( shutil.which( 'cmuts' ))
else: assert( shutil.which( 'rf-count' ) )
assert( shutil.which( 'samtools' ) )
assert( shutil.which( 'bbmerge.sh' ) )
assert( shutil.which( 'java' ) )
if args.merge_pairs_pear: assert( shutil.which('pear') )
if args.cutadapt: assert( shutil.which( 'cutadapt' ) )
if args.dedup: assert( shutil.which('umi_tools') )

assert( os.path.isfile( args.sequences_fasta ) )
assert( os.path.isfile( args.primer_barcodes_fasta ) )
assert( os.path.isfile( args.read1_fastq ) )
if args.read2_fastq: assert( os.path.isfile( args.read2_fastq ) )
read1_fastq = args.read1_fastq
read2_fastq = args.read2_fastq

# Read in FASTA files (sequences, and primer barcodes)
(sequences,headers) = read_fasta( args.sequences_fasta, args.force )
print( 'Read in %d sequences from %s.' % (len(sequences),args.sequences_fasta) )
(primer_barcodes,primer_names) = read_fasta( args.primer_barcodes_fasta, args.force )
print( 'Read in %d primer barcodes from %s.' % (len(primer_barcodes),args.primer_barcodes_fasta) )
for primer_name in primer_names: assert( primer_name.find(' ')==-1 )

if args.ultima: primer_barcodes = [ 'CTACACGACGCTCTTCCGATCT'+barcode for barcode in primer_barcodes ]

if len(sequences)>1000000:
    if not args.precomputed_bowtie_build_dir:
        print( '\nYou have a lot of sequences. It is recommended to pre-index bowtie2 build with:\n\n ubr_util.py --threads 24 --build_bowtie2_index -s %s\n\nThen re-run this script with flag --precomputed_bowtie_build_dir bowtie-build.\n(To bypass this message, use --force)\n' % args.sequences_fasta )
        if not args.force: exit()
    if not args.cmuts:
        print( '\nYou have a lot of sequences. It is recommended to use cmuts with the flag --cmuts.\n(To bypass this message, use --force)\n')
        if not args.force: exit()

print('Checking sequences...')
for sequence in sequences:
    if not check_sequence(sequence): exit('problem with sequence in sequences file: %s' % sequence )
for sequence in primer_barcodes:
    if not check_sequence(sequence): exit('problem with sequence in primer barcode file: %s' % sequence )

wd = args.outdir
if len(wd)>0:
    if wd[-1] != '/': wd += '/'
    if not os.path.isdir( wd ): os.makedirs( wd, exist_ok = True )

# excise barcodes [testing if barcodes are important]
if args.excise_barcode > 0:
    new_sequence_file = args.sequences_fasta.replace('.fa','.NO_BARCODES.fa')
    fid = open( new_sequence_file, 'w' )
    for (i,(sequence,header)) in enumerate(zip(sequences,headers)):
        sequence = sequence[:(-args.excise_barcode)]
        fid.write('>%s\n%s\n' % (header,sequence))
        sequences[i] = sequence
    print('\nEXCISE_BARCODE: Outputted %d sequences with last %d nucleotides excised in %s\n' % (len(sequences),args.excise_barcode,new_sequence_file))

# Merge
merge_pairs = not args.no_merge_pairs and args.read2_fastq != None

# Also check if sequence reads are actually long enough to overlap
if merge_pairs and not args.merge_pairs_pear and not args.force_merge_pairs:
    insert_len = len(sequences[0]) + len(primer_barcodes[0])
    read_len1 = len(os.popen( 'seqkit head -n 1 %s' % args.read1_fastq ).readlines()[1].strip())
    read_len2 = len(os.popen( 'seqkit head -n 1 %s' % args.read2_fastq ).readlines()[1].strip())
    overlap = read_len1 + read_len2 - insert_len
    merge_pairs = overlap >= 12 #minimum overlap for bbmerge.sh
    if not merge_pairs: print( 'Will not merge pairs due to insufficient overlap of read1 and read2!')
    #print( 'Insert len: ',insert_len,' Read1_len: ',read_len1,'Read2_len:',read_len2,'overlap: ',overlap,' merge_pairs? ',merge_pairs)
    #exit()

if merge_pairs:
    out_prefix = os.path.basename(args.read1_fastq).replace('.fq','').replace('.fastq','').replace('.gz','') + '_MERGED'
    merged_fastq = wd+out_prefix+'.assembled.fastq.gz'
    if os.path.isfile( merged_fastq ):
        print('Merged file already exists, skipping merge:',merged_fastq)
    else:
        if args.merge_pairs_pear:
            command = 'pear -f %s -r %s -o  %s > %s0_merge_pairs.out 2> %s0_merge_pairs.err && gzip %s' % (args.read1_fastq, args.read2_fastq, out_prefix, wd, wd, merged_fastq.replace('.gz',''))
        else:
            # Default is to use bbmerge.sh
            # turn off pigz/unpigz which does not accelerate things but spins up extra threads.
            command = 'bbmerge.sh in=%s in2=%s out=%s pigz=f unpigz=f > %s0_merge_pairs.out 2> %s0_merge_pairs.err' % (args.read1_fastq, args.read2_fastq, merged_fastq, wd, wd)
        print(command)
        os.system( command )
    assert(os.path.isfile( merged_fastq ) )
    read1_fastq = merged_fastq
    read2_fastq = None
time_merge_pairs = time.time()

# Ultraplex -- round 1 to demultiplex with respect to reverse transcription primers (e.g., RTB barcodes).
primer_barcodes_csv_file = wd + 'primer_barcodes.csv'
fid = open( primer_barcodes_csv_file, 'w' )
for (primer_name,barcode) in zip(primer_names,primer_barcodes):  fid.write( barcode+':'+primer_name+'\n' )
fid.close()

# Don't run dedup
dedup = args.dedup
barcodes_have_N = False
for barcode in primer_barcodes:
    if barcode.find('N')>-1:  barcodes_have_N = True
if dedup and not barcodes_have_N:
    print('\nAsked for --dedup but no N''s in any barcodes. Will not deduplicate.')
    dedup = False

outdir = wd + '1_ultraplex/'
os.makedirs(outdir,exist_ok = True)

# Did we already run ultraplex?
print()
any_ultraplex_out_files = False
for primer_barcode,primer_name in zip(primer_barcodes,primer_names):
    if read2_fastq and not merge_pairs:
        i1 = wd + '1_ultraplex/ultraplex_demux_%s_Fwd.fastq.gz'  % primer_name
        i2 = wd + '1_ultraplex/ultraplex_demux_%s_Rev.fastq.gz'  % primer_name
    else:
        i1 = wd + '1_ultraplex/ultraplex_demux_%s.fastq.gz'  % primer_name
        i2 = ''

    if os.path.isfile(i1) or os.path.isfile(i2):
        any_ultraplex_out_files = True
        break

# has FASTQ file already been demultiplexed? Look for the barcode string like ACCAGGCGCTGG in the filename.
skip_ultraplex = args.skip_ultraplex
for (primer_barcode,primer_name) in zip(primer_barcodes,primer_names):
    if read1_fastq.find( primer_barcode ) > -1 or read1_fastq.find( primer_name ) > -1:
        print('Detected primer %s (%s) in name of FASTQ file %s. Assuming FASTQ has already been demultiplexed.' % (primer_name,primer_barcode,read1_fastq))
        i1 = wd + '1_ultraplex/ultraplex_demux_%s.fastq.gz'  % primer_name
        if args.cutadapt:
            command = 'cutadapt --trimmed-only %s -a  AGATCGGAAGAGCACA -o %s' % (read1_fastq,i1)
        else:
            command = 'rsync -avL %s %s' % (read1_fastq,i1)
        print(command)
        os.system(command)
        skip_ultraplex = True
        assert( not read2_fastq )

if args.cutadapt and not skip_ultraplex: print("--cutadapt option is only for files that have been pre-demultiplexed by ultima.")

if not skip_ultraplex and (not any_ultraplex_out_files or args.overwrite):
    extra_flags = ''
    fastq_flags = ' -i %s' % (read1_fastq)
    if read2_fastq:  fastq_flags += ' -i2 %s' % (read2_fastq)
    if args.ultima: extra_flags +=' -a AGATCGGAAGAGCACA'
    command = 'ultraplex%s%s -b %s  -d %s  --dont_build_reference --ignore_no_match --threads %d > %s1_ultraplex/1_ultraplex.out 2> %s1_ultraplex/1_ultraplex.err' % (fastq_flags, extra_flags, primer_barcodes_csv_file, outdir, args.threads, wd, wd)
    print(command)
    os.system( command )
else:
    print('\nSkipping %s' % outdir)
time_ultraplex1 = time.time()

# Excise barcode segments from FASTQ
if args.excise_barcode:
    assert( merge_pairs )
    for primer_barcode,primer_name in zip(primer_barcodes,primer_names):
        i1 = wd + '1_ultraplex/ultraplex_demux_%s.fastq.gz'  % primer_name
        i1_save = i1.replace('.fastq.gz','.SAVE.fastq.gz')
        if not os.path.isfile( i1_save ):
            assert( os.path.isfile( i1 ) )
            shutil.move( i1, i1_save )
            command = 'seqkit subseq -r %d:-1 %s -o %s' % ((args.excise_barcode+1),i1_save,i1)
            print(command)
            os.system( command )

# Bowtie2
print()
bowtie_build_dir = wd + '2_bowtie2/bowtie-build'
os.makedirs(bowtie_build_dir,exist_ok = True)
bt2_prefix = '%s/seq'% (bowtie_build_dir)
bt2_file = '%s.1.bt2' % bt2_prefix

# Need to make sure we have DNA versions
seq_file = create_seq_fasta( sequences, headers, wd )

for primer_barcode,primer_name in zip(primer_barcodes,primer_names):
    if read2_fastq:
        i1 = wd + '1_ultraplex/ultraplex_demux_%s_Fwd.fastq.gz'  % primer_name
        i2 = wd + '1_ultraplex/ultraplex_demux_%s_Rev.fastq.gz'  % primer_name
        if not os.path.isfile(i1): continue
        if not os.path.isfile(i2): continue
        fastq_flags = ' -1 %s -2 %s' % (i1,i2)
    else:
        i1 = wd + '1_ultraplex/ultraplex_demux_%s.fastq.gz'  % primer_name
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
        if (not os.path.isfile(bt2_file) or args.overwrite) and not args.minimap2:
            if len(args.precomputed_bowtie_build_dir)>0:
                #os.symlink( os.path.abspath(args.precomputed_bowtie_build_dir), bowtie_build_dir,target_is_directory=True )
                #print( 'Symlinked %s to %s.' % (args.precomputed_bowtie_build_dir, bowtie_build_dir) )
                command = 'rm -rf %s && ln -fs %s %s' % (bowtie_build_dir, os.path.abspath(args.precomputed_bowtie_build_dir),bowtie_build_dir)
            else:
                command = 'bowtie2-build %s %s --threads %d > %s/bt2.out 2> %s/bt2.err' % (seq_file,bt2_prefix,args.threads,bowtie_build_dir,bowtie_build_dir)
            print(command)
            os.system( command )

        if args.minimap2:
            # a bit of a hack -- really should change names of directories and files to be minimap2 instead of bowtie2
            # note: minimap2 also allows precompilation of index and special options customized for illumina or pacbio -- would be worth trying.
            assert(not read2_fastq)
            minimap2_sam_file = '%s/minimap2.sam' % outdir
            command = 'minimap2 -c -a %s %s > %s  2> %s/minimap2.err' % (seq_file,i1,minimap2_sam_file,outdir)
            print(command)
            os.system( command )

            # need CIGAR strings, so use samtools calmd and output to 'bowtie2.sam'
            command = 'samtools calmd %s %s > %s 2> %s/calmd.err' % (minimap2_sam_file,seq_file,sam_file,outdir)
            print(command)
            os.system( command )
        else:
            # previously used --local --sensitive-local, but bowtie2 would punt on aligning 3' ends and misalign reads to some short parasite replicons.
            command = 'bowtie2 --end-to-end --sensitive --maxins=800 --ignore-quals --no-unal --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 -x %s %s -S %s --threads %d %s > %s/bowtie2.out 2> %s/bowtie2.err' % (bt2_prefix,fastq_flags,sam_file,args.threads,extra_flags,outdir,outdir)
            print(command)
            os.system( command )
    else:
        print( 'Skipping bowtie2-align for %s' % outdir )

# Run sorts -- useful for all later steps.
print()
for (n,primer_name) in enumerate(primer_names):
    bowtie_align_dir = wd + '2_bowtie2/%s' % (primer_name)
    sam_file = '%s/bowtie2.sam' % bowtie_align_dir
    if not os.path.isfile( sam_file ): continue

    sorted_bam_file = '%s/bowtie2.sorted.bam' % bowtie_align_dir
    if args.overwrite or not os.path.isfile( sorted_bam_file ):
        command = 'samtools sort %s -o %s > %s/sort.log 2> %s/sort.err' % (sam_file, sorted_bam_file, bowtie_align_dir, bowtie_align_dir)
        print(command)
        os.system( command )

if dedup:
    for (n,primer_name) in enumerate(primer_names):
        bowtie_align_dir = wd + '2_bowtie2/%s' % (primer_name)
        sorted_bam_file = '%s/bowtie2.sorted.bam' % bowtie_align_dir
        if not os.path.isfile( sam_file ): continue

        sorted_bam_index_file = '%s.bai' % sorted_bam_file
        if args.overwrite or not os.path.isfile( sorted_bam_index_file ):
            command = 'samtools index %s > %s/bai.log 2> %s/bai.err' % (sorted_bam_file,bowtie_align_dir,bowtie_align_dir)
            print(command)
            os.system( command )

        dedup_bam_file = '%s/bowtie2.sorted.dedup.bam' % bowtie_align_dir
        if args.overwrite or not os.path.isfile( dedup_bam_file ):
            dedup_log_file = '%s/bowtie2.sorted.dedup.log' % bowtie_align_dir
            command = "umi_tools dedup --stdin=%s --log=%s --umi-separator=':'  > %s" % (sorted_bam_file,dedup_log_file,dedup_bam_file)
            print(command)
            os.system( command )


time_bowtie2 = time.time()

if args.cmuts:
    cmuts_dir = wd + '3_cmuts'
    os.makedirs(cmuts_dir,exist_ok = True)
    print()

    for (n,primer_name) in enumerate(primer_names):
        bowtie_align_dir = wd + '2_bowtie2/%s/' % (primer_name)
        outdir = '%s/%s' % (cmuts_dir,primer_name)
        os.makedirs(outdir,exist_ok = True)

        sorted_bam_file = '%s/bowtie2.sorted.bam' % bowtie_align_dir
        if dedup: sorted_bam_file = '%s/bowtie2.sorted.dedup.bam' % bowtie_align_dir

        sorted_bam_index_file = '%s.bai' % sorted_bam_file
        if args.overwrite or not os.path.isfile( sorted_bam_index_file ):
            command = 'samtools index %s > %s/bai.log 2> %s/bai.err' % (sorted_bam_file,bowtie_align_dir,bowtie_align_dir)
            print(command)
            os.system( command )

        fasta_index_file = '%s.fai' % seq_file
        if args.overwrite or not os.path.isfile( fasta_index_file ):
            command = 'samtools faidx %s > %s/faidx.log 2> %s/faidx.err' % (seq_file,bowtie_align_dir,bowtie_align_dir)
            print(command)
            os.system( command )

        hdf5_file = wd + '%s.hdf5' % primer_name
        hdf5_outfile = hdf5_file
        lockfiles = glob.glob( hdf5_file+'-*.lock' )
        if args.overwrite or not os.path.isfile( hdf5_file ) or len( lockfiles ) > 0:
            for lockfile in lockfiles: os.remove( lockfile )
            if os.path.isfile( hdf5_file ): os.remove( hdf5_file )
            if args.use_tmp_dir: # output to /tmp/ to mitigate cost of disk output on cluster
                rand_dir = '/tmp/%d' % random.getrandbits(64)
                if os.path.isdir( rand_dir ): shutil.rmtree( rand_dir )
                os.makedirs(rand_dir,exist_ok = True)
                hdf5_outfile = '%s/%s.hdf5' % (rand_dir,primer_name)
                # also consider rsync of sorted_bam_file to /tmp/ to reduce *input* from disk
            command = 'cmuts --fasta=%s --overwrite  --min-cov-base-quality=0 --output %s %s > %s.cmuts.log 2> %s.cmuts.err' % (seq_file,hdf5_outfile,sorted_bam_file,outdir,outdir)
            print(command)
            os.system( command )
            if args.use_tmp_dir:
                if os.path.isfile( hdf5_outfile ):
                    command = 'rsync -avL %s %s' % (hdf5_outfile, hdf5_file)
                    print(command)
                    os.system(command)
                else:
                    printf('Job failed! Could not find: %s!' % hdf5_outfile )
                    continue
                shutil.rmtree(rand_dir)


else: # use RNAFramework
    # Some deprecated options
    if args.length_cutoff:
        min_seq_length = min( map( lambda x:len(x), sequences ) )
        MIN_READ_LENGTH = int(min_seq_length * 0.92)
        print( '\nrf-count will throw out reads with length smaller than 92%% of minimal sequence length: %d' % MIN_READ_LENGTH )
    else:
        MIN_READ_LENGTH = 1
    max_seq_length = max( map( lambda x:len(x), sequences ) )

    # rf-count to assign muts/dels/inserts
    print()
    rf_count_dir = wd + '3_rf_count'
    os.makedirs(rf_count_dir,exist_ok = True)
    bowtie2_tag = 'bowtie2.sorted'
    if dedup:  bowtie2_tag = 'bowtie2.sorted.dedup'
    for (n,primer_name) in enumerate(primer_names):
        bowtie_align_dir = wd + '2_bowtie2/%s/' % (primer_name)
        bam_file = '%s/%s.bam' % (bowtie_align_dir,bowtie2_tag)
        if not os.path.isfile( bam_file ): continue

        outdir = '%s/%s' % (rf_count_dir,primer_name)
        if args.overwrite or not os.path.isfile(outdir+'/'+bowtie2_tag+'.rc'):
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
            cc_option = ' -cc'
            if args.no_collapse: cc_option = ''

            command = 'rf-count --processors %d -r -wt 1 -fast -f %s -m%s -rd -ni -ds %d %s -o %s %s >> %s 2>> %s' % (args.threads, seq_file, cc_option, MIN_READ_LENGTH, extra_flags, outdir, bam_file, rf_count_outfile, rf_count_errfile)
            #command = 'rf-count --processors 1 --working-threads %d -fast -f %s -m -cc -rd -ni -ds %d %s -o %s %s >> %s 2>> %s' % (args.threads, seq_file, MIN_READ_LENGTH, extra_flags, outdir, sam_file, rf_count_outfile, rf_count_errfile)
            print(command)
            os.system( command )

            if not args.no_output_raw_counts:
                raw_counts_file = rf_count_dir+'/%s/raw_counts/%s.txt' % (primer_name,bowtie2_tag)
                assert( os.path.isfile(raw_counts_file) )
                command = 'gzip %s' % raw_counts_file
                print(command)
                os.system( command )

        else:
            print( 'Skipping rf-count into: %s' % outdir )

    time_rf_count = time.time()

    # rf-rctools view to create an easy to parse .csv
    print()
    for primer_name in primer_names:
        rc_file = wd + '3_rf_count/%s/%s.rc' % (primer_name,bowtie2_tag)

        outdir = wd + '4_rctools/%s' % (primer_name)
        rf_count_file = outdir + '/rf_count.csv'
        rf_count_file_gz = rf_count_file + '.gz'
        if os.path.isfile(rc_file) and (args.overwrite or not os.path.isfile(rf_count_file_gz)):
            os.makedirs(outdir,exist_ok = True)
            command = 'rf-rctools view %s > %s && gzip %s' % (rc_file, rf_count_file, rf_count_file)
            print(command)
            os.system( command )
            assert( len(gzip.open(rf_count_file_gz).readlines()) == 5 * len(sequences) )
        else:
            print( 'Skipping rf-rctools into: %s' % outdir )
    time_rctools = time.time()

    # Compile information into final .txt files.
    print()
    Npos = 0
    N = 0
    for primer_name in primer_names:
        outfile_counts = wd + '%s.muts.txt.gz' % primer_name
        outfile_coverage = wd + '%s.coverage.txt.gz' % primer_name
        fid_counts = gzip.open(outfile_counts,'wt')
        fid_coverage = gzip.open(outfile_coverage,'wt')

        infile = wd + '4_rctools/%s/rf_count.csv.gz' % (primer_name)
        total_coverage = 0
        tally_coverage = False
        if os.path.isfile( infile ):
            lines = gzip.open( infile, 'rt' ).readlines()
            N = int(len(lines)/5)
            all_muts = []
            all_coverage = []
            tally_coverage = (N <= 50000) # takes long time for long libraries
            for n in range( N ):
                muts     = lines[5*n+2].strip('\n')
                coverage = lines[5*n+3].strip('\n')
                pad_string = ',0'*(max_seq_length-len(sequences[n]))
                fid_counts.write( muts + pad_string +'\n' )
                fid_coverage.write( coverage + pad_string + '\n' )
                if tally_coverage: total_coverage += max(list(map( lambda x:int(x), coverage.split(',') )))
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
    for (idx,header) in enumerate(headers): design_name_idx[ header.strip().split()[0].replace('/','_') ] = idx
    assert( len(design_name_idx) == N or N == 0)

    mut_types = ['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del']
    outdir = wd+'raw_counts/'
    os.makedirs(outdir,exist_ok = True)
    for primer_name in primer_names:
        infile = wd + '3_rf_count/%s/raw_counts/%s.txt.gz' % (primer_name,bowtie2_tag)
        if os.path.isfile( infile ):
            lines = gzip.open(infile,'rt').readlines()
            N = int(len(lines)/16)
            outfiles_raw_counts = []
            for (k,mut_type) in enumerate(mut_types):
                raw_count_lines = [','.join( ['0']*max_seq_length )]*len(design_name_idx)
                for n in range( N ):
                    design_name = lines[16*n].strip('\n')
                    if design_name not in design_name_idx: print('ERROR COULD NOT FIND design_name! '+design_name)
                    assert( design_name in design_name_idx )
                    idx = design_name_idx[ design_name ]
                    cols = lines[16*n+1+k].strip('\n').split()
                    assert( cols[0] == mut_type)
                    pad_string = ',0'*(max_seq_length-len(cols[1].split(',')))
                    raw_count_lines[idx] = cols[1] + pad_string

                outfile_raw_counts = outdir + '%s.%s.txt.gz' % (primer_name,mut_type)
                outfiles_raw_counts.append(outfile_raw_counts)
                fid_raw_counts = gzip.open(outfile_raw_counts,'wt')
                for line in raw_count_lines: fid_raw_counts.write( line+'\n' )
                fid_raw_counts.close()
            print( 'Created %s for %d sequences (found %d designs)' % (','.join(outfiles_raw_counts),len(design_name_idx),N) )
        else:
            print( 'WARNING! Could not find %s' % infile )

time_end=time.time()

if len(wd) == 0: wd = './'
ubr_check_stats.check_stats( wd )

print( '\nTimings:')
print( '0_merge_pairs ' + time.strftime("%H:%M:%S",time.gmtime(time_merge_pairs-time_start) ) )
print( '1_ultraplex   ' + time.strftime("%H:%M:%S",time.gmtime(time_ultraplex1-time_merge_pairs) ) )
print( '2_bowtie2     ' + time.strftime("%H:%M:%S",time.gmtime(time_bowtie2-time_ultraplex1) ) )
if args.cmuts:
    print( '3_cmuts    ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_bowtie2) ) )
else:
    print( '3_rf_count    ' + time.strftime("%H:%M:%S",time.gmtime(time_rf_count-time_bowtie2) ) )
    print( '4_rctools     ' + time.strftime("%H:%M:%S",time.gmtime(time_rctools-time_rf_count) ) )
    print( 'Final merge   ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_rctools) ) )

print( '\nTotal time: ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_start) ) )

