#!/usr/bin/env python3

import argparse
import os
import shutil
import time
import glob
from ubr_util import read_fasta,check_sequence

parser = argparse.ArgumentParser(
                    prog = 'ubr_split.py',
                    description = 'Get ready for ubr_run.py.',
                    epilog = 'Split FASTQ files and prepare independent job directories and command-lines.')

parser.add_argument('-s','--sequences_fasta', required=True, help='FASTA of RNA sequences')
parser.add_argument('-b','--primer_barcodes_fasta', required=True, help='FASTA of primer barcodes, first nucleotides of Read 1, prepended to cDNA')
parser.add_argument('-1','--read1_fastq', required=True, help='FASTQ (can be gzipped) of Illumina run')
parser.add_argument('-2','--read2_fastq', help='FASTQ (can be gzipped) of Read 2')
parser.add_argument('-n','--nsplits', default=0, type=int, help='number of separate partitions' )
parser.add_argument('-q','--sequences_per_partition', default=0, type=int, help='number of sequences in each partition. overrides -n/--nsplits.' )
parser.add_argument('-j','--jobs_per_slurm_node', default=24,type=int )
parser.add_argument('-ow','--overwrite',action = 'store_true', help='overwrite all previous files')
parser.add_argument('-nmp','--no_merge_pairs',action = 'store_true',help='do not merge paired end reads before Bowtie2' )
parser.add_argument('--cmuts',action = 'store_true',help='use cmuts instead of RNAFramework' )
parser.add_argument('--precomputed_bowtie_build_dir',default='',help=argparse.SUPPRESS) # precomputed bowtie-build directory, which otherwise takes forever to generate on the fly for >1M seqs
parser.add_argument('-f','--force',action='store_true',help='override some warnings' )

# Deprecated
parser.add_argument('--skip_gzip',action='store_true',help=argparse.SUPPRESS) # if FASTQ is not gzipped, leave it gzipped
parser.add_argument('--excise_barcode',default=0,type=int,help=argparse.SUPPRESS) # 'remove this many nucleotides from merged FASTQ and sequences' )
parser.add_argument('-nm','--no_mixed',action = 'store_true',help=argparse.SUPPRESS )# 'No mixed reads in Bowtie2')
parser.add_argument('-sm','--score_min',help=argparse.SUPPRESS )#'minimum score for Bowtie2')
parser.add_argument('-mq','--map_quality',default=10,type=int,help=argparse.SUPPRESS )#help='minimum Bowtie2 MAPQ to consider read')
parser.add_argument('-lc','--length_cutoff',action = 'store_true',help=argparse.SUPPRESS )#help='Use length cutoff of 0.92 length for RNAFramework')
parser.add_argument('-norc','--no_output_raw_counts',action = 'store_true',help=argparse.SUPPRESS )#help='do not output raw counts from RNAFramework')
parser.add_argument('-me','--max_edit_distance',default=0.0,type=float,help=argparse.SUPPRESS )#help='max edit distance for RNAFramework (0.15)')
parser.add_argument('-mpp','--merge_pairs_pear',action = 'store_true',help=argparse.SUPPRESS)
parser.add_argument('--force_merge_pairs',action = 'store_true',help=argparse.SUPPRESS) # force merge pairs (don't bother to check for overlap)
parser.add_argument('--cutadapt',action = 'store_true',help=argparse.SUPPRESS) # force cutadapt trimming of Read2 side for pre-demuxed Ultima
parser.add_argument('--use_tmp_dir',action = 'store_true',help=argparse.SUPPRESS) # For cmuts, run job in /tmp/ to try to reduce disk i/o
parser.add_argument('--ultima',action='store_true',help=argparse.SUPPRESS) # recognize Ultima adapter in ultraplex

args = parser.parse_args()

if ( args.nsplits == 0 and args.sequences_per_partition == 0 ):
    print( "Must specify either -n/--nsplits or -q/--sequences_per_partition" )
    exit()
assert( not( args.no_merge_pairs and args.merge_pairs_pear ) )

time_start = time.time()

# Check for executable!
assert( shutil.which( 'seqkit' ) )
assert( shutil.which( 'ultraplex' ) )
assert( shutil.which( 'bowtie2-build' ) )
assert( shutil.which( 'bowtie2' ) )
assert( shutil.which( 'rf-count' ) )
assert( shutil.which( 'samtools' ) )
assert( shutil.which( 'bbmerge.sh' ) )
assert( shutil.which( 'java' ) )
if args.merge_pairs_pear: assert( shutil.which('pear') )

# Check for files
assert( os.path.isfile( args.sequences_fasta ) )
assert( os.path.isfile( args.primer_barcodes_fasta ) )
assert( os.path.isfile( args.read1_fastq ) )
assert( args.read2_fastq == None or os.path.isfile( args.read2_fastq ) )
assert( args.read1_fastq != args.read2_fastq )

(sequences,headers) = read_fasta( args.sequences_fasta )
print( 'Read in %d sequences from %s.' % (len(sequences),args.sequences_fasta) )
(primer_barcodes,primer_names) = read_fasta( args.primer_barcodes_fasta )
print( 'Read in %d primer barcodes from %s.\n' % (len(primer_barcodes),args.primer_barcodes_fasta) )
for primer_name in primer_names: assert( primer_name.find(' ')==-1 )

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

# Do the split
split_dir = 'UBR'

already_did_split = True
nsplits = args.nsplits
if args.sequences_per_partition > 0: nsplits = 1

# Figure out the fastq tag
fastq_tag = '.fastq.gz'
if args.read1_fastq.find( fastq_tag ) < 0: fastq_tag = '.fastq'
if args.read1_fastq.find( fastq_tag ) < 0: fastq_tag = '.fq.gz'
if args.read1_fastq.find( fastq_tag ) < 0: fastq_tag = '.fq'
assert( args.read1_fastq.find( fastq_tag ) >= 0 )

fastq_tag_out = fastq_tag
if not args.skip_gzip and fastq_tag[-3:] != '.gz': fastq_tag_out += '.gz'

# Check if splits exist...
for i in range(1,nsplits+1):
    index_tag = '%03d' % i
    if i>999: index_tag = '%d' % i
    outdir = '%s/%s' % (split_dir,index_tag)
    f1 = '%s.part_%s%s' % ( os.path.basename( args.read1_fastq ).replace(fastq_tag,''), index_tag, fastq_tag_out )
    if args.read2_fastq: f2 = '%s.part_%s%s' % ( os.path.basename( args.read2_fastq ).replace(fastq_tag,''), index_tag, fastq_tag_out )
    else: f2 = f1
    if not os.path.isdir( outdir ) or not os.path.isfile( outdir+'/'+f1 ) or not os.path.isfile( outdir+'/'+f2 ):
        already_did_split = False
        break

if already_did_split and not args.overwrite:
    print( 'Skipping seqkit split2 into %s!' % split_dir )
else:
    fastq2_tag = ''
    if args.read2_fastq: fastq2_tag = ' -2 %s' % args.read2_fastq
    if args.sequences_per_partition > 0:
        command = 'seqkit split2 -s %d -1 %s%s -O %s --threads 12' % (args.sequences_per_partition,args.read1_fastq,fastq2_tag,split_dir)
        print(command)
    else:
        command = 'seqkit split2 -p %d -1 %s%s -O %s --threads 12' % (nsplits,args.read1_fastq,fastq2_tag,split_dir)
        print(command)
    if not args.skip_gzip: command += ' -e .gz'
    os.system( command )

time_seqkit = time.time()

# Create job directories and compile all commands
all_commands_file = 'all_commands.sh'
fid_all = open( all_commands_file, 'w' )

slurm_file_dir = 'slurm_files'
os.makedirs(slurm_file_dir, exist_ok=True )

slurm_file_count = 1
fid_slurm = open( '%s/run_slurm_%03d.sh' % (slurm_file_dir, slurm_file_count), 'w' )
sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_run\n#SBATCH --output=ubr_run.o%%j\n#SBATCH --error=ubr_run.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=48:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n#SBATCH --mem=%dG\n\n' % (args.jobs_per_slurm_node,8*args.jobs_per_slurm_node)
fid_slurm.write( sbatch_preface )
fid_sbatch_commands = open( 'sbatch_commands.sh', 'w')
ubr_run_sh_name = 'ubr_run.sh'

if nsplits == 1:
    # We're probably running with sequences_per_partition specified. Let's actually count splits.
    nsplits = len(glob.glob( '%s/*/%s.part_*%s' % ( split_dir, os.path.basename( args.read1_fastq ).replace(fastq_tag,''), fastq_tag_out ) ))
    nsplits -= len(glob.glob( '%s/*/%s.part*MERGED*%s' % ( split_dir,os.path.basename( args.read1_fastq ).replace(fastq_tag,''), fastq_tag_out ) ))
    if nsplits == 0:
        nsplits = len(glob.glob( '%s/%s.part*%s' % ( split_dir,os.path.basename( args.read1_fastq ).replace(fastq_tag,''), fastq_tag_out ) ))
    print('Found: %d directories.' % nsplits)
    assert( nsplits > 0 )

for i in range(1,nsplits+1):
    index_tag = '%03d' % i
    if i>999: index_tag = '%d' % i
    f1 = '%s.part_%s%s' % ( os.path.basename( args.read1_fastq ).replace(fastq_tag,''), index_tag, fastq_tag_out )
    if args.read2_fastq: f2 = '%s.part_%s%s' % ( os.path.basename( args.read2_fastq ).replace(fastq_tag,''), index_tag, fastq_tag_out )
    else: f2 = f1

    outdir = '%s/%s' % (split_dir,index_tag)
    os.makedirs(outdir, exist_ok=True )

    if not os.path.isfile( outdir + '/' + f1 ):
        assert( os.path.isfile( split_dir + '/' + f1 ) )
        shutil.move( split_dir + '/' +  f1, outdir + '/' + f1 )

    if not os.path.isfile( outdir + '/' + f2 ):
        assert( os.path.isfile( split_dir + '/' + f2 ) )
        shutil.move( split_dir + '/' +  f2, outdir + '/' + f2 )

    if not os.path.isfile( outdir + '/' + os.path.basename(args.sequences_fasta) ):
        shutil.copy( args.sequences_fasta, outdir + '/' + os.path.basename(args.sequences_fasta) )

    if not os.path.isfile( outdir + '/' + os.path.basename(args.primer_barcodes_fasta) ):
        shutil.copy( args.primer_barcodes_fasta, outdir + '/' + os.path.basename(args.primer_barcodes_fasta) )

    assert( os.path.isfile(outdir+'/'+os.path.basename(args.sequences_fasta)) )
    assert( os.path.isfile(outdir+'/'+os.path.basename(args.sequences_fasta)) )
    assert( os.path.isfile(outdir+'/'+f1) )
    assert( os.path.isfile(outdir+'/'+f2) )

    extra_flags = ''
    if args.no_output_raw_counts: extra_flags += ' --no_output_raw_counts'
    if args.length_cutoff:  extra_flags += ' --length_cutoff'
    if args.no_mixed:  extra_flags += ' --no_mixed'
    if args.map_quality != 10:  extra_flags += ' --map_quality %d' % args.map_quality
    if args.max_edit_distance > 0:  extra_flags += ' --max_edit_distance %f' % args.max_edit_distance
    if args.score_min != None:  extra_flags += ' --score_min %s' % args.score_min
    if args.no_merge_pairs:  extra_flags += ' --no_merge_pairs'
    if args.merge_pairs_pear:  extra_flags += ' --merge_pairs_pear'
    if args.force_merge_pairs:  extra_flags += ' --force_merge_pairs'
    if args.ultima:  extra_flags += ' --ultima'
    if args.cmuts:  extra_flags += ' --cmuts'
    if args.cutadapt:  extra_flags += ' --cutadapt'
    if args.use_tmp_dir:  extra_flags += ' --use_tmp_dir'
    if len(args.precomputed_bowtie_build_dir)>0: extra_flags += ' --precomputed_bowtie_build_dir %s' % args.precomputed_bowtie_build_dir
    if args.excise_barcode>0:  extra_flags += ' --excise_barcode %d' % args.excise_barcode

    fid = open( outdir + '/'+ubr_run_sh_name, 'w' )
    fastq2_tag = ''
    if args.read2_fastq: fastq2_tag = ' -2 %s' % os.path.basename(f2)
    fid.write( 'ubr_run.py -s %s -b %s -1 %s%s%s > ubr_run.out 2> ubr_run.err & \n' % ( os.path.basename(args.sequences_fasta), os.path.basename(args.primer_barcodes_fasta), os.path.basename(f1),fastq2_tag,extra_flags) )
    fid.close()

    fastq2_tag = ''
    if args.read2_fastq: fastq2_tag = ' -2 %s/%s' % (outdir,f2)
    command = 'ubr_run.py -s %s/%s -b %s/%s -1 %s/%s%s -O %s%s > %s/ubr_run.out 2> %s/ubr_run.err &' % ( outdir, args.sequences_fasta, outdir, args.primer_barcodes_fasta, outdir, f1, fastq2_tag, outdir, extra_flags, outdir, outdir)
    fid_all.write( command + '\n' )

    fid_slurm.write( command +'\n' )
    if (i % args.jobs_per_slurm_node == 0 or i == nsplits ):
        fid_sbatch_commands.write('sbatch %s\n' % fid_slurm.name )
        fid_slurm.write('\nwait\necho "DONE"\n')
        fid_slurm.close()
        slurm_file_count += 1
        if i < nsplits:
            fid_slurm = open( '%s/run_slurm_%03d.sh' % (slurm_file_dir, slurm_file_count), 'w' )
            fid_slurm.write( sbatch_preface )

time_end = time.time()

# excise barcodes [testing if barcodes are important]
if args.excise_barcode > 0:
    new_sequence_file = args.sequences_fasta.replace('.fa','.NO_BARCODES.fa')
    fid = open( new_sequence_file, 'w' )
    for (i,(sequence,header)) in enumerate(zip(sequences,headers)):
        sequence = sequence[:(-args.excise_barcode)]
        fid.write('>%s\n%s\n' % (header,sequence))
        sequences[i] = sequence
    print('\nEXCISE_BARCODE: Outputted %d sequences with last %d nucleotides excised in %s\n' % (len(sequences),args.excise_barcode,new_sequence_file))

print( '\nTimings:')
print( 'seqkit split2 ' + time.strftime("%H:%M:%S",time.gmtime(time_seqkit-time_start) ) )
print( 'file mv/cp    ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_seqkit) ) )

print( '\nTotal time: ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_start) ) )

print( '\n\nAll %d commands can be run with:\n source %s' % (nsplits,fid_all.name) )
print( "\nOr you can go into each subdirectory of %s/ and run:\n source %s" % (split_dir, ubr_run_sh_name ) )
print( "\nOr to queue up %d slurm jobs on sherlock you can run:\n source %s\n" % (slurm_file_count-1,fid_sbatch_commands.name) )

fid_all.close()
fid_sbatch_commands.close()

