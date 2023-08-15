#!/usr/bin/env python3

import argparse
import os
import shutil
import time
import glob

parser = argparse.ArgumentParser(
                    prog = 'ubr_split.py',
                    description = 'Get ready for ubr_run.py.',
                    epilog = 'Split FASTQ files and prepare independent job directories and command-lines.')

parser.add_argument('-s','--sequences_fasta', required=True)
parser.add_argument('-b','--primer_barcodes_fasta', required=True)
parser.add_argument('-1','--read1_fastq', required=True)
parser.add_argument('-2','--read2_fastq', required=True)
parser.add_argument('-n','--nsplits', default=0, type=int, help='number of separate partitions' )
parser.add_argument('-q','--sequences_per_partition', default=0, type=int, help='number of sequences in each partition. overrides -n/--nsplits.' )
parser.add_argument('-j','--jobs_per_slurm_node', default=24,type=int )
parser.add_argument('-ow','--overwrite',action = 'store_true')
parser.add_argument('-mp','--merge_pairs',action = 'store_true',help='Merge paired reads')

parser.add_argument('-nm','--no_mixed',action = 'store_true',help='No mixed reads in Bowtie2')
parser.add_argument('-sm','--score_min',help='minimum score for Bowtie2')
parser.add_argument('-mq','--map_quality',default=10,type=int,help='minimum Bowtie2 MAPQ to consider read')

parser.add_argument('-lc','--length_cutoff',action = 'store_true',help='Use length cutoff of 0.92 length for RNAframework')
parser.add_argument('-nlc','--no_length_cutoff',action = 'store_true')
parser.add_argument('-norc','--no_output_raw_counts',action = 'store_true',help='do not output raw counts from RNAframework')
parser.add_argument('-orc','--output_raw_counts',action = 'store_true')
parser.add_argument('-me','--max_edit_distance',default=0.0,type=float,help='max edit distance for RNAFramework (0.15)')

args = parser.parse_args()

if ( args.nsplits == 0 and args.sequences_per_partition == 0 ):
    print( "Must specify either -n/--nsplits or -q/--sequences_per_partition" )
    exit()
if args.no_length_cutoff:  print( '\n--no_length_cutoff is on by default now!\n' )
if args.output_raw_counts: print( '\n--output_raw_counts is on by default now!\n' )

time_start = time.time()

# Check for executable!
assert( shutil.which( 'seqkit' ) )
assert( shutil.which( 'ultraplex' ) )
assert( shutil.which( 'bowtie2-build' ) )
assert( shutil.which( 'bowtie2' ) )
assert( shutil.which( 'rf-count' ) )
assert( shutil.which( 'samtools' ) )

# Check for files
assert( os.path.isfile( args.sequences_fasta ) )
assert( os.path.isfile( args.primer_barcodes_fasta ) )
assert( os.path.isfile( args.read1_fastq ) )
assert( os.path.isfile( args.read2_fastq ) )
assert( args.read1_fastq != args.read2_fastq )

# Check FASTA
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
print( 'Read in %d sequences from %s.' % (len(sequences),args.sequences_fasta) )
(primer_barcodes,primer_names) = read_fasta( args.primer_barcodes_fasta )
print( 'Read in %d primer barcodes from %s.\n' % (len(primer_barcodes),args.primer_barcodes_fasta) )

def check_sequence(sequence):
    for c in sequence:
        if c not in 'ACGTU': return False
    return True

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

# Check if splits exist...
for i in range(1,nsplits+1):
    index_tag = '%03d' % i
    if i>999: index_tag = '%d' % i
    outdir = '%s/%s' % (split_dir,index_tag)
    f1 = '%s.part_%s%s' % ( os.path.basename( args.read1_fastq ).replace(fastq_tag,''), index_tag, fastq_tag )
    f2 = '%s.part_%s%s' % ( os.path.basename( args.read2_fastq ).replace(fastq_tag,''), index_tag, fastq_tag )
    if not os.path.isdir( outdir ) or not os.path.isfile( outdir+'/'+f1 ) or not os.path.isfile( outdir+'/'+f2 ):
        already_did_split = False
        break

if already_did_split and not args.overwrite:
    print( 'Skipping seqkit split2 into %s!' % split_dir )
else:
    if args.sequences_per_partition > 0:
        command = 'seqkit split2 -s %d -1 %s -2 %s -O %s --threads 12' % (args.sequences_per_partition,args.read1_fastq,args.read2_fastq,split_dir)
        print(command)
    else:
        command = 'seqkit split2 -p %d -1 %s -2 %s -O %s --threads 12' % (nsplits,args.read1_fastq,args.read2_fastq,split_dir)
        print(command)
    os.system( command )

time_seqkit = time.time()

# Create job directories and compile all commands
all_commands_file = 'all_commands.sh'
fid_all = open( all_commands_file, 'w' )

slurm_file_dir = 'slurm_files'
os.makedirs(slurm_file_dir, exist_ok=True )

slurm_file_count = 1
fid_slurm = open( '%s/run_slurm_%03d.sh' % (slurm_file_dir, slurm_file_count), 'w' )
sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_run\n#SBATCH --output=ubr_run.o%%j\n#SBATCH --error=ubr_run.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=48:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n#SBATCH --mem=%dG\n\n' % (args.jobs_per_slurm_node,2*args.jobs_per_slurm_node)
fid_slurm.write( sbatch_preface )
fid_sbatch_commands = open( 'sbatch_commands.sh', 'w')
ubr_run_sh_name = 'ubr_run.sh'

if nsplits == 1:
    # We're probably running with sequences_per_partition specified. Let's actually count splits.
    nsplits = len(glob.glob( '%s/*/%s.part_*%s' % ( split_dir, os.path.basename( args.read1_fastq ).replace(fastq_tag,''), fastq_tag ) ))
    if nsplits == 0:
        nsplits = len(glob.glob( '%s/%s.part_*%s' % ( split_dir,os.path.basename( args.read1_fastq ).replace(fastq_tag,''), fastq_tag ) ))
    assert( nsplits > 0 )

for i in range(1,nsplits+1):
    index_tag = '%03d' % i
    if i>999: index_tag = '%d' % i
    f1 = '%s.part_%s%s' % ( os.path.basename( args.read1_fastq ).replace(fastq_tag,''), index_tag, fastq_tag )
    f2 = '%s.part_%s%s' % ( os.path.basename( args.read2_fastq ).replace(fastq_tag,''), index_tag, fastq_tag )

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
    if args.merge_pairs:  extra_flags += ' --merge_pairs'

    fid = open( outdir + '/'+ubr_run_sh_name, 'w' )
    fid.write( 'ubr_run.py -s %s -b %s -1 %s -2 %s%s > ubr_run.out 2> ubr_run.err & \n' % ( os.path.basename(args.sequences_fasta), os.path.basename(args.primer_barcodes_fasta), os.path.basename(f1), os.path.basename(f2),extra_flags) )
    fid.close()

    command = 'ubr_run.py -s %s/%s -b %s/%s -1 %s/%s -2 %s/%s -O %s%s > %s/ubr_run.out 2> %s/ubr_run.err &' % ( outdir, args.sequences_fasta, outdir, args.primer_barcodes_fasta, outdir, f1, outdir, f2, outdir, extra_flags, outdir, outdir)
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

print( '\nTimings:')
print( 'seqkit split2 ' + time.strftime("%H:%M:%S",time.gmtime(time_seqkit-time_start) ) )
print( 'file mv/cp    ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_seqkit) ) )

print( '\nTotal time: ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_start) ) )

print( '\n\nAll %d commands can be run with:\n source %s' % (nsplits,fid_all.name) )
print( "\nOr you can go into each subdirectory of %s/ and run:\n source %s" % (split_dir, ubr_run_sh_name ) )
print( "\nOr to queue up %d slurm jobs on sherlock you can run:\n source %s\n" % (slurm_file_count-1,fid_sbatch_commands.name) )

fid_all.close()
fid_sbatch_commands.close()

