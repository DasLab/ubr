#!/usr/bin/env python3

import argparse
import os
import glob
import shutil
import time
import gzip
import pandas as pd

parser = argparse.ArgumentParser(
                    prog = 'ubr_subdivide.py',
                    description = 'Separate out sequence & sublibraries from UBR',
                    epilog = 'Subdivide .muts.txt and .coverage.txt output files from split ubr job directory.')

parser.add_argument('--sequences_fasta', required=True, help='FASTA of RNA sequences, headers should include tags like sublibrary:my_library1.')
parser.add_argument( '--subdivide_files' ,nargs='+',help='Filenames to find and subdivide')
parser.add_argument( 'ubr_dir',default='./',help='name of directory with files *.mut.txt, raw_counts/, etc.' )
parser.add_argument( '-s','--setup_slurm', action='store_true',help='Instead of running, create sbatch script' )
parser.add_argument('-j','--jobs_per_slurm_node', default=16,type=int )

args = parser.parse_args()
subdivide_files = args.subdivide_files
ubr_dir = args.ubr_dir

if subdivide_files == None:
    # TODO unify with same code block in ubr_subdivide.py
    globfiles = []
    for tag in ['coverage','muts']:
        # Find all the relevant files
        globfiles += sorted(glob.glob('%s/*.%s.txt' % (ubr_dir,tag)  ))
        globfiles += sorted(glob.glob('%s/*.%s.txt.gz' % (ubr_dir,tag)  ))

    mut_types = ['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del']
    for tag in mut_types:
        # Find all the relevant files
        globfiles += sorted(glob.glob('%s/raw_counts/*.%s.txt' % (ubr_dir,tag)  ))
        globfiles += sorted(glob.glob('%s/raw_counts/*.%s.txt.gz' % (ubr_dir,tag)  ))

    globfiles = [x.replace('.gz','') for x in globfiles]
    unique_files = sorted(list(set(map( os.path.basename, globfiles ))))
    subdivide_files = unique_files

print('\nWill subdivide:')
for subdivide_file in subdivide_files:
    print(subdivide_file)

if args.setup_slurm:
    # TODO: unify with same code block in ubr_merge.py
    if len(subdivide_files) < args.jobs_per_slurm_node:
        sbatch_file = 'run_ubr_subdivide.sh'
        fid_slurm = open( sbatch_file, 'w' )
        sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_subdivide\n#SBATCH --output=ubr_subdivide.o%%j\n#SBATCH --error=ubr_subdivide.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=8:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n\n' % len(subdivide_files)
        fid_slurm.write( sbatch_preface )
        for subdivide_file in subdivide_files:
            command = 'ubr_subdivide.py %s --subdivide_files %s --sequences_fasta %s &' % ( ubr_dir, subdivide_file, args.sequences_fasta )
            fid_slurm.write( command +'\n' )
        fid_slurm.write('\nwait\n')
        fid_slurm.close()
        print( '\nCreated %s with %d commands. Run:\n sbatch %s\n\n' % (sbatch_file,len(subdivide_files),sbatch_file) )
    else:
        slurm_file_dir = 'slurm_files'
        os.makedirs(slurm_file_dir, exist_ok=True )

        slurm_file_count = 1
        fid_slurm = open( '%s/run_ubr_subdivide_%03d.sh' % (slurm_file_dir, slurm_file_count), 'w' )
        sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_subdivide\n#SBATCH --output=ubr_subdivide.o%%j\n#SBATCH --error=ubr_subdivide.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=2:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n\n' % args.jobs_per_slurm_node
        fid_slurm.write( sbatch_preface )
        fid_sbatch_commands = open( 'sbatch_subdivide_commands.sh', 'w')

        for (i,subdivide_file) in enumerate(subdivide_files):
            command = 'ubr_subdivide.py %s --subdivide_files %s  --sequences_fasta %s &' % ( ubr_dir, subdivide_file, args.sequences_fasta )
            fid_slurm.write( command +'\n' )
            if ( (i+1) % args.jobs_per_slurm_node == 0 ) or i == len(subdivide_files)-1:
                fid_sbatch_commands.write('sbatch %s\n' % fid_slurm.name )
                fid_slurm.write('\nwait\necho "DONE"\n')
                fid_slurm.close()
                if i < len(subdivide_files)-1:
                    slurm_file_count += 1
                    fid_slurm = open( '%s/run_ubr_subdivide_%03d.sh' % (slurm_file_dir, slurm_file_count), 'w' )
                    fid_slurm.write( sbatch_preface )
        fid_sbatch_commands.close()
        print( '\nCreated %s slurm files containing %d commands. Run:\n source %s\n\n' % (slurm_file_count,len(subdivide_files),fid_sbatch_commands.name) )
    exit(0)


# Check FASTA
# TODO: use biopython or at least shared util.py
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
print( '\nRead in %d sequences from %s.' % (len(sequences),args.sequences_fasta) )

time_start = time.time()

sublibrary_idx = {}
sublibrary_tag = 'sublibrary:'
for (i,header) in enumerate(headers):
    sublibrary_tag_idx = header.find(sublibrary_tag)
    if sublibrary_tag_idx < 0:
        sublibrary = 'unassigned'
    else:
        sublibrary = header[sublibrary_tag_idx+len(sublibrary_tag):].split(', \t')[0]
    if sublibrary not in sublibrary_idx:
        sublibrary_idx[sublibrary] = []
    sublibrary_idx[sublibrary].append(i)

if len(sublibrary_idx.keys()) == 1 and 'unassigned' in sublibrary_idx:
    print( '\nDid not find any sublibrary: tags in %s. No subdividing to do! Exiting...\n')
    exit()

for sublibrary in sublibrary_idx.keys():
    print( 'Sublibrary: %30s with %6d sequences' % (sublibrary,len(sublibrary_idx[sublibrary])) )

time_sequence_readin = time.time()

# Do the subdivides of sequence files
def write_fasta( fasta_file, sequences, headers, idx ):
    fid = open(fasta_file,'w')
    for i in idx:
        sequence = sequences[i]
        header = headers[i]
        fid.write('>%s\n%s\n' % (header,sequence) )
    fid.close()
    print('Outputted %d sequences to %s' % (len(idx),fasta_file))

sublibrary_dir = 'SUBLIBRARIES'
print('\nChecking sequence files...')
for sublibrary in sublibrary_idx:
    dirname = '%s/%s' % (sublibrary_dir,sublibrary)
    os.makedirs( dirname, exist_ok = True )
    #outfile = dirname + '/' + os.path.basename(args.sequences_fasta).replace('.fa', '_%s.fa' % sublibrary)
    outfile = '%s/%s.fa' % (dirname,sublibrary)
    if os.path.isfile(outfile): continue
    write_fasta( outfile, sequences, headers, sublibrary_idx[sublibrary] )

# Do the subdivides of data files
time_readin = 0
time_output = 0
for filename in subdivide_files:
    time_startfile = time.time()

    outdir = ''
    infile = '%s/%s%s' % (ubr_dir,outdir,filename)
    if not os.path.isfile(infile):
        infile = '%s/%s%s.gz' % (ubr_dir,outdir,filename)
    if not os.path.isfile(infile):
        outdir = 'raw_counts/'
        infile = '%s/%s%s' % (ubr_dir,outdir,filename)
        if not os.path.isfile(infile):
            infile = '%s/%s%s.gz' % (ubr_dir,outdir,filename)

    # Check if we need to make any files:
    need_to_make_outfile = 0
    for sublibrary in sublibrary_idx:  # list(sublibrary_idx.keys())[:1]:
        dirname = '%s/%s/%s' % (sublibrary_dir,sublibrary,outdir)
        if not os.path.isdir(dirname): os.makedirs( dirname )
        outfile= '%s%s' % (dirname,filename)
        if os.path.isfile(outfile): continue
        need_to_make_outfile = 1
    if not need_to_make_outfile: continue
    print()

    # need to figure out separator, b/c pandas doesn't do a good job.
    sepchar = ','
    if infile.find('.gz')>-1:
        fid = gzip.open(infile,'rt')
    else:
        fid = open(infile)
    if fid.readline().find(sepchar) == -1: sepchar = ' '
    fid.close()

    # OK read it in!
    df = pd.read_table(infile,sep=sepchar,header=None)

    time_after_infile = time.time()

    time_readin += time_after_infile - time_startfile

    for sublibrary in sublibrary_idx:  # list(sublibrary_idx.keys())[:1]:
        dirname = '%s/%s/%s' % (sublibrary_dir,sublibrary,outdir)
        os.makedirs( dirname, exist_ok = True )
        outfile= '%s%s.gz' % (dirname,filename)
        if os.path.isfile(outfile): continue
        print( 'Creating... %s' % outfile )
        df_sub = df[df.index.isin( sublibrary_idx[sublibrary] )]
        df_sub.to_csv(outfile,sep=',',header=None,index=None)
    time_output += (time.time()-time_after_infile)

time_end = time.time()

print( '\nTimings:')
print( 'Read in data: ' + time.strftime("%H:%M:%S",time.gmtime(time_sequence_readin-time_start) ) )
print( 'Read in data: ' + time.strftime("%H:%M:%S",time.gmtime(time_readin) ) )
print( 'Output  data: ' + time.strftime("%H:%M:%S",time.gmtime(time_output) ) )

print( '\nTotal time: ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_start) ) )

