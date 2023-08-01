#!/usr/bin/env python3

import argparse
import os
import glob
import shutil
import time
import pandas as pd

parser = argparse.ArgumentParser(
                    prog = 'ubr_merge.py',
                    description = 'Merge outputs from split ubr_run directories',
                    epilog = 'Merge .muts.txt and .coverage.txt output files from split ubr_run job directories.')

args = parser.add_argument( 'split_dir',default='UBR',help='name of directory with 001/,002/ job directories (default UBR)' )
args = parser.add_argument( '--merge_files' ,nargs='+',help='Filenames to find and merge')
args = parser.add_argument( '-s','--setup_slurm', action='store_true',help='Instead of running, create sbatch script' )
args = parser.add_argument('-j','--jobs_per_slurm_node', default=16,type=int )
args = parser.parse_args()

merge_files = args.merge_files

if merge_files == None:
    unique_files = []
    for tag in ['coverage','muts']:
        # Find all the relevant files
        globfiles = sorted(glob.glob('%s/*/*.%s.txt' % (args.split_dir,tag)  ))
        # Get unique names
        unique_files = unique_files + sorted(list(set(map( os.path.basename, globfiles ))))

    mut_types = ['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del']
    for tag in mut_types:
        # Find all the relevant files
        globfiles = sorted(glob.glob('%s/*/raw_counts/*.%s.txt' % (args.split_dir,tag)  ))
        # Get unique names
        unique_files = unique_files + sorted(list(set(map( os.path.basename, globfiles ))))

    merge_files = unique_files

print( merge_files )

if args.setup_slurm:
    if len(merge_files) < args.jobs_per_slurm_node:
        sbatch_file = 'run_ubr_merge.sh'
        fid_slurm = open( sbatch_file, 'w' )
        sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_merge\n#SBATCH --output=ubr_merge.o%%j\n#SBATCH --error=ubr_merge.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=8:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n\n' % len(merge_files)
        fid_slurm.write( sbatch_preface )
        for merge_file in merge_files:
            command = 'ubr_merge.py %s --merge_files %s &' % ( args.split_dir, merge_file )
            fid_slurm.write( command +'\n' )
        fid_slurm.write('\nwait\n')
        fid_slurm.close()
        print( '\nCreated %s with %d commands. Run:\n sbatch %s\n\n' % (sbatch_file,len(merge_files),sbatch_file) )
    else:
        slurm_file_dir = 'slurm_files'
        os.makedirs(slurm_file_dir, exist_ok=True )

        slurm_file_count = 1
        fid_slurm = open( '%s/run_ubr_merge_%03d.sh' % (slurm_file_dir, slurm_file_count), 'w' )
        sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_merge\n#SBATCH --output=ubr_merge.o%%j\n#SBATCH --error=ubr_merge.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=2:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n\n' % args.jobs_per_slurm_node
        fid_slurm.write( sbatch_preface )
        fid_sbatch_commands = open( 'sbatch_merge_commands.sh', 'w')

        for (i,merge_file) in enumerate(merge_files):
            command = 'ubr_merge.py %s --merge_files %s &' % ( args.split_dir, merge_file )
            fid_slurm.write( command +'\n' )
            if ( (i+1) % args.jobs_per_slurm_node == 0 ) or i == len(merge_files)-1:
                fid_sbatch_commands.write('sbatch %s\n' % fid_slurm.name )
                fid_slurm.write('\nwait\necho "DONE"\n')
                fid_slurm.close()
                if i < len(merge_files)-1:
                    slurm_file_count += 1
                    fid_slurm = open( '%s/run_ubr_merge_%03d.sh' % (slurm_file_dir, slurm_file_count), 'w' )
                    fid_slurm.write( sbatch_preface )
        fid_sbatch_commands.close()
        print( '\nCreated %s slurm files containing %d commands. Run:\n source %s\n\n' % (slurm_file_count,len(merge_files),fid_sbatch_commands.name) )
    exit(0)

# Do the merges
time_start = time.time()
time_readin = 0
time_output = 0
for filename in merge_files:
    time_startfile = time.time()

    infiles = sorted(glob.glob('%s/*/%s' % (args.split_dir,filename) ))
    outdir = ''
    if len(infiles) == 0:
        infiles = sorted(glob.glob('%s/*/raw_counts/%s' % (args.split_dir,filename) ))
        outdir = 'raw_counts/'
        os.makedirs( outdir, exist_ok = True )

    counts = []
    df = None
    for infile in infiles:
        df_infile = pd.read_table(infile,sep=' ',header=None)
        if df is None:
            df = df_infile
        else:
            df = df.add(df_infile,fill_value = 0)
    nseq = len(df)

    time_after_infile = time.time()

    time_readin += time_after_infile - time_startfile

    outfile = outdir+filename
    tot_counts = df.max(axis=1).sum()
    df.to_csv(outfile,sep=' ',header=None,index=None)

    print( 'Compiled %8d total counts for %6d sequences from %6d files into: %s' % (tot_counts,nseq,len(infiles),outfile) )

    time_output += (time.time()-time_after_infile)

time_end = time.time()

print( '\nTimings:')
print( 'Read in data: ' + time.strftime("%H:%M:%S",time.gmtime(time_readin) ) )
print( 'Output  data: ' + time.strftime("%H:%M:%S",time.gmtime(time_output) ) )

print( '\nTotal time: ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_start) ) )

