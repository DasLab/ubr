#!/usr/bin/env python3

import argparse
import os
import os.path
import glob
import shutil
import time
import gzip
import pandas as pd

parser = argparse.ArgumentParser(
                    prog = 'ubr_merge.py',
                    description = 'Merge outputs from split ubr_run directories',
                    epilog = 'Merge .muts.txt and .coverage.txt output files from split ubr_run job directories.')

args = parser.add_argument( 'split_dir',default='UBR/',help='name of directory with 001/,002/ job subdirectories; if multiple directories, just look directly inside directories.',nargs='+' )
args = parser.add_argument( '--merge_files' ,nargs='+',help='Filenames to find and merge')
args = parser.add_argument( '-s','--setup_slurm', action='store_true',help='Instead of running, create sbatch script' )
args = parser.add_argument('-j','--jobs_per_slurm_node', default=16,type=int )
args = parser.add_argument( '--no_overwrite', action='store_true',help='Do not merge if merged files already exist' )
args = parser.parse_args()

split_dirs = args.split_dir
for (i,split_dir) in enumerate(split_dirs):
    if split_dir[:-1] != '/': split_dirs[i]=split_dir+'/'
if len(split_dirs) == 1: split_dirs[0] = split_dirs[0]+'*/'
print("Merging from the following directories: ")
for split_dir in split_dirs: print(' ',split_dir)

merge_files = args.merge_files

if merge_files == None:
    unique_files = []
    for split_dir in split_dirs:
        subdir = ''
        for tag in ['coverage','muts']:
            # Find all the relevant files
            globfiles = sorted(glob.glob('%s/*.%s.txt' % (split_dir,tag)  ))
            globfiles += sorted(glob.glob('%s/*.%s.txt.gz' % (split_dir,tag)  ))
            # Get unique names
            unique_files += sorted(list(set(map( os.path.basename, globfiles ))))

        mut_types = ['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del']
        for tag in mut_types:
            # Find all the relevant files
            globfiles = sorted(glob.glob('%s/raw_counts/*.%s.txt' % (split_dir,tag)  ))
            globfiles += sorted(glob.glob('%s/raw_counts/*.%s.txt.gz' % (split_dir,tag)  ))
            # Get unique names
            unique_files += sorted(list(set(map( os.path.basename, globfiles ))))

    merge_files = sorted(list(set(unique_files)))

if args.no_overwrite:
    merge_files = [x for x in merge_files if not os.path.isfile(x) and not os.path.isfile('raw_counts/%s' % x) ]

print('\nWill merge:')
for merge_file in merge_files:
    print(merge_file)

if args.setup_slurm:
    if len(merge_files) < args.jobs_per_slurm_node:
        sbatch_file = 'run_ubr_merge.sh'
        fid_slurm = open( sbatch_file, 'w' )
        sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_merge\n#SBATCH --output=ubr_merge.o%%j\n#SBATCH --error=ubr_merge.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=8:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n#SBATCH --mem=%dG\n\n' % (len(merge_files),4*len(merge_files))
        fid_slurm.write( sbatch_preface )
        for merge_file in merge_files:
            command = 'ubr_merge.py %s --merge_files %s &' % ( ' '.join(args.split_dir), merge_file )
            fid_slurm.write( command +'\n' )
        fid_slurm.write('\nwait\n')
        fid_slurm.close()
        print( '\nCreated %s with %d commands. Run:\n sbatch %s\n\n' % (sbatch_file,len(merge_files),sbatch_file) )
    else:
        slurm_file_dir = 'slurm_files'
        os.makedirs(slurm_file_dir, exist_ok=True )

        slurm_file_count = 1
        fid_slurm = open( '%s/run_ubr_merge_%03d.sh' % (slurm_file_dir, slurm_file_count), 'w' )
        sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_merge\n#SBATCH --output=ubr_merge.o%%j\n#SBATCH --error=ubr_merge.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=8:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n#SBATCH --mem=%dG\n\n' % (args.jobs_per_slurm_node,4*args.jobs_per_slurm_node)
        fid_slurm.write( sbatch_preface )
        fid_sbatch_commands = open( 'sbatch_merge_commands.sh', 'w')

        for (i,merge_file) in enumerate(merge_files):
            command = 'ubr_merge.py %s --merge_files %s &' % ( ' '.join(args.split_dir), merge_file )
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

    infiles = []
    for split_dir in split_dirs:
        infiles += sorted(glob.glob('%s/%s' % (split_dir,filename) ))
        infiles += sorted(glob.glob('%s/%s.gz' % (split_dir,filename) ))
        outdir = ''

    if len(infiles) == 0:
        for split_dir in split_dirs:
            infiles += sorted(glob.glob('%s/raw_counts/%s' % (split_dir,filename) ))
            infiles += sorted(glob.glob('%s/raw_counts/%s.gz' % (split_dir,filename) ))
            outdir = 'raw_counts/'
            os.makedirs( outdir, exist_ok = True )

    counts = []
    df = None
    df_init = False
    numfiles = 0
    for (i,infile) in enumerate(infiles):

        # need to figure out separator
        sepchar = ','
        if infile.find('.gz')>-1:
            fid = gzip.open(infile,'rt')
        else:
            fid = open(infile)
        if fid.readline().find(sepchar) == -1: sepchar = ' '
        fid.close()

        # OK read it in!
        try:
            # uint32 needed to prevent memory overflow. Note maximum is 4Gb, so this will soon be issue, and need to shift to streaming.
            df_infile = pd.read_table(infile,sep=sepchar,header=None,dtype='uint32')
            df_infile_read_correctly = True
        except pd.errors.EmptyDataError:
            print('Note: %s was empty. Skipping.' % infile )
            df_infile_read_correctly = False
            continue # will skip the rest of the block and move to next file
        if not df_init and df_infile_read_correctly:
            df = df_infile
            df_init = True
        elif df_infile_read_correctly:
            df = df.add(df_infile,fill_value = 0)
        numfiles += 1

    time_after_infile = time.time()

    time_readin += time_after_infile - time_startfile

    if len(infiles) == 0 or not df_init:
        print( 'Did not find any files to merge: %s' % filename )
    else:
        outfile = outdir+filename
        df.to_csv(outfile,sep=',',header=None,index=None)

        nseq = len(df)
        tot_counts = df.max(axis=1).sum()
        print( 'Compiled %8d total counts for %6d sequences from %6d of %6d files into: %s' % (tot_counts,nseq,numfiles,len(infiles),outfile) )

    time_output += (time.time()-time_after_infile)

time_end = time.time()

print( '\nTimings:')
print( 'Read in data: ' + time.strftime("%H:%M:%S",time.gmtime(time_readin) ) )
print( 'Output  data: ' + time.strftime("%H:%M:%S",time.gmtime(time_output) ) )

print( '\nTotal time: ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_start) ) )

