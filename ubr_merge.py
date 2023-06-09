#!/usr/bin/env python3

import argparse
import os
import glob
import shutil

parser = argparse.ArgumentParser(
                    prog = 'ubr_merge.py',
                    description = 'Merge outputs from split ubr_run directories',
                    epilog = 'Merge .muts.txt and .coverage.txt output files from split ubr_run job directories.')

args = parser.add_argument( 'split_dir',default='UBR',help='name of directory with 001/,002/ job directories (default UBR)' )
args = parser.add_argument( '--merge_files' ,nargs='+',help='Filenames to find and merge')
args = parser.add_argument( '-s','--setup_slurm', action='store_true',help='Instead of running, create sbatch script' )
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
    # Do the merges
    for filename in merge_files:
        infiles = sorted(glob.glob('%s/*/%s' % (args.split_dir,filename) ))
        outdir = ''
        if len(infiles) == 0:
            infiles = sorted(glob.glob('%s/*/raw_counts/%s' % (args.split_dir,filename) ))
            outdir = 'raw_counts/'
            os.makedirs( outdir, exist_ok = True )

        counts = []
        for infile in infiles:
            lines = open(infile).readlines()
            for (i,line) in enumerate(lines):
                if len(counts)<=i: counts.append([])
                for (j,count) in enumerate(line.split()):
                    if len(counts[i])<=j: counts[i].append(0)
                    counts[i][j] += int(count)
        nseq = len(counts)

        # pad all counts lines to same length
        max_seq_length = max( map( len, counts ) )
        for i in range(nseq):
            for j in range(len(counts[i]),max_seq_length):
                counts[i].append(0)

        tot_counts = 0
        outfile = filename
        fid = open( outdir+outfile, 'w' )
        for line in counts:
            tot_counts += max(line)
            for (j,count) in enumerate(line):
                fid.write('%d' % count)
                if j == len(line)-1: fid.write('\n')
                else: fid.write(' ')
        print( 'Compiled %8d total counts for %6d sequences from %6d files into: %s' % (tot_counts,nseq,len(infiles),fid.name) )
        fid.close()
