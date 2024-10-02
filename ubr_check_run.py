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
                    prog = 'ubr_check_run.py',
                    description = 'Check outputs from UBR',
                    epilog = 'Look for non-empty error or empty out files. Enable clean and rerun of those.')


args = parser.add_argument( '--clean', action='store_true',help='Get ready to clean out problem directories and prepare job submission file to rerun.' )
args = parser.parse_args()

ubr_dir = 'UBR'
assert(os.path.isdir(ubr_dir))

dirs = os.popen('ls -d %s/*' % ubr_dir).readlines()
numdirs = len(dirs)
no_run_dirs = []
for subdir in dirs:
    #lines = os.popen('ls %s/ubr_run.out' % subdir.strip() ).readlines()
    if not os.path.isfile( '%s/ubr_run.out' % subdir.strip() ): no_run_dirs.append(subdir.strip())

lines = os.popen('ls -alS %s/*/ubr_run.err' % ubr_dir).readlines()
err_dirs = []
for line in lines: 
    cols = line.split()
    if int(cols[4])>0 and len(cols)>8: err_dirs.append( cols[8] )

lines = os.popen('ls -alS %s/*/ubr_run.out'  % ubr_dir).readlines()
empty_out_dirs = []
for line in lines: 
    cols = line.split()
    if int(cols[4])==0 and len(cols)>8: empty_out_dirs.append( cols[8] )

lines = os.popen('ls -alS %s/*/ubr_run.out' % ubr_dir).readlines()
count = 0
num_empty_outs = 0
for line in lines:
    if int(line.split()[4])==0: num_empty_outs += 1

problem_dirs = ['/'.join(x.split('/')[0:2]) for x in list(set(no_run_dirs + err_dirs + empty_out_dirs)) ]

print('Reviewing run in UBR/')
print('Number of subdirectories                 : %5d' % numdirs)
print('Number of subdirs with jobs run          : %5d' % numdirs)
print()
print('Number of jobs without jobs run          : %5d' % len(no_run_dirs))
print('Number of jobs with non-empty ubr_run.err: %5d' % len(err_dirs))
print('Number of jobs with     empty ubr_run.out: %5d' % len(empty_out_dirs))
print('Number of problem subdirs:               : %5d' % len(problem_dirs) )

run_slurm_lines = []
for problem_dir in problem_dirs:       
    # look for slurm files that need to be rerun
    grepline = os.popen('grep %s/ slurm_files/run_slurm_*.sh' % problem_dir).readlines() 
    run_slurm_line = grepline[0].split(':')[0]
    run_slurm_lines.append( run_slurm_line )
run_slurm_lines = list( set( run_slurm_lines ) )
print( '\nNumber of slurm jobs that need to be re-run: %5d' % len(run_slurm_lines) )

if args.clean:
   for problem_dir in problem_dirs:
       command = 'rm -rf %s/?_*  %s/ubr_run.[eo]* %s/*MERGE*' % (problem_dir,problem_dir,problem_dir)
       print(command)
       os.system(command)
   rerun_commands_file = 'rerun_sbatch_commands.sh'
   fid = open( rerun_commands_file, 'w' )
   for run_slurm_line in run_slurm_lines: fid.write('sbatch %s\n' % run_slurm_line )
   fid.close()
   print( '\n\nTo kick off new jobs type:  source rerun_sbatch_commands.sh\n' )
else:
   print( '\nRerun with --clean option to actually clean out directories and prepare for re-run.' )
		


