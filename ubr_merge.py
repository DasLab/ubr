#!/usr/bin/env python3

import argparse
import os
import os.path
import glob
import shutil
import time
import gzip
import math

parser = argparse.ArgumentParser(
                    prog = 'ubr_merge.py',
                    description = 'Merge outputs from split ubr_run directories',
                    epilog = 'Merge .muts.txt and .coverage.txt output files from split ubr_run job directories.')

args = parser.add_argument( 'split_dir',default='UBR/',help='name of directory with 001/,002/ job subdirectories; if multiple directories, just look directly inside directories.',nargs='+' )
args = parser.add_argument( '--merge_files' ,nargs='+',help='Filenames to find and merge')
args = parser.add_argument( '-s','--setup_slurm', action='store_true',help='Instead of running, create sbatch script' )
args = parser.add_argument('-j','--jobs_per_slurm_node', default=24,type=int )
args = parser.add_argument('--start_seq', default=0,type=int,help='which sequence to start with [default 1]' )
args = parser.add_argument('--end_seq', default=0,type=int,help='which sequence to end on [default Nseq]' )
args = parser.add_argument('-n','--nsplits', default=0, type=int, help='number of separate partitions & threads for hdf5 split; creates SLURM job' )
args = parser.add_argument( '--no_compression', action='store_true',help='Do not compress hdf5 files' )
args = parser.add_argument( '--no_overwrite', action='store_true',help='Do not merge if merged files already exist' )
args = parser.add_argument( '--save_split_files', action='store_true',help='Save hdf5,log,err files from splits instead of removing' )
args = parser.add_argument('-O','--outdir',default='',help='output directory for files')
args = parser.parse_args()

split_dirs = args.split_dir
for (i,split_dir) in enumerate(split_dirs):
    if split_dir[:-1] != '/': split_dirs[i]=split_dir+'/'
if len(split_dirs) == 1: split_dirs[0] = split_dirs[0]+'*/'
print("Merging from the following directories: ")
for split_dir in split_dirs: print(' ',split_dir)
if len(args.outdir)>0 and args.outdir[-1] != '/': args.outdir += '/'

merge_files = args.merge_files

if merge_files == None:
    unique_files = []
    for split_dir in split_dirs:

        # cmuts output is .hdf5
        globfiles = sorted(glob.glob('%s/*.hdf5' % (split_dir)  ))
        # Get unique names
        unique_files += sorted(list(set(map( os.path.basename, globfiles ))))

        # ubr/RNAFramework output is .txt.gz
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

has_hdf5 = False
print('\nWill merge:')
for merge_file in merge_files:
    print(merge_file)
    if merge_file.find( '.hdf5')>-1: has_hdf5 = True

if args.setup_slurm and args.nsplits==0:
    if len(merge_files) < args.jobs_per_slurm_node:
        sbatch_file = 'run_ubr_merge.sh'
        fid_slurm = open( sbatch_file, 'w' )
        sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_merge\n#SBATCH --output=ubr_merge.o%%j\n#SBATCH --error=ubr_merge.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=8:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n#SBATCH --mem=%dG\n\n' % (len(merge_files),4*len(merge_files))
        fid_slurm.write( sbatch_preface )
        if has_hdf5: fid_slurm.write('ml py-h5py/3.10.0_py312\n')
        outdir_tag = ''
        if len(args.outdir)>0: outdir_tag = ' --outdir %s' % args.outdir
        for merge_file in merge_files:
            command = 'ubr_merge.py %s --merge_files %s %s&' % ( ' '.join(args.split_dir), merge_file, outdir_tag)
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
        if has_hdf5: fid_slurm.write('ml py-h5py/3.10.0_py312\n')
        fid_sbatch_commands = open( 'sbatch_merge_commands.sh', 'w')

        for (i,merge_file) in enumerate(merge_files):
            outdir_tag = ''
            if len(args.outdir)>0: outdir_tag = ' --outdir %s' % args.outdir
            command = 'ubr_merge.py %s --merge_files %s %s&' % ( ' '.join(args.split_dir), merge_file, outdir_tag)
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

# separately process chunks of hdf5 files in different processes
if args.nsplits > 0:
    slurm_file_dir = 'slurm_files'
    os.makedirs(slurm_file_dir, exist_ok=True )
    slurm_file_count = 0
    slurm_command_count = 0
    fid_sbatch_commands = open( 'sbatch_merge_commands.sh', 'w')
    slurm_combine_file_count = 0
    fid_sbatch_combine_commands = open( 'sbatch_combine_commands.sh', 'w')

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
        outdir = args.outdir

    if len(infiles) == 0:
        for split_dir in split_dirs:
            infiles += sorted(glob.glob('%s/raw_counts/%s' % (split_dir,filename) ))
            infiles += sorted(glob.glob('%s/raw_counts/%s.gz' % (split_dir,filename) ))
            outdir = args.outdir + 'raw_counts/'
            os.makedirs( outdir, exist_ok = True )

    assert(len(infiles)>0)
    if len(filename) > 4 and filename[-5:]=='.hdf5': # handle cmuts output (HDF5 format)
        import h5py
        assert( shutil.which('ubr_combine_hdf5_splits.py') ) # found in ubr/devel/

        ds_all = []
        f_all = []
        if args.nsplits > 0: # NEED TO UNIFY WITH SLURM CODE BLOCK ABOVE!
            nseq = 0
            for infile in infiles: # figure out sequence length
                try:
                    f = h5py.File(infile,'r')
                    dataname = list(f.keys())[0]
                    nseq = f[dataname]['mutations'].shape[0]
                except:
                    continue
                if nseq > 0: break

            num_slurm_files = math.ceil(args.nsplits/args.jobs_per_slurm_node)
            chunk_size = math.ceil(nseq/args.nsplits)
            out_tags = []
            for k in range(num_slurm_files):
                slurm_file_start_seq = k * args.jobs_per_slurm_node * chunk_size + 1
                slurm_file_end_seq   = min(slurm_file_start_seq + args.jobs_per_slurm_node * chunk_size - 1, nseq)
                num_jobs_for_slurm_file = math.ceil( (slurm_file_end_seq - slurm_file_start_seq + 1) / chunk_size )
                slurm_file_count += 1
                fid_slurm = open( '%s/run_ubr_merge_%03d.sh' % (slurm_file_dir, slurm_file_count), 'w' )
                sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_merge\n#SBATCH --output=ubr_merge.o%%j\n#SBATCH --error=ubr_merge.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=24:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n#SBATCH --mem=%dG\n\n' % (num_jobs_for_slurm_file,4*num_jobs_for_slurm_file)
                fid_slurm.write( sbatch_preface )
                fid_slurm.write('ml py-h5py/3.10.0_py312\n')
                fid_slurm.write( 'date\n' )
                fid_slurm.write( 'echo "Kicking off %d jobs"\n' % args.nsplits )
                tmp_dir = args.outdir + 'tmp_merge/'
                os.makedirs( tmp_dir, exist_ok = True )
                for n in range(num_jobs_for_slurm_file):
                    slurm_command_count += 1
                    start_seq = slurm_file_start_seq + n*chunk_size
                    end_seq   = min( start_seq + chunk_size - 1, nseq)
                    out_file = tmp_dir+filename
                    out_tag = out_file.replace('.hdf5','.%07d_%07d.hdf5' % (start_seq,end_seq))
                    out_tags.append( out_tag )
                    command = 'ubr_merge.py %s --merge_files %s --start_seq %s --end_seq %s --no_compression --outdir %s > %s.out 2> %s.err &' % \
                        ( ' '.join(args.split_dir), filename, start_seq,end_seq,tmp_dir,out_tag,out_tag )
                    fid_slurm.write( command +'\n' )
                fid_slurm.write('\nwait\n')
                fid_slurm.write('\necho "DONE"\n')
                fid_slurm.write( 'date\n' )
                fid_slurm.close()
                fid_sbatch_commands.write('sbatch %s\n' % fid_slurm.name )

            # also prepare a job to run_ubr_combine.
            slurm_combine_file_count += 1
            fid_slurm = open( '%s/run_ubr_combine_%03d.sh' % (slurm_file_dir, slurm_combine_file_count), 'w' )
            num_jobs_per_slurm_file = 1
            sbatch_preface = '#!/bin/bash\n#SBATCH --job-name=ubr_combine\n#SBATCH --output=ubr_combine.o%%j\n#SBATCH --error=ubr_combine.e%%j\n#SBATCH --partition=biochem,owners\n#SBATCH --time=24:00:00\n#SBATCH -n %d\n#SBATCH -N 1\n#SBATCH --mem=%dG\n\n' % (num_jobs_for_slurm_file,4*num_jobs_for_slurm_file)
            fid_slurm.write('ubr_combine_hdf5_splits.py %s%s\n' % (tmp_dir,filename.replace('.hdf5','.*_*.hdf5') ))
            fid_slurm.write( sbatch_preface )
            fid_slurm.write('ml py-h5py/3.10.0_py312\n')
            fid_slurm.write( 'date\n' )
            if not args.save_split_files:
                for out_tag in out_tags:
                    fid_slurm.write('rm  %s %s.out %s.err\n' % (out_tag,out_tag,out_tag))
            fid_slurm.write('\necho "DONE"\n')
            fid_slurm.write( 'date\n' )
            fid_slurm.close()
            fid_sbatch_combine_commands.write('sbatch %s\n' % fid_slurm.name )

            continue

        print('Reading in h5py files')
        for infile in infiles:
            try:
                f = h5py.File(infile,'r')
                dataname = list(f.keys())[0]
                ds_all.append( f[dataname] )
                f_all.append(f)
            except:
                continue
        numfiles = len(ds_all)
        ds = ds_all[0]
        nseq = ds['mutations'].shape[0]

        filename_out = filename
        if args.start_seq > 0: filename_out = filename_out.replace('.hdf5','.%07d_%07d.hdf5' % (args.start_seq,args.end_seq))
        outfile = outdir + filename_out
        f_out = h5py.File(outfile,'w')

        compression = 'gzip'
        if args.no_compression: compression = None
        for ds_type in ['mutations','insertions']:
            f_out.create_dataset( '%s/%s' % (filename_out,ds_type), ds[ds_type].shape,
                                  dtype=ds[ds_type].dtype,chunks=ds[ds_type].chunks,compression=compression )
        ds_out = f_out[filename_out]

        tot_counts = 0
        for ds_type in ['mutations','insertions']:
            print('Adding up %s by chunk from %d files' % (ds_type, numfiles))
            chunks = ds_out[ds_type].iter_chunks()
            if args.start_seq > 0: chunks = [range(args.start_seq-1,args.end_seq)]
            for s in chunks:
                time_start_read = time.time()

                chunk = ds_all[0][ds_type][s]
                for (q,ds) in enumerate(ds_all[1:]):
                    print( 'Adding chunk %d of %d',q,len(ds_all) )
                    chunk += ds[ds_type][s]
                    if (q % 10) == 0: ds_out[ds_type][s] = chunk # to allow tracking of progress...
                time_end_read = time.time()

                if ds_type=='mutations': tot_counts += chunk.max(axis=1).sum()
                ds_out[ds_type][s] = chunk
                time_end_write = time.time()

                time_readin += (time_end_read - time_start_read)
                time_output += (time_end_write - time_end_read)

        print( 'Compiled %8d total counts for %6d sequences from %6d of %6d files into: %s' % (tot_counts,nseq,numfiles,len(infiles),outfile) )

        for f in f_all: f.close()
        f_out.close()

    else: # handle .txt.gz files with pandas dataframes
        import pandas as pd
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

# Finish up separately process chunks of hdf5 files in different processes
if args.nsplits > 0:
    fid_sbatch_commands.close()
    print( '\nCreated %s slurm files containing %d commands. Run:\n source %s\n' % (slurm_file_count,slurm_command_count,fid_sbatch_commands.name) )
    fid_sbatch_commands.close()
    print( '\n[ Then run: source %s ]\n' % (fid_sbatch_combine_commands.name) )
    exit(0)

print( '\nTimings:')
print( 'Read in data: ' + time.strftime("%H:%M:%S",time.gmtime(time_readin) ) )
print( 'Output  data: ' + time.strftime("%H:%M:%S",time.gmtime(time_output) ) )

print( '\nTotal time: ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_start) ) )

