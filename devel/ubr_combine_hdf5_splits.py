#!/usr/bin/env python3

import argparse
import os
import os.path
import glob
import shutil
import time
import gzip
import math
import h5py

parser = argparse.ArgumentParser(
                    prog = 'ubr_combine_hdf5_splits.py',
                    description = 'Merge hdf5 splits directories',
                    epilog = 'Merge files with names like RTB000.1_100000.hdf5, etc.' )

args = parser.add_argument( 'split_files',help='hdf5 files to merge',nargs='+' )
args = parser.add_argument( '--no_compression', action='store_true',help='Do not compress hdf5 files' )
args = parser.add_argument( '--no_overwrite', action='store_true',help='Do not merge if merged files already exist' )
args = parser.parse_args()

time_start = time.time()
time_readin = 0
time_output = 0

split_files = args.split_files
tags = []
start_seqs = []
end_seqs = []
ds_all = []
f_all = []
for (i,split_file) in enumerate(split_files):
    # name assumed to be RTB000_Marathon_Bicine_3pct_DMS.1_125000.hdf5
    # later could look inside hdf5 attributes?
    cols = split_file.split('.')
    assert( cols[-1] == 'hdf5' )
    tags.append( '.'.join(cols[:-2]) )
    start_seqs.append( int(cols[-2].split('_')[0]) )
    end_seqs.append( int(cols[-2].split('_')[1]) )

    try:
        f = h5py.File(split_file,'r')
        dataname = list(f.keys())[0]
        ds_all.append( f[dataname] )
        f_all.append(f)
    except:
        pass

assert( len(set(tags)) == 1 )
numfiles = len(ds_all)

tag = tags[0]
outfile = tag + '.hdf5'
f_out = h5py.File(outfile,'w')

ds = ds_all[0]
nseq = ds['mutations'].shape[0]
compression = 'gzip'
if args.no_compression: compression = None
for ds_type in ['mutations','insertions']:
    f_out.create_dataset( '%s/%s' % (tag,ds_type), ds[ds_type].shape,
                          dtype=ds[ds_type].dtype,chunks=ds[ds_type].chunks,compression=compression )
ds_out = f_out[tag]

tot_counts = 0
for (f,ds,split_file,start_seq,end_seq) in zip(f_all,ds_all,split_files,start_seqs,end_seqs):
    print( 'Reading in from: ',split_file )
    for ds_type in ['mutations','insertions']:
        time_start_read = time.time()
        s = range(start_seq-1,end_seq)
        chunk = ds[ds_type][s]
        time_end_read = time.time()

        if ds_type=='mutations': tot_counts += chunk.max(axis=1).sum()
        ds_out[ds_type][s] = chunk
        time_end_write = time.time()

        time_readin += (time_end_read - time_start_read)
        time_output += (time_end_write - time_end_read)


time_end = time.time()

print( 'Compiled %8d total counts for %6d sequences from %6d of %6d files into: %s' % (tot_counts,nseq,numfiles,len(split_files),outfile) )


print( '\nTimings:')
print( 'Read in data: ' + time.strftime("%H:%M:%S",time.gmtime(time_readin) ) )
print( 'Output  data: ' + time.strftime("%H:%M:%S",time.gmtime(time_output) ) )

print( '\nTotal time: ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_start) ) )

