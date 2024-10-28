#!/usr/bin/env python3
import argparse
import h5py
import os
import time
import shutil
import gzip
import numpy as np
import string

parser = argparse.ArgumentParser(
                    prog = 'cmuts_to_ubr.py',
                    description = 'Convert cmuts hdf5 output to text files',
                    epilog = 'cmuts to hdf5' )

parser.add_argument('cmuts_file', help='HDF5-formatted file output by cmuts')
parser.add_argument('-o','--out_prefix', default='out', help='Name of output prefix')

args = parser.parse_args()
cmuts_file = args.cmuts_file
prefix = args.out_prefix

f = h5py.File(cmuts_file,'r')
dataname = list(f.keys())[0]
ds = f[dataname]
N = ds['mutations'].shape[0]

wd = './'


def writeline(fid,outline):
    for k in range(len(outline)):
       fid.write( '%d' % outline[k] )
       if k<(len(outline)-1): fid.write(',')
    fid.write('\n')
    #fid.write( np.array2string(outline) + '\n')

def writelines(outfile,outdata):
    fid = gzip.open(outfile,'wt')
    N = outdata.shape[0]
    for n in range( N ): writeline(fid,outdata[n])
    fid.close()
    print( 'Created: %s for %d sequences with coverage %d' % (outfile,N,outdata.max(axis=1).sum()) )


outfile_coverage = wd + '%s.coverage.txt.gz' % prefix
coverage = ds['mutations'][...].sum(2).sum(2)
writelines(outfile_coverage,coverage)

outdir = wd+'raw_counts/'
os.makedirs(outdir,exist_ok = True)

nts = 'ACGT'
mut_types = ['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del']
muts = np.zeros_like( coverage )
for mut_type in mut_types[:12]:
    outfile_raw_counts = outdir + '%s.%s.txt.gz' % (prefix,mut_type)
    i = nts.find(mut_type[0])
    j = nts.find(mut_type[1])
    writelines(outfile_raw_counts, ds['mutations'][:,:,i,j] )
    muts += ds['mutations'][:,:,i,j]

outfile_raw_counts = outdir + '%s.%s.txt.gz' % (prefix,'del')
dels = ds['mutations'][:,:,:,4].sum(axis=2)
writelines(outfile_raw_counts, dels )
muts += dels

outfile_raw_counts = outdir + '%s.%s.txt.gz' % (prefix,'ins')
writelines(outfile_raw_counts, ds['insertions'][:,:,:].sum(axis=2) )

outfile_muts = wd + '%s.muts.txt.gz' % prefix
writelines(outfile_muts,muts)



