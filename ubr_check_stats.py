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
                    prog = 'ubr_check_stats.py',
                    description = 'Check statistics of reads from an ubr_run.py run',
                    epilog = 'Look through logs of ubr_run.py and figure out reads and efficiency.')


args = parser.add_argument( '--fast', action='store_true',help='Only look through logs, no actual data analysis' )
args = parser.parse_args()


# Merge
numreads_into_merge = numreads_merged = 0
logfiles = glob.glob('0_merge_pairs.err')
if len( logfiles ) > 0:
    assert( len(logfiles) == 1 )
    logfile = logfiles[0]
    lines = open( logfile ).readlines()
    for line in lines:
        if line.find('Pairs:')>-1: numreads_into_merge = int(line.split()[1])
        if line.find('Joined:')>-1: numreads_merged = int(line.split()[1])

# ultraplex
numreads_into_demux = numreads_demultiplexed = 0
logfiles = glob.glob('1_ultraplex/ultraplex*log')
if len( logfiles ) > 0:
    assert( len(logfiles) == 1 )
    logfile = logfiles[0]
    lines = open( logfile ).readlines()
    for line in lines:
        if line.find('Demultiplexing')>-1: numreads_into_demux = int(line.split()[5])
        if line.find('correctly')>-1: numreads_demultiplexed = int(line.split()[3])

# bowtie2
numreads_into_bt2 = numreads_align = 0
logfiles = glob.glob('2_bowtie2/*/bowtie2.err')
if len( logfiles ) > 0:
    for logfile in logfiles:
        lines = open( logfile ).readlines()
        numreads = 0
        for line in lines:
            if line.find('reads;')>-1: numreads_into_bt2 += int(line.split()[0])
            if line.find('exactly')>-1: numreads_align += int(line.split()[0])
            if line.find('>1 times')>-1: numreads_align += int(line.split()[0])

# rf-count
numreads_into_rfcount = numreads_counted = 0
logfiles = glob.glob('3_rf_count/*out')
if len( logfiles ) > 0:
    for logfile in logfiles:
        lines = open( logfile ).readlines()
        numreads = 0
        for line in lines:
            if line.find('covered')>-1:
                cols = line.split()[15].split('/')
                numreads_into_rfcount += int(cols[1])
                numreads_counted += int(cols[0])

def writeout(tag,numreads):
    if numreads>0: print( '%16s: %9d' % (tag,numreads) )

print( os.getcwd() )
print()
print( 'Reads table' )
writeout( 'Into bbmerge', numreads_into_merge )
writeout( 'Merged pairs', numreads_merged )
writeout( 'Into ultraplex', numreads_into_demux )
writeout( 'Demultiplexed', numreads_demultiplexed )
writeout( 'Into bowtie2', numreads_into_bt2 )
writeout( 'Aligned bowtie2', numreads_align )
writeout( 'Into rf-count', numreads_into_rfcount )
writeout( 'Counted rf-count', numreads_counted )
print()

