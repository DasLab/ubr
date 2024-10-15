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

args = parser.add_argument( 'ubr_run_dirs',default=['./'],nargs='*',help='Directories in which ubr_run was run')
#args = parser.add_argument( '--fast', action='store_true',help='Only look through logs, no actual data analysis' )
args = parser.parse_args()


def check_stats( ubr_run_dir ):
    # Merge
    numreads_into_merge = numreads_merged = 0
    logfiles = glob.glob('%s/0_merge_pairs.err' % ubr_run_dir)
    if len( logfiles ) > 0:
        assert( len(logfiles) == 1 )
        logfile = logfiles[0]
        lines = open( logfile ).readlines()
        for line in lines:
            if line.find('Pairs:')>-1: numreads_into_merge = int(line.split()[1])
            if line.find('Joined:')>-1: numreads_merged = int(line.split()[1])

    # ultraplex
    numreads_into_demux = numreads_demultiplexed = 0
    logfiles = glob.glob('%s/1_ultraplex/ultraplex*log' % ubr_run_dir)
    if len( logfiles ) > 0:
        assert( len(logfiles) == 1 )
        logfile = logfiles[0]
        lines = open( logfile ).readlines()
        for line in lines:
            if line.find('Demultiplexing')>-1: numreads_into_demux = int(line.split()[5])
            if line.find('correctly')>-1: numreads_demultiplexed = int(line.split()[3])

    # bowtie2
    numreads_into_bt2 = numreads_align = 0
    logfiles = glob.glob('%s/2_bowtie2/*/bowtie2.err' % ubr_run_dir)
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
    logfiles = glob.glob('%s/3_rf_count/*out' % ubr_run_dir)
    if len( logfiles ) > 0:
        for logfile in logfiles:
            lines = open( logfile ).readlines()
            numreads = 0
            for line in lines:
                if line.find('covered')>-1:
                    cols = line.split()
                    if len(cols)>15:
                        cols = line.split()[15].split('/')
                        numreads_into_rfcount += int(cols[1])
                        numreads_counted += int(cols[0])

    def writeout(tag,numreads,numreads_prev = 0 ):
        if numreads>0:
            out_line = '%16s: %9d' % (tag,numreads)
            if numreads_prev > 0:
                out_line += '   (%5.1f%%)' % (100 * numreads/numreads_prev)
            print(out_line)


    print()
    print( 'Reads table for', os.path.abspath(ubr_run_dir) )
    writeout( 'Into bbmerge', numreads_into_merge )
    writeout( 'Merged pairs', numreads_merged, numreads_into_merge )
    writeout( 'Into ultraplex', numreads_into_demux, numreads_merged )
    writeout( 'Demultiplexed', numreads_demultiplexed, numreads_into_demux )
    writeout( 'Into bowtie2', numreads_into_bt2, numreads_demultiplexed )
    writeout( 'Aligned bowtie2', numreads_align, numreads_into_bt2 )
    writeout( 'Into rf-count', numreads_into_rfcount, numreads_align )
    writeout( 'Counted rf-count', numreads_counted, numreads_into_rfcount )
    print()



for ubr_run_dir in args.ubr_run_dirs:
    check_stats( ubr_run_dir )
