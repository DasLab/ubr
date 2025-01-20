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

    #minimap2 (alternative to bowtie2)
    logfiles = glob.glob('%s/2_bowtie2/*/minimap2.err' % ubr_run_dir)
    numreads_align_minimap2 = 0
    if len( logfiles ) > 0:
        for logfile in logfiles:
            lines = open( logfile ).readlines()
            numreads = 0
            for line in lines:
                if line.find('mapped')>-1: numreads_align_minimap2 += int(line.split()[2])

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

    # cmuts (alternative to rf-count)
    numreads_into_cmuts = numreads_counted_cmuts = 0
    logfiles = glob.glob('%s/3_cmuts/*log' % ubr_run_dir)
    if len( logfiles ) > 0:
        for logfile in logfiles:
            lines = open( logfile ).readlines()
            numreads = 0
            for line in lines:
                if line.find('Total')>-1 and line.find('reads')>-1:
                    cols = line.split()
                    numreads_into_cmuts += int(cols[-1].replace(',',''))
                if line.find('Mapped')>-1:
                    cols = line.split()
                    numreads_counted_cmuts += int(cols[-1].replace(',',''))

    def writeout(tag,numreads,numreads_prev, allreads):
        if numreads>0:
            out_line = '%16s: %9d' % (tag,numreads)
            if numreads_prev > 0:
                out_line += '   (%5.1f%%)' % (100 * numreads/numreads_prev)
            print(out_line)
            allreads.append(numreads)

    allreads = []
    print()
    print( 'Reads table for', os.path.abspath(ubr_run_dir) )
    writeout( 'Into bbmerge', numreads_into_merge, 0, allreads )
    writeout( 'Merged pairs', numreads_merged, numreads_into_merge, allreads )
    writeout( 'Into ultraplex', numreads_into_demux, numreads_merged, allreads )
    writeout( 'Demultiplexed', numreads_demultiplexed, numreads_into_demux, allreads )
    writeout( 'Into bowtie2', numreads_into_bt2, numreads_demultiplexed, allreads )
    writeout( 'Aligned bowtie2', numreads_align, numreads_into_bt2, allreads )
    writeout( 'Aligned minimap2', numreads_align_minimap2, numreads_demultiplexed, allreads )
    writeout( 'Into rf-count', numreads_into_rfcount, numreads_align, allreads )
    writeout( 'Counted rf-count', numreads_counted, numreads_into_rfcount, allreads )
    if numreads_align_minimap2>0:
        writeout( 'Into cmuts', numreads_into_cmuts, numreads_align_minimap2, allreads )
    else:
        writeout( 'Into cmuts', numreads_into_cmuts, numreads_align, allreads )
    writeout( 'Counted cmuts', numreads_counted_cmuts, numreads_into_cmuts, allreads )
    if len(allreads)>0: print( 'Overall: %9d/%9d = %6.2f%% ' % (min(allreads),max(allreads),100*min(allreads)/max(allreads)))
    print()


if __name__ == "__main__":
    args = parser.parse_args()
    for ubr_run_dir in args.ubr_run_dirs:
        check_stats( ubr_run_dir )
